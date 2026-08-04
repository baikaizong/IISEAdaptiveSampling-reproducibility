!==============================================================
!  Module: Experiment_Runner_mod
!  Purpose: Define POS experiment configurations and execution flow.
!==============================================================
module Experiment_Runner_mod
    use mpi
    use GlobalSettings_mod
    use POS_mod            
    use utils_mod
    use SamplesGeneration_mod
    use performance_evaluation 
    use Data_Loader 
    use RuntimeConfig_mod, only: cache_path, ensure_directory
    implicit none
    
    public :: Run_POS_Bspline_Test
    public :: Run_POS_Kernel_CircleConv_Test

contains

!==============================================================
!  Subroutine: Run_POS_Bspline_Test
!  Description: Standard test with B-Spline Background & Anomaly
!==============================================================
subroutine Run_POS_Bspline_Test(base_output_dir)
    implicit none
    
    ! Input
    character(len=*), intent(in) :: base_output_dir

    ! Variables
    integer :: ierr, ranks, num_procs, stat, i
    
    ! Simulation results
    real(dp) :: BLimitL, BLimitR, Limit, arl, std
    real(dp) :: max_arl, q1_arl, q2_arl, q3_arl
    real(dp) :: cur_amp
    
    ! Fixed Parameters
    integer :: noise_type, back_index, anomaly_type
    
    ! Spatial data
    real(dp), allocatable :: Nodesset(:,:)
    real(dp), allocatable :: NeighborDis(:,:), Kernel_NeighborDis(:,:)
    integer, allocatable :: NeighborId(:,:), Kernel_NeighborId(:,:)
    type(Bbasis_sparse_type) :: B0_all, B1_all, CovMat
    
    ! File handling
    character(len=512) :: output_dir
    character(len=512) :: file_limit, file_res_amp, file_config
    
    ! MPI Info
    call MPI_Comm_rank(MPI_COMM_WORLD, ranks, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)

    ! --- Construct Paths ---
    output_dir = trim(base_output_dir) // '/POS/Bspline'
    call ensure_directory(output_dir)
    
    ! Create Directory (Rank 0 only)
    if (ranks == 0) then
        print *, "=========================================================="
        print *, " STARTING POS EXPERIMENT: B-Spline Back + B-Spline Ano"
        print *, " Output Dir: ", trim(output_dir)
        print *, "=========================================================="
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    file_limit     = trim(output_dir) // '/Limit_BsplineBack_Gauss.txt'
    file_res_amp   = trim(output_dir) // '/OC_BsplineAno_VaryAmp.txt'
    file_config    = trim(output_dir) // '/Parameter_Settings.txt'

    ! --- 1. Settings ---
    call set_generalparams()
    
    ! Set POS Specific Parameters
    ! miu_min: min shift size, theta_1/theta_2: allocation control
    call set_POSparams(omega_POS_in=5, statisticwindowsize_POS_in=100, dim_POS_in=2, &
                       c1_POS_in=1.5d0, c2_POS_in=0.1d0, &
                       num_tkelnear_POS_in=500, num_ekelnear_POS_in=500)

    noise_type   = NOISE_GAUSSIAN
    back_index   = BACK_BSPLINE
    anomaly_type = TYPE_BSPLINE

    ! --- 2. Write Config Log ---
    if (ranks == 0) then
        open(unit=50, file=file_config, status='replace', action='write', iostat=stat)
        if (stat == 0) then
            write(50, '(A)') "########################################################"
            write(50, '(A)') "#           POS SIMULATION PARAMETER LOG             #"
            write(50, '(A)') "########################################################"
            write(50, *)
            write(50, '(A,I8)') " MPI Processes: ", num_procs
            write(50, '(A,I8)') " Simu Runs:     ", simu
            write(50, '(A,F12.4)') " Target ARL:    ", IcArl
            write(50, *)
            write(50, '(A)') "--- POS Parameters ---"
            call get_POSparams(50)
            write(50, *)
            write(50, '(A)') "--- Experiment Config ---"
            write(50, '(A,I4)') " Noise Type: ", noise_type
            write(50, '(A,I4)') " Back Index: ", back_index
            write(50, '(A,I4)') " Ano Type:   ", anomaly_type
            close(50)
        end if
    end if

    ! --- 3. Load Data & Init ---
    call Prepare_Simulation_Data(ranks, cache_path('Prepared_settings'), &
        Nodesset=Nodesset, B0_out=B0_all, B1_out=B1_all)
    call Prepare_POS_Kernels(ranks, cache_path('Prepared_settings'), Nodesset, &
                             NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId)

    ! POS needs both Kernel Neighbors (for CUSUM) and Standard Neighbors (for Sampling)
    call Init_Performance_Evaluation( &
        NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId, &
        Nodesset, B0_all, B1_all)

    ! --- 4. Limit Search ---
    call Set_SampleParams(noise_type_in=noise_type, back_index_in=back_index)
    call POS_Calibrate(noise_type, back_index)

    if (ranks == 0) print *, ">> Step 1: Estimating Limit Range (POS)..."
    Limit = 15.0d0 ! Initial Guess
    call POS_EstimateLimitRange(noise_type, back_index, Limit, BLimitL, BLimitR)
    
    if (ranks == 0) then
        print *, ">> Step 2: Exact Limit Search..."
        open(10, file=file_limit, status='replace')
    end if
    
    call POS_limitSearch(noise_type, back_index, &
                           BLimitL, BLimitR, Limit, arl, std, fid=10)
    
    if (ranks == 0) then
        write(10, *) "Final Limit:", Limit, " ARL:", arl
        close(10)
        print *, ">> Limit Found:", Limit
        
        open(unit=50, file=file_config, status='old', position='append', action='write')
        write(50, *)
        write(50, '(A,F12.6)') " Found Limit: ", Limit
        write(50, '(A,F12.6)') " IC ARL:      ", arl
        close(50)
    end if
    
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    ! --- 5. OCARL Experiment ---
    if (ranks == 0) then
        print *, ">> Step 3: Running OCARL (Varying Amplitude)..."
        open(20, file=file_res_amp, status='replace')
        write(20, '(A)') "Amp  ARL  Std  Q1  Med  Q3  Max"
    end if

    do i = 1, 20
        cur_amp = 0.5d0 * dble(i) 
        
        call POS_OCARL(noise_type, back_index, anomaly_type, &
                         ! Inputs
                         Limit=Limit, amplitude=cur_amp, &
                         ! Outputs
                         arl=arl, std=std, max_arl=max_arl, &
                         q1_arl=q1_arl, q2_arl=q2_arl, q3_arl=q3_arl)
        
        if (ranks == 0) then
            write(20, '(F6.2, 6(1X,F10.4))') cur_amp, arl, std, q1_arl, q2_arl, q3_arl, max_arl
            print '(A,F6.2,A,F8.2)', "  Amp=", cur_amp, " -> ARL=", arl
        end if
    end do
    
    if (ranks == 0) close(20)

    ! --- 6. Cleanup ---
    call Clean_Performance_Evaluation()
    if (allocated(Nodesset)) deallocate(Nodesset)
    if (allocated(NeighborDis)) deallocate(NeighborDis)
    if (allocated(NeighborId))  deallocate(NeighborId)
    if (allocated(Kernel_NeighborDis)) deallocate(Kernel_NeighborDis)
    if (allocated(Kernel_NeighborId))  deallocate(Kernel_NeighborId)
    call free_Bbasis_sparse(B0_all)
    call free_Bbasis_sparse(B1_all)
    call free_Bbasis_sparse(CovMat)

    if (ranks == 0) print *, "=== POS B-Spline Experiment Completed ==="

end subroutine Run_POS_Bspline_Test


!==============================================================
!  Subroutine: Run_POS_Kernel_CircleConv_Test
!  Description: Test with Kernel Random Background & Circle Conv Anomaly
!==============================================================
subroutine Run_POS_Kernel_CircleConv_Test(base_output_dir)
    implicit none
    
    ! Input
    character(len=*), intent(in) :: base_output_dir

    ! Variables
    integer :: ierr, ranks, num_procs, stat, i
    
    ! Results
    real(dp) :: BLimitL, BLimitR, Limit, arl, std
    real(dp) :: max_arl, q1_arl, q2_arl, q3_arl
    real(dp) :: cur_amp, cur_rad
    
    ! Params
    integer :: noise_type, back_index, anomaly_type
    real(dp) :: back_bw, ano_bw
    
    ! Data
    real(dp), allocatable :: Nodesset(:,:)
    real(dp), allocatable :: NeighborDis(:,:), Kernel_NeighborDis(:,:)
    integer, allocatable :: NeighborId(:,:), Kernel_NeighborId(:,:)
    type(Bbasis_sparse_type) :: B0_all, B1_all, CovMat
    
    ! Files
    character(len=512) :: output_dir
    character(len=512) :: file_limit, file_res_rad, file_res_amp, file_config
    
    call MPI_Comm_rank(MPI_COMM_WORLD, ranks, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)

    ! --- Construct Paths ---
    output_dir = trim(base_output_dir) // '/POS/Kernel_CircleConv'
    call ensure_directory(output_dir)
    
    if (ranks == 0) then
        print *, "=========================================================="
        print *, " STARTING POS EXPERIMENT: Kernel Back + CircleConv"
        print *, " Output Dir: ", trim(output_dir)
        print *, "=========================================================="
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    file_limit     = trim(output_dir) // '/Limit_KernelBack_Gauss.txt'
    file_res_rad   = trim(output_dir) // '/OC_VaryRadius_FixShift2.0.txt'
    file_res_amp   = trim(output_dir) // '/OC_VaryShift_FixRadius0.05.txt'
    file_config    = trim(output_dir) // '/Parameter_Settings.txt'

    ! --- 1. Settings ---
    call set_generalparams()
    call set_POSparams(omega_POS_in=5, statisticwindowsize_POS_in=100, dim_POS_in=2, &
                       c1_POS_in=1.5d0, c2_POS_in=0.1d0, &
                       num_tkelnear_POS_in=500, num_ekelnear_POS_in=500)

    noise_type   = NOISE_GAUSSIAN
    back_index   = BACK_KERNEL_RANDOM
    anomaly_type = TYPE_CIRCLE_CONV
    
    back_bw      = 1.0d0
    ano_bw       = 1.0d0

    ! --- 2. Write Config Log ---
    if (ranks == 0) then
        open(unit=50, file=file_config, status='replace', action='write', iostat=stat)
        if (stat == 0) then
            write(50, '(A)') "########################################################"
            write(50, '(A)') "#           POS KERNEL TEST PARAMETER LOG            #"
            write(50, '(A)') "########################################################"
            write(50, *)
            write(50, '(A,I8)') " MPI Processes: ", num_procs
            write(50, '(A,I8)') " Simu Runs:     ", simu
            write(50, '(A,F12.4)') " Target ARL:    ", IcArl
            write(50, *)
            write(50, '(A)') "--- POS Parameters ---"
            call get_POSparams(50)
            write(50, *)
            write(50, '(A)') "--- Experiment Config ---"
            write(50, '(A,I4)') " Noise Type: ", noise_type
            write(50, '(A,I4)') " Back Index: ", back_index
            write(50, '(A,I4)') " Ano Type:   ", anomaly_type
            write(50, '(A,F12.4)') " Back BW:    ", back_bw
            write(50, '(A,F12.4)') " Ano BW:     ", ano_bw
            close(50)
        end if
    end if

    ! --- 3. Load Data & Init ---
    call Prepare_Simulation_Data(ranks, cache_path('Prepared_settings'), &
        Nodesset=Nodesset, B0_out=B0_all, B1_out=B1_all)
    call Prepare_POS_Kernels(ranks, cache_path('Prepared_settings'), Nodesset, &
                             NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId)

    call Init_Performance_Evaluation( &
        NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId, &
        Nodesset, B0_all, B1_all, CovMat)

    ! --- 4. Limit Search ---
    ! Pass background bandwidth to SampleParams for correct generation
    call Set_SampleParams(noise_type_in=noise_type, &
                          back_index_in=back_index, &
                          back_bandwidth_in=back_bw)
    call POS_Calibrate(noise_type, back_index, back_bandwidth=back_bw)

    if (ranks == 0) print *, ">> Step 1: Estimating Limit Range (POS)..."
    Limit = 15.0d0 
    
    call POS_EstimateLimitRange(noise_type, back_index, &
                                    Limit, BLimitL, BLimitR, &
                                    back_bandwidth=back_bw)
    
    if (ranks == 0) then
        print *, ">> Step 2: Exact Limit Search..."
        open(10, file=file_limit, status='replace')
    end if
    
    call POS_limitSearch(noise_type, back_index, &
                           BLimitL, BLimitR, Limit, arl, std, fid=10, &
                           back_bandwidth=back_bw)
    
    if (ranks == 0) then
        write(10, *) "Final Limit:", Limit, " ARL:", arl
        close(10)
        print *, ">> Limit Found:", Limit
        
        open(unit=50, file=file_config, status='old', position='append', action='write')
        write(50, *)
        write(50, '(A,F12.6)') " Found Limit: ", Limit
        write(50, '(A,F12.6)') " IC ARL:      ", arl
        close(50)
    end if
    
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    ! --- 5. Experiment A: Varying Radius ---
    if (ranks == 0) then
        print *, ">> Step 3: Running OCARL (Varying Radius)..."
        open(20, file=file_res_rad, status='replace')
        write(20, '(A)') "Radius  ARL  Std  Q1  Med  Q3  Max"
    end if

    do i = 1, 20
        cur_rad = 0.01d0 * dble(i)
        cur_amp = 2.0d0

        call POS_OCARL(noise_type, back_index, anomaly_type, &
                         ! Inputs (Arrays + Params)
                         Limit=Limit, &
                         shiftradius=cur_rad, &
                         amplitude=cur_amp, &
                         anomaly_bandwidth=ano_bw, &
                         back_bandwidth=back_bw, &
                         ! Outputs
                         arl=arl, std=std, max_arl=max_arl, &
                         q1_arl=q1_arl, q2_arl=q2_arl, q3_arl=q3_arl)
        
        if (ranks == 0) then
            write(20, '(F6.3, 6(1X,F10.4))') cur_rad, arl, std, q1_arl, q2_arl, q3_arl, max_arl
            print '(A,F6.3,A,F8.2)', "  Rad=", cur_rad, " -> ARL=", arl
        end if
    end do
    
    if (ranks == 0) close(20)

    ! --- 6. Experiment B: Varying Amplitude ---
    if (ranks == 0) then
        print *, ">> Step 4: Running OCARL (Varying Shift)..."
        open(30, file=file_res_amp, status='replace')
        write(30, '(A)') "Shift  ARL  Std  Q1  Med  Q3  Max"
    end if

    do i = 1, 20
        cur_rad = 0.05d0
        cur_amp = 0.2d0 * dble(i) 

        call POS_OCARL(noise_type, back_index, anomaly_type, &
                         Limit=Limit, &
                         shiftradius=cur_rad, &
                         amplitude=cur_amp, &
                         anomaly_bandwidth=ano_bw, &
                         back_bandwidth=back_bw, &
                         arl=arl, std=std, max_arl=max_arl, &
                         q1_arl=q1_arl, q2_arl=q2_arl, q3_arl=q3_arl)
        
        if (ranks == 0) then
            write(30, '(F6.2, 6(1X,F10.4))') cur_amp, arl, std, q1_arl, q2_arl, q3_arl, max_arl
            print '(A,F6.2,A,F8.2)', "  Shift=", cur_amp, " -> ARL=", arl
        end if
    end do

    if (ranks == 0) close(30)

    ! --- 7. Cleanup ---
    call Clean_Performance_Evaluation()
    
    if (allocated(Nodesset)) deallocate(Nodesset)
    if (allocated(NeighborDis)) deallocate(NeighborDis)
    if (allocated(NeighborId))  deallocate(NeighborId)
    if (allocated(Kernel_NeighborDis)) deallocate(Kernel_NeighborDis)
    if (allocated(Kernel_NeighborId))  deallocate(Kernel_NeighborId)
    
    call free_Bbasis_sparse(B0_all)
    call free_Bbasis_sparse(B1_all)
    call free_Bbasis_sparse(CovMat)

    if (ranks == 0) print *, "=== POS Kernel Test Completed ==="

end subroutine Run_POS_Kernel_CircleConv_Test

subroutine Prepare_POS_Kernels(ranks, output_dir, Nodesset, &
                               Kernelhe_Dis, Kernelhe_Id, Kernelht_Dis, Kernelht_Id)
    integer, intent(in) :: ranks
    character(len=*), intent(in) :: output_dir
    real(dp), intent(in) :: Nodesset(:,:)
    real(dp), allocatable, intent(out) :: Kernelhe_Dis(:,:), Kernelht_Dis(:,:)
    integer, allocatable, intent(out) :: Kernelhe_Id(:,:), Kernelht_Id(:,:)
    character(len=1024) :: he_dis_file, he_id_file, ht_dis_file, ht_id_file
    logical :: files_exist(4)
    integer :: ierr

    allocate(Kernelhe_Dis(num_allnodes, num_ekelnear_POS))
    allocate(Kernelhe_Id(num_allnodes, num_ekelnear_POS))
    allocate(Kernelht_Dis(num_allnodes, num_tkelnear_POS))
    allocate(Kernelht_Id(num_allnodes, num_tkelnear_POS))

    he_dis_file = trim(output_dir) // '/POS_Kernelhe_Dis.dat'
    he_id_file  = trim(output_dir) // '/POS_Kernelhe_Id.dat'
    ht_dis_file = trim(output_dir) // '/POS_Kernelht_Dis.dat'
    ht_id_file  = trim(output_dir) // '/POS_Kernelht_Id.dat'

    if (ranks == 0) then
        inquire(file=he_dis_file, exist=files_exist(1))
        inquire(file=he_id_file, exist=files_exist(2))
        inquire(file=ht_dis_file, exist=files_exist(3))
        inquire(file=ht_id_file, exist=files_exist(4))
        if (all(files_exist)) then
            call load_dense_array(he_dis_file, Kernelhe_Dis)
            call load_dense_array_int(he_id_file, Kernelhe_Id)
            call load_dense_array(ht_dis_file, Kernelht_Dis)
            call load_dense_array_int(ht_id_file, Kernelht_Id)
        else
            call EpaKernelTopK(Nodesset, he_POS, num_ekelnear_POS, Kernelhe_Dis, Kernelhe_Id)
            call EpaKernelTopK(Nodesset, ht_POS, num_tkelnear_POS, Kernelht_Dis, Kernelht_Id)
            call save_dense_array(he_dis_file, Kernelhe_Dis)
            call save_dense_array_int(he_id_file, Kernelhe_Id)
            call save_dense_array(ht_dis_file, Kernelht_Dis)
            call save_dense_array_int(ht_id_file, Kernelht_Id)
        end if
    end if

    call MPI_Bcast(Kernelhe_Dis, size(Kernelhe_Dis), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(Kernelhe_Id, size(Kernelhe_Id), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(Kernelht_Dis, size(Kernelht_Dis), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(Kernelht_Id, size(Kernelht_Id), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
end subroutine Prepare_POS_Kernels

end module Experiment_Runner_mod
