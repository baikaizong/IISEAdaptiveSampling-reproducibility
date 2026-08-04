!==============================================================
!  Module: Experiment_Runner_mod
!  Purpose: Define mMKSTD experiment configurations and execution flow.
!==============================================================
module Experiment_Runner_mod
    use ExperimentSupport_mod, only: Init_Experiment_Common, Print_Banner, &
        Write_OCARL_Header, Write_OCARL_Row
    use mpi
    use GlobalSettings_mod
    use mMSTD_mod
    use utils_mod
    use SamplesGeneration_mod
    use performance_evaluation 
    use Data_Loader 
    use RuntimeConfig_mod, only: cache_path, ensure_directory
    USE ANORIN_INT
    implicit none
    
    public

    private :: Calc_mMSTD_Params, Execute_Limit_Search_Workflow, &
               Write_OCARL_Header, Write_OCARL_Row, Init_Experiment_Common

contains

!==============================================================
!  (1) B-Spline Test 
!==============================================================
subroutine Run_mMSTD_Bspline_Test(base_output_dir)
    implicit none
    character(len=*), intent(in) :: base_output_dir
    
    ! Locals
    integer :: ranks, num_procs, k
    real(dp) :: Limit, cur_amp
    type(ARL_Stats_type) :: res_OC, res_Exp, res_Exploit
    
    ! Config Params
    integer :: noise_type, back_index, anomaly_type
    real(dp) :: tem_K, P_val, noise_var
    character(len=1024) :: output_dir, file_res, str_K, str_P

    ! Spatial data
    real(dp), allocatable :: Nodesset(:,:), NeighborDis(:,:), Kernel_NeighborDis(:,:)
    integer, allocatable :: NeighborId(:,:), Kernel_NeighborId(:,:)
    type(Bbasis_sparse_type) :: B0_all, B1_all
    real(dp), parameter :: shift_size(2)    = (/8.0d0,10.0d0 /)
    ! --- 1. Init ---
    call Init_Experiment_Common(ranks, num_procs)
    output_dir = trim(base_output_dir) // '/mMKSTD/BsplineNoBack-2'
    call ensure_directory(output_dir)
    if (ranks == 0) call Print_Banner("B-Spline Test", output_dir)

    ! --- 2. Setup Parameters ---
    call set_generalparams(num_samplingnodes_in=30)
    noise_type = NOISE_GAUSSIAN; back_index = BACK_NONE; anomaly_type = TYPE_BSPLINE
    tem_K = 5.0d0; P_val = 0.90d0
    write(str_K, '(F4.1)') tem_K; str_K = trim(adjustl(str_K))
    write(str_P, '(F4.2)') P_val; str_P = trim(adjustl(str_P))
    noise_var=SQRT(sigma_noise**2.0+sigma_ground**2.0)
    call Calc_mMSTD_Params(K_in=tem_K, P_in=P_val, nodes_in=num_samplingnodes, &
                           Epah_in=0.02d0, lambda_in=0.90d0, sigma_in=0.1d0, noise_var=noise_var, &
                           rank=ranks, log_dir=output_dir, tag_K=str_K, tag_P=str_P)
    ! --- 3. Load Data ---
    call Prepare_Simulation_Data(ranks, cache_path('Prepared_settings-1'), &
            Nodesset=Nodesset, NeighborDis=NeighborDis, NeighborId=NeighborId, &
            Kernel_NeighborDis=Kernel_NeighborDis, Kernel_NeighborId=Kernel_NeighborId, &
            B0_out=B0_all, B1_out=B1_all)
    call Init_Performance_Evaluation(NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId, Nodesset, B0_all, B1_all)

    ! --- 4. In-control charting limit ---
    Limit = 4.17287964774770d0
    ! --- 5. OCARL Experiment ---
    file_res = trim(output_dir) // '/OC_Result-1-add.txt'
    if (ranks == 0) call Write_OCARL_Header(file_res, "Amp")

    do k = 1, SIZE(shift_size)
        cur_amp = shift_size(k) 
        call mMSTD_OCARL(noise_type, back_index, anomaly_type, Limit=Limit, amplitude=cur_amp, &
                         res_OC_ARL=res_OC, res_Exploration_ARL=res_Exp, res_Exploitation_Num=res_Exploit)
        
        if (ranks == 0) call Write_OCARL_Row(file_res, (/ cur_amp /), res_OC, res_Exp, res_Exploit)
    end do

    ! --- 6. Cleanup ---
    call Clean_Performance_Evaluation()
    call free_Bbasis_sparse(B0_all); call free_Bbasis_sparse(B1_all)
    if (allocated(Nodesset)) deallocate(Nodesset)
    if (allocated(NeighborDis)) deallocate(NeighborDis)
    if (allocated(NeighborId)) deallocate(NeighborId)
    if (allocated(Kernel_NeighborDis)) deallocate(Kernel_NeighborDis)
    if (allocated(Kernel_NeighborId)) deallocate(Kernel_NeighborId)
    
    if (ranks == 0) print *, "=== B-Spline Completed ==="

end subroutine Run_mMSTD_Bspline_Test


!==============================================================
!  (2) Kernel CircleConv Test
!==============================================================
subroutine Run_mMSTD_Kernel_CircleConv_Test(base_output_dir)
    implicit none
    character(len=*), intent(in) :: base_output_dir

    integer :: ierr, ranks, num_procs, i, k, s, shift_pointnum
    real(dp) :: Limit, arl, cur_amp
    type(ARL_Stats_type) :: res_OC, res_Exp, res_Exploit
    
    integer :: noise_type, back_index, anomaly_type
    real(dp) :: tem_K, P_val, noise_var
    character(len=1024) :: output_dir, file_limit, file_res_rad, file_res_amp, file_ICp, str_K, str_P
    
    real(dp), allocatable :: Nodesset(:,:), NeighborDis(:,:), Kernel_NeighborDis(:,:)
    integer, allocatable :: NeighborId(:,:), Kernel_NeighborId(:,:)
    real(dp):: falserate
    
    real(dp), parameter :: shift_size(6)    = (/ 5.0d0, 10.0d0, 20.0d0,30.0d0, 50.0d0, 100.0d0 /)
    integer, parameter :: shift_num(3) = (/ 1, 2, 3 /)
    
    call Init_Experiment_Common(ranks, num_procs)
    output_dir = trim(base_output_dir) // '/mMKSTD/Kernel_CircleConv'
    call ensure_directory(output_dir)
    if (ranks == 0) call Print_Banner("Kernel + CircleConv", output_dir)
    
    call set_generalparams(num_x_in=20, num_y_in=20, num_disnearspatial_in=50, num_samplingnodes_in=10)
    noise_type = NOISE_GAUSSIAN; back_index = BACK_KERNEL_RANDOM; anomaly_type = TYPE_RANDOM_POINTS_CONV
    tem_K = 1.0d0; P_val = 0.95d0
    write(str_K, '(F4.1)') tem_K; str_K = trim(adjustl(str_K))
    write(str_P, '(F4.2)') P_val; str_P = trim(adjustl(str_P))    
    call compute_conv_stddev( noise_var)
    call Calc_mMSTD_Params(K_in=tem_K, P_in=P_val, nodes_in=num_samplingnodes, &
                           Epah_in=0.1d0, lambda_in=0.90d0, sigma_in=0.1d0, noise_var=noise_var, &
                           rank=ranks, log_dir=output_dir, tag_K=str_K, tag_P=str_P)

    call Prepare_Simulation_Data(ranks, cache_path('Prepared_settings-2'), &
        Nodesset, NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId)
    call Init_Performance_Evaluation(NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId, Nodesset)

    ! --- Limit Search ---
    file_limit = trim(output_dir) // '/Limit_KernelBack_Gauss-11.txt'
    ! Using auto-range (.true.) for limit search
    call Execute_Limit_Search_Workflow(ranks, noise_type, back_index, &
                                       20.0d0, 30.0d0, Limit, arl, file_limit, .false.)
    ! --- Exp A: Vary Radius ---
    file_res_rad = trim(output_dir) // '/OC_Arl-11.txt'
    if (ranks == 0) call Write_OCARL_Header(file_res_rad, "Amp num_point")

    do s =  SIZE(shift_num), 1, -1
        shift_pointnum = shift_num(s)
        do k = SIZE(shift_size),1,-1
            cur_amp = shift_size(k)
            call mMSTD_OCARL(noise_type, back_index, anomaly_type, Limit=Limit, shiftvalue=cur_amp, shift_pointnum=shift_pointnum,&
                                res_OC_ARL=res_OC, res_Exploration_ARL=res_Exp, res_Exploitation_Num=res_Exploit)
            if (ranks == 0) call Write_OCARL_Row(file_res_rad, (/ cur_amp, 1.0d0*shift_pointnum /), res_OC, res_Exp, res_Exploit)
        end do
    end do

    call Clean_Performance_Evaluation()
    if (allocated(Nodesset)) deallocate(Nodesset)
    if (allocated(NeighborDis)) deallocate(NeighborDis)
    if (allocated(NeighborId)) deallocate(NeighborId)
    if (allocated(Kernel_NeighborDis)) deallocate(Kernel_NeighborDis)
    if (allocated(Kernel_NeighborId)) deallocate(Kernel_NeighborId)
    if (ranks == 0) print *, "=== Kernel Test Completed ==="

end subroutine Run_mMSTD_Kernel_CircleConv_Test

!==============================================================
!  (3) Circle Test 
!==============================================================
subroutine Run_mMSTD_Circle_comparison(base_output_dir)
    implicit none
    character(len=*), intent(in) :: base_output_dir
    
    integer :: ranks, num_procs, d, i, k, s, ierr
    real(dp) :: Limit, arl, cur_amp, cur_rad
    type(ARL_Stats_type) :: res_OC, res_Exp, res_Exploit
    
    integer :: noise_type, back_index, anomaly_type
    real(dp) :: tem_K, P_val, noise_var
    character(len=1024) :: output_dir, file_limit, file_res, file_ICp
    character(len=32) :: str_K, str_P
    real(dp):: falserate

    real(dp), parameter :: shift_size(3)    = (/ 0.5d0, 1.0d0, 3.0d0 /)
    real(dp), parameter :: shift_radius(2) = (/ 0.1d0, 0.3d0 /)
    
    real(dp), allocatable :: Nodesset(:,:), NeighborDis(:,:), Kernel_NeighborDis(:,:)
    integer, allocatable :: NeighborId(:,:), Kernel_NeighborId(:,:)

    call Init_Experiment_Common(ranks, num_procs)
    output_dir = trim(base_output_dir) // '/mMKSTD/ExactAlgorithm'
    call ensure_directory(output_dir)
    if (ranks == 0) call Print_Banner("Gauss + Circle", output_dir)

    call set_generalparams(num_x_in=10, num_y_in=10, num_disnearspatial_in=50, num_samplingnodes_in=5)
    noise_type = NOISE_GAUSSIAN; back_index = BACK_NONE; anomaly_type = TYPE_CIRCLE
    
    call Prepare_Simulation_Data(ranks, cache_path('Prepared_settings-4'), &
        Nodesset=Nodesset, NeighborDis=NeighborDis, NeighborId=NeighborId, &
        Kernel_NeighborDis=Kernel_NeighborDis, Kernel_NeighborId=Kernel_NeighborId)
    call Init_Performance_Evaluation(NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId, Nodesset)

    ! --- Outer Loops ---
        P_val = 0.90
        write(str_P, '(F4.2)') P_val; str_P = trim(adjustl(str_P))
        tem_K = 5.0
        write(str_K, '(F4.1)') tem_K; str_K = trim(adjustl(str_K))
            
        call Calc_mMSTD_Params(K_in=tem_K, P_in=P_val, nodes_in=num_samplingnodes, &
                               Epah_in=0.2d0, lambda_in=0.90d0, sigma_in=0.1d0, noise_var=sigma_noise, &
                               rank=ranks, log_dir=output_dir, tag_K=str_K, tag_P=str_P)

        file_limit = trim(output_dir) // '/Limit_K-'//trim(str_K)//'_P-'//trim(str_P)//'-4'//'.txt'
        call Execute_Limit_Search_Workflow(ranks, noise_type, back_index, &
                                            3.55d0, 3.70d0, Limit, arl, file_limit, .false.)
        file_res = trim(output_dir) // '/OC_K-'//trim(str_K)//'_P-'//trim(str_P)//'-4'//'.txt'
        if (ranks == 0) call Write_OCARL_Header(file_res, "Amp Radius")
        do s = 1,  SIZE(shift_radius)
            cur_rad = shift_radius(s)
            do k = 1, SIZE(shift_size)
                cur_amp = shift_size(k)
                call mMSTD_OCARL(noise_type, back_index, anomaly_type, Limit=Limit, shiftradius=cur_rad, shiftvalue=cur_amp, &
                                    res_OC_ARL=res_OC, res_Exploration_ARL=res_Exp, res_Exploitation_Num=res_Exploit)
                if (ranks == 0) call Write_OCARL_Row(file_res, (/ cur_amp, cur_rad /), res_OC, res_Exp, res_Exploit)
            end do
        end do
        call MPI_Barrier(MPI_COMM_WORLD, ierr) 

    
    call Clean_Performance_Evaluation()
    if (allocated(Nodesset)) deallocate(Nodesset)
    if (allocated(NeighborDis)) deallocate(NeighborDis)
    if (allocated(NeighborId)) deallocate(NeighborId)
    if (allocated(Kernel_NeighborDis)) deallocate(Kernel_NeighborDis)
    if (allocated(Kernel_NeighborId)) deallocate(Kernel_NeighborId)

    if (ranks == 0) print *, "=== Circle Test Completed ==="

end subroutine Run_mMSTD_Circle_comparison

!==============================================================
!  (3) Circle Test 
!==============================================================
subroutine Run_mMSTD_Circle_Test(base_output_dir)
    implicit none
    character(len=*), intent(in) :: base_output_dir
    
    integer :: ranks, num_procs, d, i, k, s, ierr
    real(dp) :: Limit, arl, cur_amp, cur_rad
    type(ARL_Stats_type) :: res_OC, res_Exp, res_Exploit
    
    integer :: noise_type, back_index, anomaly_type
    real(dp) :: tem_K, P_val
    character(len=1024) :: output_dir, file_limit, file_res, file_ICp
    character(len=32) :: str_K, str_P
    real(dp):: falserate

    real(dp), parameter :: Ksize(5)    = (/  5.0d0,  10.0d0,  15.0d0 , 20.0d0,  25.0d0 /)
    real(dp), parameter :: P_vector(2) = (/ 0.85d0, 0.95d0 /)
        real(dp), parameter :: shift_size(6)    = (/ 0.5d0, 1.0d0, 1.5d0, 2.0d0, 2.5d0, 3.0d0 /)
    real(dp), parameter :: shift_radius(4) = (/ 0.02d0, 0.03d0,0.04d0, 0.05d0 /)
    real(dp), allocatable :: Nodesset(:,:), NeighborDis(:,:), Kernel_NeighborDis(:,:)
    integer, allocatable :: NeighborId(:,:), Kernel_NeighborId(:,:)

    call Init_Experiment_Common(ranks, num_procs)
    output_dir = trim(base_output_dir) // '/mMKSTD/Circle-17'
    call ensure_directory(output_dir)
    if (ranks == 0) call Print_Banner("Gauss + Circle", output_dir)

    call set_generalparams(num_samplingnodes_in=50)
    noise_type = NOISE_GAUSSIAN; back_index = BACK_NONE; anomaly_type = TYPE_CIRCLE
    
    call Prepare_Simulation_Data(ranks, cache_path('Prepared_settings-1'), &
        Nodesset=Nodesset, NeighborDis=NeighborDis, NeighborId=NeighborId, &
        Kernel_NeighborDis=Kernel_NeighborDis, Kernel_NeighborId=Kernel_NeighborId)
    call Init_Performance_Evaluation(NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId, Nodesset)

    ! --- Outer Loops ---
    do d = 1, 1
        P_val = P_vector(d)
        write(str_P, '(F4.2)') P_val; str_P = trim(adjustl(str_P))
        
        do i = 3, size(Ksize)
            tem_K = Ksize(i)
            write(str_K, '(F4.1)') tem_K; str_K = trim(adjustl(str_K))
            
        call Calc_mMSTD_Params(K_in=tem_K, P_in=P_val, nodes_in=num_samplingnodes, &
                               Epah_in=0.02d0, lambda_in=0.90d0, sigma_in=0.1d0, noise_var=sigma_noise, &
                               rank=ranks, log_dir=output_dir, tag_K=str_K, tag_P=str_P)
            
            file_limit = trim(output_dir) // '/Limit_K-'//trim(str_K)//'_P-'//trim(str_P)//'.txt'
            call Execute_Limit_Search_Workflow(ranks, noise_type, back_index, &
                                               0.0d0, 0.0d0, Limit, arl, file_limit, .true.)
            call mMSTD_ICPerformance(noise_type, back_index, Limit, falserate)
            file_ICp = trim(output_dir)  // '/ICperformance_K-'//trim(str_K)//'_P-'//trim(str_P)//'.txt'
            if (ranks == 0) then
                open(unit=70, file=file_ICp, status='unknown', position='append', action='write')
                write(70, '(F8.6, 1X)', advance='no') falserate
                close(70)
            end if 
            call MPI_Barrier(MPI_COMM_WORLD, ierr)

            file_res = trim(output_dir) // '/OC_K-'//trim(str_K)//'_P-'//trim(str_P)//'.txt'
            if (ranks == 0) call Write_OCARL_Header(file_res, "Amp Radius")

            do s = 1,  SIZE(shift_radius)
                cur_rad = shift_radius(s)
                do k = 1, SIZE(shift_size)
                    cur_amp = shift_size(k)
                    call mMSTD_OCARL(noise_type, back_index, anomaly_type, Limit=Limit, shiftradius=cur_rad, shiftvalue=cur_amp, &
                                     res_OC_ARL=res_OC, res_Exploration_ARL=res_Exp, res_Exploitation_Num=res_Exploit)
                    if (ranks == 0) call Write_OCARL_Row(file_res, (/ cur_amp, cur_rad /), res_OC, res_Exp, res_Exploit)
                end do
            end do
            call MPI_Barrier(MPI_COMM_WORLD, ierr) 
        end do
    end do
    
    call Clean_Performance_Evaluation()
    if (allocated(Nodesset)) deallocate(Nodesset)
    if (allocated(NeighborDis)) deallocate(NeighborDis)
    if (allocated(NeighborId)) deallocate(NeighborId)
    if (allocated(Kernel_NeighborDis)) deallocate(Kernel_NeighborDis)
    if (allocated(Kernel_NeighborId)) deallocate(Kernel_NeighborId)

    if (ranks == 0) print *, "=== Circle Test Completed ==="

end subroutine Run_mMSTD_Circle_Test

subroutine Run_mMSTD_Spatial_comparison(base_output_dir)
    implicit none
    character(len=*), intent(in) :: base_output_dir
    
    integer :: ranks, num_procs, d, i, k, s, ierr
    real(dp) :: Limit, arl, cur_amp, cur_rad
    type(ARL_Stats_type) :: res_OC, res_Exp, res_Exploit
    
    integer :: noise_type, back_index, anomaly_type
    real(dp) :: tem_K, P_val
    character(len=1024) :: output_dir, file_limit, file_res, file_ICp
    character(len=32) :: str_K, str_P
    real(dp):: falserate

    real(dp), parameter :: shift_size(6)    = (/ 0.05d0, 0.10d0, 0.15d0, 0.20d0, 0.25d0, 0.30d0 /)
    
    real(dp), allocatable :: Nodesset(:,:), NeighborDis(:,:), Kernel_NeighborDis(:,:)
    integer, allocatable :: NeighborId(:,:), Kernel_NeighborId(:,:)

    call Init_Experiment_Common(ranks, num_procs)
    output_dir = trim(base_output_dir) // '/mMKSTD/Comparison_STnormal'
    call ensure_directory(output_dir)
    if (ranks == 0) call Print_Banner("Gauss + ST", output_dir)

    call set_generalparams(num_samplingnodes_in=50)
    noise_type = NOISE_GAUSSIAN; back_index = BACK_NONE; anomaly_type = TYPE_SPATIAL
    
    call Prepare_Simulation_Data(ranks, cache_path('Prepared_settings-1'), &
        Nodesset=Nodesset, NeighborDis=NeighborDis, NeighborId=NeighborId, &
        Kernel_NeighborDis=Kernel_NeighborDis, Kernel_NeighborId=Kernel_NeighborId)
    call Init_Performance_Evaluation(NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId, Nodesset)

    ! --- Outer Loops ---
        P_val = 0.90
        write(str_P, '(F4.2)') P_val; str_P = trim(adjustl(str_P))
        tem_K = 20.0
        write(str_K, '(F4.1)') tem_K; str_K = trim(adjustl(str_K))
            
        call Calc_mMSTD_Params(K_in=tem_K, P_in=P_val, nodes_in=num_samplingnodes, &
                               Epah_in=0.02d0, lambda_in=0.90d0, sigma_in=0.1d0, noise_var=sigma_noise, &
                               rank=ranks, log_dir=output_dir, tag_K=str_K, tag_P=str_P)
            
        Limit = 4.31944979939881d0
        file_res = trim(output_dir) // '/OC_K-'//trim(str_K)//'_P-'//trim(str_P)//'.txt'
        if (ranks == 0) call Write_OCARL_Header(file_res, "Amp")
        do k = 1, SIZE(shift_size)
            cur_amp = shift_size(k)
            call mMSTD_OCARL(noise_type, back_index, anomaly_type, Limit=Limit,  theta=cur_amp, &
                                res_OC_ARL=res_OC, res_Exploration_ARL=res_Exp, res_Exploitation_Num=res_Exploit)
            if (ranks == 0) call Write_OCARL_Row(file_res, (/ cur_amp /), res_OC, res_Exp, res_Exploit)
        end do
        call MPI_Barrier(MPI_COMM_WORLD, ierr) 

    
    call Clean_Performance_Evaluation()
    if (allocated(Nodesset)) deallocate(Nodesset)
    if (allocated(NeighborDis)) deallocate(NeighborDis)
    if (allocated(NeighborId)) deallocate(NeighborId)
    if (allocated(Kernel_NeighborDis)) deallocate(Kernel_NeighborDis)
    if (allocated(Kernel_NeighborId)) deallocate(Kernel_NeighborId)

    if (ranks == 0) print *, "=== Circle Test Completed ==="

end subroutine Run_mMSTD_Spatial_comparison

subroutine Run_mMSTD_ST_comparison(base_output_dir)
    implicit none
    character(len=*), intent(in) :: base_output_dir
    
    integer :: ranks, num_procs, d, i, k, s, ierr
    real(dp) :: Limit, arl, cur_amp, cur_rad
    type(ARL_Stats_type) :: res_OC, res_Exp, res_Exploit
    
    integer :: noise_type, back_index, anomaly_type
    real(dp) :: tem_K, P_val
    character(len=1024) :: output_dir, file_limit, file_res, file_ICp
    character(len=32) :: str_K, str_P
    real(dp):: falserate

    real(dp), parameter :: shift_size(3)    = (/0.1d0,0.01d0,0.001d0 /)
    real(dp), parameter :: shift_radius(3) = (/ 0.05d0, 0.005d0, 0.0005d0 /)
    
    real(dp), allocatable :: Nodesset(:,:), NeighborDis(:,:), Kernel_NeighborDis(:,:)
    integer, allocatable :: NeighborId(:,:), Kernel_NeighborId(:,:)

    call Init_Experiment_Common(ranks, num_procs)
    output_dir = trim(base_output_dir) // '/mMKSTD/Comparison_STEXP'
    call ensure_directory(output_dir)
    if (ranks == 0) call Print_Banner("EXP + ST", output_dir)

    call set_generalparams(num_samplingnodes_in=50)
    noise_type = NOISE_EXPONENTIAL; back_index = BACK_NONE; anomaly_type = TYPE_ST
    
    call Prepare_Simulation_Data(ranks, cache_path('Prepared_settings-1'), &
        Nodesset=Nodesset, NeighborDis=NeighborDis, NeighborId=NeighborId, &
        Kernel_NeighborDis=Kernel_NeighborDis, Kernel_NeighborId=Kernel_NeighborId)
    call Init_Performance_Evaluation(NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId, Nodesset)

    ! --- Outer Loops ---
        P_val = 0.90
        write(str_P, '(F4.2)') P_val; str_P = trim(adjustl(str_P))
        tem_K = 20.0
        write(str_K, '(F4.1)') tem_K; str_K = trim(adjustl(str_K))
            
        call Calc_mMSTD_Params(K_in=tem_K, P_in=P_val, nodes_in=num_samplingnodes, &
                               Epah_in=0.02d0, lambda_in=0.90d0, sigma_in=0.1d0, noise_var=sigma_noise, &
                               rank=ranks, log_dir=output_dir, tag_K=str_K, tag_P=str_P)
            
        file_limit = trim(output_dir) // '/Limit_K-'//trim(str_K)//'_P-'//trim(str_P)//'.txt'
        call Execute_Limit_Search_Workflow(ranks, noise_type, back_index, &
                                            7.5d0, 9.0d0, Limit, arl, file_limit, .false.)
        call mMSTD_ICPerformance(noise_type, back_index, Limit, falserate)
        file_ICp = trim(output_dir)  // '/ICperformance_K-'//trim(str_K)//'_P-'//trim(str_P)//'.txt'
        if (ranks == 0) then
            open(unit=70, file=file_ICp, status='unknown', position='append', action='write')
            write(70, '(F8.6, 1X)', advance='no') falserate
            close(70)
        end if 
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        file_res = trim(output_dir) // '/OC_K-'//trim(str_K)//'_P-'//trim(str_P)//'.txt'
        if (ranks == 0) call Write_OCARL_Header(file_res, "a b")
        
        do s = 1,  SIZE(shift_radius)
            cur_rad = shift_radius(s)
            do k = 1, SIZE(shift_size)
                cur_amp = shift_size(k)
                call mMSTD_OCARL(noise_type, back_index, anomaly_type, Limit=Limit, &
                                  shiftradius_zero=0.02d0, delta_radius=cur_rad, delta_time=cur_amp, &
                                  res_OC_ARL=res_OC, res_Exploration_ARL=res_Exp, &
                                  res_Exploitation_Num=res_Exploit)
                if (ranks == 0) call Write_OCARL_Row(file_res, (/ cur_amp, cur_rad /), res_OC, res_Exp, res_Exploit)
            end do
        end do
        call MPI_Barrier(MPI_COMM_WORLD, ierr) 

    
    call Clean_Performance_Evaluation()
    if (allocated(Nodesset)) deallocate(Nodesset)
    if (allocated(NeighborDis)) deallocate(NeighborDis)
    if (allocated(NeighborId)) deallocate(NeighborId)
    if (allocated(Kernel_NeighborDis)) deallocate(Kernel_NeighborDis)
    if (allocated(Kernel_NeighborId)) deallocate(Kernel_NeighborId)

    if (ranks == 0) print *, "=== Circle Test Completed ==="

end subroutine Run_mMSTD_ST_comparison

!==============================================================
!  (4) Ellipse  Test
!==============================================================
subroutine Run_mMSTD_Ellipse_Test(base_output_dir)
    implicit none
    character(len=*), intent(in) :: base_output_dir
    
    integer :: ranks, num_procs, k, s, i, ierr
    real(dp) :: Limit, arl, cur_amp, cur_area, cur_ellipse
    type(ARL_Stats_type) :: res_OC, res_Exp, res_Exploit
    
    integer :: noise_type, back_index, anomaly_type
    real(dp) :: tem_K, P_val
    character(len=1024) :: output_dir, file_limit, file_res
    character(len=32) :: str_K, str_P
    
    real(dp), parameter :: shift_value(2) = (/ 0.5d0, 2.0d0 /)
    real(dp), parameter :: shift_area(2)   = (/ 0.02d0, 0.03d0 /)

    real(dp), allocatable :: Nodesset(:,:), NeighborDis(:,:), Kernel_NeighborDis(:,:)
    integer, allocatable :: NeighborId(:,:), Kernel_NeighborId(:,:)

    call Init_Experiment_Common(ranks, num_procs)
    output_dir = trim(base_output_dir) // '/mMKSTD/ELLIPSE-2'
    call ensure_directory(output_dir)
    if (ranks == 0) call Print_Banner("Gauss + ELLIPSE", output_dir)

    call set_generalparams(num_samplingnodes_in=50)
    noise_type = NOISE_GAUSSIAN; back_index = BACK_NONE; anomaly_type = TYPE_ELLIPSE
    tem_K = 20.0d0; P_val = 0.90d0
    
    write(str_K, '(F4.1)') tem_K; str_K = trim(adjustl(str_K))
    write(str_P, '(F4.2)') P_val; str_P = trim(adjustl(str_P))

    call Calc_mMSTD_Params(K_in=tem_K, P_in=P_val, nodes_in=num_samplingnodes, &
                           Epah_in=0.02d0, lambda_in=0.90d0, sigma_in=0.1d0, noise_var=sigma_noise, &
                           rank=ranks, log_dir=output_dir, tag_K=str_K, tag_P=str_P)

    call Prepare_Simulation_Data(ranks, cache_path('Prepared_settings-1'), &
        Nodesset=Nodesset, NeighborDis=NeighborDis, NeighborId=NeighborId, &
        Kernel_NeighborDis=Kernel_NeighborDis, Kernel_NeighborId=Kernel_NeighborId)
    call Init_Performance_Evaluation(NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId, Nodesset)

    ! --- In-control charting limit ---
    Limit = 4.31769021472954d0
    ! --- OCARL Loop ---
    file_res = trim(output_dir) // '/OC_K-'//trim(str_K)//'_P_value-'//trim(str_P)//'.txt'
    if (ranks == 0) call Write_OCARL_Header(file_res, "Amp Area Ellipse")

    do k = 1, SIZE(shift_value)
        cur_amp = shift_value(k) 
        do s = 1, SIZE(shift_area)
            cur_area = shift_area(s)
            do i = 20, 1, -1
                cur_ellipse = 0.05d0 * dble(i)
                call mMSTD_OCARL(noise_type, back_index, anomaly_type, &
                                 Limit=Limit, shiftarea=cur_area, shiftvalue=cur_amp, shiftellipticity=cur_ellipse, &
                                 res_OC_ARL=res_OC, res_Exploration_ARL=res_Exp, res_Exploitation_Num=res_Exploit)
        
                if (ranks == 0) call Write_OCARL_Row(file_res, (/ cur_amp, cur_area, cur_ellipse /), &
                                                     res_OC, res_Exp, res_Exploit)
            end do
        end do
    end do
    
    call Clean_Performance_Evaluation()
    if (allocated(Nodesset)) deallocate(Nodesset)
    if (allocated(NeighborDis)) deallocate(NeighborDis)
    if (allocated(NeighborId)) deallocate(NeighborId)
    if (allocated(Kernel_NeighborDis)) deallocate(Kernel_NeighborDis)
    if (allocated(Kernel_NeighborId)) deallocate(Kernel_NeighborId)
    
    if (ranks == 0) print *, "===  Test Completed ==="

end subroutine Run_mMSTD_Ellipse_Test

!==============================================================
!  (4)Crescent Test
!==============================================================
subroutine Run_mMSTD_Crescent_Test(base_output_dir)
    implicit none
    character(len=*), intent(in) :: base_output_dir
    
    integer :: ranks, num_procs, k, s, i, ierr
    real(dp) :: Limit, arl, cur_amp, cur_area, cur_crescent
    type(ARL_Stats_type) :: res_OC, res_Exp, res_Exploit
    
    integer :: noise_type, back_index, anomaly_type
    real(dp) :: tem_K, P_val
    character(len=1024) :: output_dir, file_limit, file_res
    character(len=32) :: str_K, str_P
    
    real(dp), parameter :: shift_value(2) = (/ 0.5d0, 2.0d0 /)
    real(dp), parameter :: shift_area(2)   = (/ 0.02d0, 0.03d0 /)

    real(dp), allocatable :: Nodesset(:,:), NeighborDis(:,:), Kernel_NeighborDis(:,:)
    integer, allocatable :: NeighborId(:,:), Kernel_NeighborId(:,:)

    call Init_Experiment_Common(ranks, num_procs)
    output_dir = trim(base_output_dir) // '/mMKSTD/Crescent-2'
    call ensure_directory(output_dir)
    if (ranks == 0) call Print_Banner("Gauss + CRESCENT", output_dir)

    call set_generalparams(num_samplingnodes_in=50)
    noise_type = NOISE_GAUSSIAN; back_index = BACK_NONE; anomaly_type = TYPE_CRESCENT
    tem_K = 20.0d0; P_val = 0.90d0
    
    write(str_K, '(F4.1)') tem_K; str_K = trim(adjustl(str_K))
    write(str_P, '(F4.2)') P_val; str_P = trim(adjustl(str_P))

    call Calc_mMSTD_Params(K_in=tem_K, P_in=P_val, nodes_in=num_samplingnodes, &
                           Epah_in=0.02d0, lambda_in=0.90d0, sigma_in=0.1d0, noise_var=sigma_noise, &
                           rank=ranks, log_dir=output_dir, tag_K=str_K, tag_P=str_P)

    call Prepare_Simulation_Data(ranks, cache_path('Prepared_settings-1'), &
        Nodesset=Nodesset, NeighborDis=NeighborDis, NeighborId=NeighborId, &
        Kernel_NeighborDis=Kernel_NeighborDis, Kernel_NeighborId=Kernel_NeighborId)
    call Init_Performance_Evaluation(NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId, Nodesset)

    ! --- In-control charting limit ---
    Limit = 4.31769021472954d0
    ! --- OCARL Loop ---
    file_res = trim(output_dir) // '/OC_K-'//trim(str_K)//'_P_value-'//trim(str_P)//'.txt'
    if (ranks == 0) call Write_OCARL_Header(file_res, "Amp Area Crescent")

    do k = 1, SIZE(shift_value)
        cur_amp = shift_value(k) 
        do s = 1, SIZE(shift_area)
            cur_area = shift_area(s)
            do i = 20, 1, -1
                cur_crescent = 0.1d0 * dble(i)
                call mMSTD_OCARL(noise_type, back_index, anomaly_type, &
                                 Limit=Limit, shiftarea=cur_area, shiftvalue=cur_amp, shiftellipticity=cur_crescent, &
                                 res_OC_ARL=res_OC, res_Exploration_ARL=res_Exp, res_Exploitation_Num=res_Exploit)
        
                if (ranks == 0) call Write_OCARL_Row(file_res, (/ cur_amp, cur_area, cur_crescent /), &
                                                     res_OC, res_Exp, res_Exploit)
            end do
        end do
    end do
    
    call Clean_Performance_Evaluation()
    if (allocated(Nodesset)) deallocate(Nodesset)
    if (allocated(NeighborDis)) deallocate(NeighborDis)
    if (allocated(NeighborId)) deallocate(NeighborId)
    if (allocated(Kernel_NeighborDis)) deallocate(Kernel_NeighborDis)
    if (allocated(Kernel_NeighborId)) deallocate(Kernel_NeighborId)
    
    if (ranks == 0) print *, "===  Test Completed ==="

end subroutine Run_mMSTD_Crescent_Test

!==============================================================
!  (5) Circle Qnum Test 
!==============================================================
subroutine Run_mMSTD_Circle_QnumTest(base_output_dir)
    implicit none
    character(len=*), intent(in) :: base_output_dir

    integer :: ranks, num_procs, q, k, s, Qnum, ierr
    real(dp) :: Limit, arl, cur_amp, cur_rad
    type(ARL_Stats_type) :: res_OC, res_Exp, res_Exploit
    
    integer :: noise_type, back_index, anomaly_type
    real(dp) :: tem_K, P_val
    character(len=1024) :: output_dir, file_limit, file_res
    character(len=32) :: str_K, str_P, str_Qnum

    real(dp) :: limit_vec(20)
    real(dp), parameter :: shift_value(3) = (/ 0.5d0, 1.0d0, 2.0d0 /)
    real(dp), parameter :: shift_rad(3)   = (/ 0.02d0, 0.05d0, 0.1d0 /)

    real(dp), allocatable :: Nodesset(:,:), NeighborDis(:,:), Kernel_NeighborDis(:,:)
    integer, allocatable :: NeighborId(:,:), Kernel_NeighborId(:,:)

    ! Initialize limit vectors (Manual array constructor compatible with F90)
    limit_vec = (/ 3.60392264899071d0, 3.84183112728113d0, 3.99252399732375d0, &
        4.06959054959716d0, 4.14345469228244d0, 4.20470463590402d0, &
        4.24539689527241d0, 4.27517580996728d0, 4.30948586860874d0, &
        4.33351926007100d0, 4.35856704326198d0, 4.37585741345452d0, &
        4.39887867309339d0, 4.41564224404610d0, 4.43237024646909d0, &
        4.44407322624305d0, 4.46676576862997d0, 4.47045851890357d0, &
        4.48921058299815d0, 4.50362366699675d0 /)
    call Init_Experiment_Common(ranks, num_procs)
    output_dir = trim(base_output_dir) // '/mMKSTD/Qnum'
    call ensure_directory(output_dir)
    if (ranks == 0) call Print_Banner("Circle Qnum Test", output_dir)
    call set_generalparams(num_x_in=100, num_y_in=100, num_disnearspatial_in=100)
    call Prepare_Simulation_Data(ranks, cache_path('Prepared_settings-1'), &
        Nodesset=Nodesset, NeighborDis=NeighborDis, NeighborId=NeighborId, &
        Kernel_NeighborDis=Kernel_NeighborDis, Kernel_NeighborId=Kernel_NeighborId)
    call Init_Performance_Evaluation(NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId, Nodesset)

    ! Setup Base Params
    tem_K = 10.0d0; P_val = 0.90d0
    noise_type = NOISE_GAUSSIAN; back_index = BACK_NONE; anomaly_type = TYPE_CIRCLE
    
    write(str_K, '(F4.1)') tem_K; str_K = trim(adjustl(str_K))
    write(str_P, '(F4.2)') P_val; str_P = trim(adjustl(str_P))

    ! --- Outer Loop (Qnum) ---
    do q = 1, 20
        Qnum = 5 * q
        call set_generalparams(num_samplingnodes_in=Qnum)
        write(str_Qnum, '(I3)') Qnum; str_Qnum = trim(adjustl(str_Qnum))
        
        call Calc_mMSTD_Params(K_in=tem_K, P_in=P_val, nodes_in=num_samplingnodes, &
                               Epah_in=0.02d0, lambda_in=0.90d0, sigma_in=0.1d0, noise_var=sigma_noise, &
                               rank=ranks, log_dir=output_dir, tag_K=str_K, tag_P=str_P)
                               
        Limit = limit_vec(q)
        ! --- OCARL Loop ---
        file_res = trim(output_dir) // '/OC_K-'//trim(str_K)//'_P-'//trim(str_P)//'_Q-'//trim(str_Qnum)//'.txt'
        if (ranks == 0) call Write_OCARL_Header(file_res, "Amp Radius")

        do k = 1, SIZE(shift_value)
            cur_amp = shift_value(k) 
            do s = 1, SIZE(shift_rad)
                cur_rad = shift_rad(s)
                
                call mMSTD_OCARL(noise_type, back_index, anomaly_type, Limit=Limit, shiftradius=cur_rad, shiftvalue=cur_amp, &
                                 res_OC_ARL=res_OC, res_Exploration_ARL=res_Exp, res_Exploitation_Num=res_Exploit)
                
                if (ranks == 0) call Write_OCARL_Row(file_res, (/ cur_amp, cur_rad /), res_OC, res_Exp, res_Exploit)
            end do
        end do
        call MPI_Barrier(MPI_COMM_WORLD, ierr) 
    end do
    
    call Clean_Performance_Evaluation()
    if (allocated(Nodesset)) deallocate(Nodesset)
    if (allocated(NeighborDis)) deallocate(NeighborDis)
    if (allocated(NeighborId)) deallocate(NeighborId)
    if (allocated(Kernel_NeighborDis)) deallocate(Kernel_NeighborDis)
    if (allocated(Kernel_NeighborId)) deallocate(Kernel_NeighborId)
    
    if (ranks == 0) print *, "=== Circle Qnum Test Completed ==="

end subroutine Run_mMSTD_Circle_QnumTest


!==============================================================
!  PRIVATE HELPER SUBROUTINES
!==============================================================

! --- Helper 1: Initialize MPI Rank and Print ---
! --- Calculate mMKSTD parameters and write the optional configuration log ---
subroutine Calc_mMSTD_Params(K_in, P_in, nodes_in, Epah_in, lambda_in, sigma_in, noise_var, &
                             rank, log_dir, tag_K, tag_P)
    real(dp), intent(in) :: K_in, P_in, Epah_in, lambda_in, sigma_in,noise_var
    integer, intent(in) :: nodes_in, rank
    character(len=*), intent(in), optional :: log_dir, tag_K, tag_P
    
    ! Locals
    real(dp) :: C_val, delta_m, exp_num, rad_sq, rad, pho, alpha0, alpha1, alpha, u_t
    character(len=1024) :: file_conf
    integer :: stat

    ! 1. Math Logic
    C_val   = 1.0d0 - (1.0d0 - P_in) / 2.0d0
    delta_m = ANORIN(C_val)*noise_var
    exp_num = (P_in + (1.0_dp - P_in)*2.0_dp / K_in) * real(nodes_in, dp)
    rad_sq  = 2.0_dp / (theta_coverage * exp_num * pi_value * K_in)
    pho     = rad_sq / K_in
    
    alpha0  = (1.0_dp + sigma_in) * log((rad_sq / pho) + 1.0_dp) 
    u_t     = (1.0_dp - ((0.5)**2.0))**2.0
    alpha1  = (u_t + sigma_in)/u_t * LOG((rad_sq + pho)/(((Epah_in/2.0d0)**2.0) + pho))
    alpha   = MIN(alpha0, alpha1)/(delta_m**2.0_dp)

    if (rank == 0) print *, ">> Calc Params: K=", K_in, " P=", P_in, " Alpha=", alpha

    ! 2. Set Global Params
    call set_mMSTDparams(Epah_mM_in=Epah_in, alpha_mM_in=alpha, &
                         lambda_mM_in=lambda_in, pho_mM_in=pho, &
                         sigmaprio_mM_in=sigma_in)
                         
    ! 3. Write Config File (Optional)
    if (present(log_dir) .and. present(tag_K) .and. rank == 0) then
        file_conf = trim(log_dir) // '/Config_K-' // trim(tag_K) // '_P-' // trim(tag_P) // '.txt'
        open(unit=50, file=file_conf, status='replace', action='write', iostat=stat)
        if (stat == 0) then
            write(50, '(A)') "### AUTOMATED CONFIG LOG ###"
            write(50, '(A,F8.4)') " K: ", K_in
            write(50, '(A,F8.4)') " P: ", P_in
            write(50, '(A,F12.6)')" Pho:", pho
            write(50, '(A,F12.6)')" Alpha:", alpha
            call get_mMSTDparams(50)
            close(50)
        end if
    end if
end subroutine Calc_mMSTD_Params

! --- Helper 3: Encapsulate Limit Search Workflow ---
subroutine Execute_Limit_Search_Workflow(rank, n_type, b_idx, bound_L, bound_R, &
                                         limit_out, arl_out, filename, auto_range)
    integer, intent(in) :: rank, n_type, b_idx
    real(dp), intent(in) :: bound_L, bound_R
    real(dp), intent(out) :: limit_out, arl_out
    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: auto_range
    
    real(dp) :: L, R, std_dummy
    logical :: use_auto
    
    use_auto = .false.
    if (present(auto_range)) use_auto = auto_range
    
    if (use_auto) then
        if (rank == 0) print *, ">> Estimating Limit Range..."
        call mMSTD_EstimateLimitRange(n_type, b_idx, limit_out, L, R)
    else
        L = bound_L
        R = bound_R
    end if

    if (rank == 0) open(10, file=filename, status='replace')
    
    call mMSTD_limitSearch(n_type, b_idx, L, R, limit_out, arl_out, std_dummy, fid=10)
    
    if (rank == 0) then
        write(10, *) "Final Limit:", limit_out, " ARL:", arl_out
        close(10)
        print *, ">> Limit Found:", limit_out
    end if
end subroutine Execute_Limit_Search_Workflow

! --- Helper 4: Write OCARL Header (Unified) ---
! --- Helper 5: Write OCARL Data Row (Unified) ---
end module Experiment_Runner_mod
