!==============================================================
!  Module: Experiment_Runner_mod
!  Purpose: Define SASAM experiment configurations and execution flow.
!==============================================================
module Experiment_Runner_mod
    use ExperimentSupport_mod, only: Init_Experiment_Common, Print_Banner, &
        Write_OCARL_Header, Write_OCARL_Row
    use mpi
    use GlobalSettings_mod
    use SASAM_mod            
    use utils_mod
    use SamplesGeneration_mod
    use performance_evaluation 
    use Data_Loader 
    use RuntimeConfig_mod, only: cache_path, ensure_directory
    implicit none
    
    public

    ! Private helpers to reduce code duplication
    private :: Setup_SASAM_Config, Execute_Limit_Search_Workflow_SASAM, &
               Write_OCARL_Header, Write_OCARL_Row, Init_Experiment_Common, &
               Print_Banner

contains

!==============================================================
!  (1) B-Spline Test
!==============================================================
subroutine Run_SASAM_Bspline_Test(base_output_dir)
    implicit none
    character(len=*), intent(in) :: base_output_dir
    
    ! Locals
    integer :: ranks, num_procs, k, ierr
    real(dp) :: Limit, arl, cur_amp
    type(ARL_Stats_type) :: res_OC, res_Exp, res_Exploit
    
    ! Config
    integer :: noise_type, back_index, anomaly_type
    character(len=1024) :: output_dir, file_limit, file_res
    real(dp), parameter :: shift_size(7)    = (/3.0d0, 4.0d0, 5.0d0, 8.0d0, 10.0d0,  15.0d0,20.0d0 /)
    ! Spatial data
    real(dp), allocatable :: Nodesset(:,:), NeighborDis(:,:), Kernel_NeighborDis(:,:)
    integer, allocatable :: NeighborId(:,:), Kernel_NeighborId(:,:)
    type(Bbasis_sparse_type) :: B0_all, B1_all

    ! --- 1. Init ---
    call Init_Experiment_Common(ranks, num_procs)
    output_dir = trim(base_output_dir) // '/SASAM/Bspline'
    call ensure_directory(output_dir)
    if (ranks == 0) call Print_Banner("SASAM: B-Spline Test", output_dir)

    ! --- 2. Setup Params ---
    noise_type = NOISE_GAUSSIAN; back_index = BACK_BSPLINE; anomaly_type = TYPE_BSPLINE
    
    call set_generalparams(num_samplingnodes_in=30)
    call Setup_SASAM_Config(ranks, output_dir, "Standard", &
                            0.02d0, SQRT(sigma_noise**2.0+sigma_ground**2.0), 0.1d0, 0.7d0, 1, &
                            noise_type, back_index, anomaly_type)

    ! --- 3. Load Data ---
    call Prepare_Simulation_Data(ranks, cache_path('Prepared_settings-1'), &
            Nodesset=Nodesset, Kernel_NeighborDis=Kernel_NeighborDis, Kernel_NeighborId=Kernel_NeighborId, &
            B0_out=B0_all, B1_out=B1_all)
    call Prepare_Spatial_Neighbors(ranks, cache_path('Prepared_settings-1'), &
            NeighborDis=NeighborDis, NeighborId=NeighborId, Nodesset=Nodesset)

    call Init_Performance_Evaluation(NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId, Nodesset, B0_all, B1_all)

    ! --- 4. Limit Search ---
    call Set_SampleParams(noise_type_in=noise_type, back_index_in=back_index)
    file_limit = trim(output_dir) // '/Limit_BsplineBack_Gauss.txt'
    
    call Execute_Limit_Search_Workflow_SASAM(ranks, noise_type, back_index, &
                                            60.0d0, 75.0d0, Limit, arl, file_limit, .false.) 
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    ! --- 5. OCARL Experiment ---
    file_res = trim(output_dir) // '/OC_BsplineAno_VaryAmp.txt'
    if (ranks == 0) call Write_OCARL_Header(file_res, "Amp")

    do k = 1, SIZE(shift_size)
        cur_amp = shift_size(k) 
        
        call SASAM_OCARL(noise_type, back_index, anomaly_type, &
                         Limit=Limit, amplitude=cur_amp, &
                         res_OC_ARL=res_OC, res_Exploration_ARL=res_Exp, res_Exploitation_Num=res_Exploit)
        
        if (ranks == 0) call Write_OCARL_Row(file_res, (/ cur_amp /), res_OC, res_Exp, res_Exploit)
    end do

    ! --- 6. Cleanup ---
    call Clean_Performance_Evaluation()
    call free_Bbasis_sparse(B0_all); call free_Bbasis_sparse(B1_all)
    if (allocated(Nodesset)) deallocate(Nodesset)
    if (ranks == 0) print *, "=== SASAM B-Spline Completed ==="

end subroutine Run_SASAM_Bspline_Test


!==============================================================
!  (2) Kernel CircleConv Test
!==============================================================
subroutine Run_SASAM_Kernel_CircleConv_Test(base_output_dir)
    implicit none
    character(len=*), intent(in) :: base_output_dir

    integer :: ranks, num_procs, i, ierr, k, s, shift_pointnum
    real(dp) :: Limit, arl, cur_amp, miu_min
    type(ARL_Stats_type) :: res_OC, res_Exp, res_Exploit
    
    integer :: noise_type, back_index, anomaly_type
    character(len=1024) :: output_dir, file_limit, file_res_rad, file_res_amp
    
    real(dp), allocatable :: Nodesset(:,:), NeighborDis(:,:), Kernel_NeighborDis(:,:)
    integer, allocatable :: NeighborId(:,:), Kernel_NeighborId(:,:)
    
        
    real(dp), parameter :: shift_size(6)    =  (/ 5.0d0, 10.0d0, 20.0d0, 30.0d0, 50.0d0,100.0d0 /)
    integer, parameter :: shift_num(3) = (/ 1, 2, 3 /)

    call Init_Experiment_Common(ranks, num_procs)
    output_dir = trim(base_output_dir) // '/SASAM/Kernel_CircleConv'
    call ensure_directory(output_dir)
    if (ranks == 0) call Print_Banner("SASAM: Kernel + CircleConv", output_dir)

    ! --- Setup ---
    noise_type = NOISE_GAUSSIAN; back_index = BACK_KERNEL_RANDOM; anomaly_type = TYPE_RANDOM_POINTS_CONV
    call set_generalparams(num_x_in=20, num_y_in=20, num_disnearspatial_in=50, num_samplingnodes_in=10)
    call compute_conv_stddev(miu_min)
    call Setup_SASAM_Config(ranks, output_dir, "KernelTest", &
                            0.1d0, miu_min, 0.1d0, 0.7d0, 1, &
                            noise_type, back_index, anomaly_type)

    ! --- Load Data ---
    call Prepare_Simulation_Data(ranks, cache_path('Prepared_settings-2'), &
            Nodesset=Nodesset, Kernel_NeighborDis=Kernel_NeighborDis, Kernel_NeighborId=Kernel_NeighborId)
    call Prepare_Spatial_Neighbors(ranks, cache_path('Prepared_settings-2'), &
            NeighborDis=NeighborDis, NeighborId=NeighborId, Nodesset=Nodesset)
    call Init_Performance_Evaluation(NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId, Nodesset)

    ! --- Limit Search ---
    call Set_SampleParams(noise_type_in=noise_type, back_index_in=back_index)
    file_limit = trim(output_dir) // '/Limit_KernelBack_Gauss-4.txt'
    
    call Execute_Limit_Search_Workflow_SASAM(ranks, noise_type, back_index, &
                                             200.0d0, 300.0d0, Limit, arl, file_limit, .false.)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    ! --- Exp A: Vary Radius ---
    file_res_rad = trim(output_dir) // '/OC_arl-4.txt'
    if (ranks == 0) call Write_OCARL_Header(file_res_rad, "Amp num_point")

    do s =  SIZE(shift_num), 1, -1
        shift_pointnum = shift_num(s)
        do k = SIZE(shift_size),1,-1
            cur_amp = shift_size(k)
            call SASAM_OCARL(noise_type, back_index, anomaly_type, Limit=Limit, shiftvalue=cur_amp, shift_pointnum=shift_pointnum,&
                                res_OC_ARL=res_OC, res_Exploration_ARL=res_Exp, res_Exploitation_Num=res_Exploit)
            if (ranks == 0) call Write_OCARL_Row(file_res_rad, (/ cur_amp, 1.0d0*shift_pointnum /), res_OC, res_Exp, res_Exploit)
        end do
    end do

    call Clean_Performance_Evaluation()
    if (ranks == 0) print *, "=== SASAM Kernel Test Completed ==="

end subroutine Run_SASAM_Kernel_CircleConv_Test


!==============================================================
!  (3) Ellipse  Test
!==============================================================
subroutine Run_SASAM_Ellipse_Test(base_output_dir)
    implicit none
    character(len=*), intent(in) :: base_output_dir
    
    integer :: ranks, num_procs, k, s, i, ierr
    real(dp) :: Limit, arl, cur_amp, cur_area, cur_ellipse
    type(ARL_Stats_type) :: res_OC, res_Exp, res_Exploit
    
    integer :: noise_type, back_index, anomaly_type
    character(len=1024) :: output_dir, file_limit, file_res
    
    real(dp), parameter :: shift_value(2) = (/ 0.5d0, 2.0d0 /)
    real(dp), parameter :: shift_area(2)   = (/ 0.02d0, 0.03d0 /)

    real(dp), allocatable :: Nodesset(:,:), NeighborDis(:,:), Kernel_NeighborDis(:,:)
    integer, allocatable :: NeighborId(:,:), Kernel_NeighborId(:,:)

    call Init_Experiment_Common(ranks, num_procs)
    output_dir = trim(base_output_dir) // '/SASAM/ELLIPSE-1'
    call ensure_directory(output_dir)
    if (ranks == 0) call Print_Banner("SASAM: Gauss + Ellipse", output_dir)

    call set_generalparams(num_samplingnodes_in=50)
    noise_type = NOISE_GAUSSIAN; back_index = BACK_NONE; anomaly_type = TYPE_ELLIPSE

    call Setup_SASAM_Config(ranks, output_dir, "Ellipse", &
                            0.02d0, 1.0d0, 0.1d0, 0.7d0, 1, &
                            noise_type, back_index, anomaly_type)

    ! --- Load Data ---
    call Prepare_Simulation_Data(ranks, cache_path('Prepared_settings-1'), &
            Nodesset=Nodesset, Kernel_NeighborDis=Kernel_NeighborDis, Kernel_NeighborId=Kernel_NeighborId)
    call Prepare_Spatial_Neighbors(ranks, cache_path('Prepared_settings-1'), &
            NeighborDis=NeighborDis, NeighborId=NeighborId, Nodesset=Nodesset)
    call Init_Performance_Evaluation(NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId, Nodesset)

    ! --- Charting limit ---
    file_limit = trim(output_dir) // '/Limit_Gauss.txt'
    
    ! The bounded search determines the charting limit.
    call Execute_Limit_Search_Workflow_SASAM(ranks, noise_type, back_index, &
                                             5.0d0, 8.0d0, Limit, arl, file_limit, .false.) 
    call MPI_Barrier(MPI_COMM_WORLD, ierr)    
    ! --- OCARL ---
    file_res = trim(output_dir) // '/OC_ELLIPSE.txt'
    if (ranks == 0) call Write_OCARL_Header(file_res, "Amp Area Ellipse")

    do k = 1, SIZE(shift_value)
        cur_amp = shift_value(k) 
        do s = 1, SIZE(shift_area)
            cur_area = shift_area(s)
            do i = 20, 1, -1
                cur_ellipse = 0.05d0 * dble(i)
                
                call SASAM_OCARL(noise_type, back_index, anomaly_type, &
                                 Limit=Limit, shiftarea=cur_area, shiftvalue=cur_amp, shiftellipticity=cur_ellipse, &
                                 res_OC_ARL=res_OC, res_Exploration_ARL=res_Exp, res_Exploitation_Num=res_Exploit)
                
                if (ranks == 0) call Write_OCARL_Row(file_res, (/ cur_amp, cur_area, cur_ellipse /), &
                                                     res_OC, res_Exp, res_Exploit)
            end do
        end do
    end do

    call Clean_Performance_Evaluation()
    if (ranks == 0) print *, "=== SASAM Ellipse Test Completed ==="

end subroutine Run_SASAM_Ellipse_Test

!==============================================================
!  (3) Crescent Test
!==============================================================
subroutine Run_SASAM_CRESCENT_Test(base_output_dir)
    implicit none
    character(len=*), intent(in) :: base_output_dir
    
    integer :: ranks, num_procs, k, s, i, ierr
    real(dp) :: Limit, arl, cur_amp, cur_area, cur_CRESCENT
    type(ARL_Stats_type) :: res_OC, res_Exp, res_Exploit
    
    integer :: noise_type, back_index, anomaly_type
    character(len=1024) :: output_dir, file_limit, file_res
    
    real(dp), parameter :: shift_value(2) = (/ 0.5d0, 2.0d0 /)
    real(dp), parameter :: shift_area(2)   = (/ 0.02d0, 0.03d0 /)

    real(dp), allocatable :: Nodesset(:,:), NeighborDis(:,:), Kernel_NeighborDis(:,:)
    integer, allocatable :: NeighborId(:,:), Kernel_NeighborId(:,:)

    call Init_Experiment_Common(ranks, num_procs)
    output_dir = trim(base_output_dir) // '/SASAM/CRESCENT-1'
    call ensure_directory(output_dir)
    if (ranks == 0) call Print_Banner("SASAM: Gauss + CRESCENT", output_dir)

    call set_generalparams(num_samplingnodes_in=50)
    noise_type = NOISE_GAUSSIAN; back_index = BACK_NONE; anomaly_type = TYPE_CRESCENT

    call Setup_SASAM_Config(ranks, output_dir, "CRESCENT", &
                            0.02d0, 1.0d0, 0.1d0, 0.7d0, 1, &
                            noise_type, back_index, anomaly_type)

    ! --- Load Data ---
    call Prepare_Simulation_Data(ranks, cache_path('Prepared_settings-1'), &
            Nodesset=Nodesset, Kernel_NeighborDis=Kernel_NeighborDis, Kernel_NeighborId=Kernel_NeighborId)
    call Prepare_Spatial_Neighbors(ranks, cache_path('Prepared_settings-1'), &
            NeighborDis=NeighborDis, NeighborId=NeighborId, Nodesset=Nodesset)
    call Init_Performance_Evaluation(NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId, Nodesset)

    ! --- Charting limit ---
    Limit = 6.68750000000000d0
    ! --- OCARL ---
    file_res = trim(output_dir) // '/OC_CRESCENT.txt'
    if (ranks == 0) call Write_OCARL_Header(file_res, "Amp Area CRESCENT")

    do k = 1, SIZE(shift_value)
        cur_amp = shift_value(k) 
        do s = 1, SIZE(shift_area)
            cur_area = shift_area(s)
            do i = 20, 1, -1
                cur_CRESCENT = 0.1d0 * dble(i)
                
                call SASAM_OCARL(noise_type, back_index, anomaly_type, &
                                 Limit=Limit, shiftarea=cur_area, shiftvalue=cur_amp, shiftellipticity=cur_CRESCENT, &
                                 res_OC_ARL=res_OC, res_Exploration_ARL=res_Exp, res_Exploitation_Num=res_Exploit)
                
                if (ranks == 0) call Write_OCARL_Row(file_res, (/ cur_amp, cur_area, cur_CRESCENT /), &
                                                     res_OC, res_Exp, res_Exploit)
            end do
        end do
    end do

    call Clean_Performance_Evaluation()
    if (ranks == 0) print *, "===  Test Completed ==="

end subroutine Run_SASAM_CRESCENT_Test

!==============================================================
!  (4) Exploration Test (Note: This uses different logic)
!==============================================================
subroutine Run_SASAM_Exploration(base_output_dir)
    implicit none
    character(len=*), intent(in) :: base_output_dir
    
    integer :: ranks, num_procs, ierr
    integer :: t_start, t_end, t_rate
    real(dp) :: t_duration, Limit
    integer, parameter :: FID_DATA = 30, FID_CONF = 50
    character(len=1024) :: output_dir, file_dis, file_config
    
    ! Spatial Data
    real(dp), allocatable :: Nodesset(:,:), NeighborDis(:,:), Kernel_NeighborDis(:,:), NeighborDis_m0(:,:)
    integer, allocatable :: NeighborId(:,:), Kernel_NeighborId(:,:), NeighborId_m0(:,:)

    call Init_Experiment_Common(ranks, num_procs)
    output_dir = trim(base_output_dir) // '/SASAM/Exploration-I'
    call ensure_directory(output_dir)
    if (ranks == 0) call Print_Banner("SASAM: Exploration", output_dir)

    file_dis = trim(output_dir) // '/Exploration_Random.txt'
    file_config = trim(output_dir) // '/Parameter_Settings.txt'

    call system_clock(count_rate=t_rate)
    call set_generalparams(num_samplingnodes_in=10)
    
    ! Config
    call Setup_SASAM_Config(ranks, output_dir, "Exploration", &
                            0.02d0, 1.0d0, 0.1d0, 0.7d0, 1, &
                            NOISE_GAUSSIAN, BACK_NONE, TYPE_CIRCLE)

    ! Load Data
    call Prepare_Simulation_Data(ranks, cache_path('Prepared_settings-1'), &
            Nodesset=Nodesset, NeighborDis=NeighborDis, NeighborId=NeighborId, &
            Kernel_NeighborDis=Kernel_NeighborDis, Kernel_NeighborId=Kernel_NeighborId)
    call Prepare_Spatial_Neighbors(ranks, cache_path('Prepared_settings-1'), &
            NeighborDis_m0, NeighborId_m0, Nodesset)
    
    call Init_Performance_Evaluation(NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId, Nodesset)
    call Init_Exploration_Data(NeighborDis_m0, NeighborId_m0, Nodesset)

    ! Execution
    Limit = 5.11381900097257d0
    
    if (ranks == 0) then
        open(unit=FID_DATA, file=file_dis, status='replace')    
        call system_clock(count=t_start)
    end if
            
    call Exploration_Evaluation(Limit=Limit, data_fid=FID_DATA, config_fid=FID_CONF)
    
    if (ranks == 0) then
        call system_clock(count=t_end)
        t_duration = real(t_end - t_start, dp) / real(t_rate, dp)
        write(FID_CONF, '(A, F10.2, A)') " Wall Time: ", t_duration, " sec"
        close(FID_DATA); close(FID_CONF)
    end if
    
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    
    call Clean_Performance_Evaluation()
    call Clean_Exploration_Data()
end subroutine Run_SASAM_Exploration


!==============================================================
!  (5) Circle Test (Varying Theta)
!==============================================================
subroutine Run_SASAM_Circle_Test(base_output_dir)
    implicit none
    character(len=*), intent(in) :: base_output_dir
    
    integer :: ranks, num_procs, i, d, k, s, ierr
    real(dp) :: Limit, arl, cur_amp, cur_rad
    type(ARL_Stats_type) :: res_OC, res_Exp, res_Exploit
    
    real(dp) :: temtheta_1, temtheta_2
    real(dp), parameter :: theta_1(4) = (/ 0.05d0, 0.10d0, 0.15d0, 0.20d0 /)
    real(dp), parameter :: theta_2(4) = (/ 0.10d0, 0.20d0, 0.50d0, 0.70d0 /)
    
    character(len=1024) :: output_dir, file_limit, file_res
    character(len=32) :: str_t1, str_t2
    
    real(dp), allocatable :: Nodesset(:,:), NeighborDis(:,:), Kernel_NeighborDis(:,:)
    integer, allocatable :: NeighborId(:,:), Kernel_NeighborId(:,:)

    call Init_Experiment_Common(ranks, num_procs)
    output_dir = trim(base_output_dir) // '/SASAM/Circle-1'
    call ensure_directory(output_dir)
    if (ranks == 0) call Print_Banner("SASAM: Gauss + Circle (Theta)", output_dir)

    call Prepare_Simulation_Data(ranks, cache_path('Prepared_settings-1'), &
            Nodesset=Nodesset, Kernel_NeighborDis=Kernel_NeighborDis, Kernel_NeighborId=Kernel_NeighborId)
    call Prepare_Spatial_Neighbors(ranks, cache_path('Prepared_settings-1'), &
            NeighborDis=NeighborDis, NeighborId=NeighborId, Nodesset=Nodesset)
    call Init_Performance_Evaluation(NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId, Nodesset)

    ! --- Outer Loops ---
    do i = 1, size(theta_1)
        temtheta_1 = theta_1(i)
        write(str_t1, '(F4.1)') temtheta_1; str_t1 = trim(adjustl(str_t1))
        
        do d = 1, size(theta_2) 
            temtheta_2 = theta_2(d)
            write(str_t2, '(F4.2)') temtheta_2; str_t2 = trim(adjustl(str_t2))
            
            call set_generalparams(num_samplingnodes_in=10)
            call Setup_SASAM_Config(ranks, output_dir, "Theta", &
                                    0.02d0, 1.0d0, temtheta_1, temtheta_2, 1, &
                                    NOISE_GAUSSIAN, BACK_NONE, TYPE_CIRCLE, &
                                    tag_custom="Config_t1-"//trim(str_t1)//"_t2-"//trim(str_t2)//".txt")

            ! Limit Search
            file_limit = trim(output_dir) // '/Limit_theta1-'//trim(str_t1)//'_theta2-'//trim(str_t2)//'.txt'
            call Execute_Limit_Search_Workflow_SASAM(ranks, NOISE_GAUSSIAN, BACK_NONE, &
                                                     3.0d0, 8.0d0, Limit, arl, file_limit, .false.)
            call MPI_Barrier(MPI_COMM_WORLD, ierr)

            ! OCARL
            file_res = trim(output_dir) // '/OC_theta1-'//trim(str_t1)//'_theta2-'//trim(str_t2)//'.txt'
            if (ranks == 0) call Write_OCARL_Header(file_res, "Amp Radius")

            do k = 10, 1, -1
                cur_amp = 0.2d0 * dble(k) 
                do s = 10, 1, -1
                    cur_rad = 0.02d0 * dble(s)
                    call SASAM_OCARL(NOISE_GAUSSIAN, BACK_NONE, TYPE_CIRCLE, &
                                     Limit=Limit, shiftradius=cur_rad, shiftvalue=cur_amp, &
                                     res_OC_ARL=res_OC, res_Exploration_ARL=res_Exp, res_Exploitation_Num=res_Exploit)
                    if (ranks == 0) call Write_OCARL_Row(file_res, (/ cur_amp, cur_rad /), res_OC, res_Exp, res_Exploit)
                end do
            end do
            call MPI_Barrier(MPI_COMM_WORLD, ierr) 
        end do 
    end do
    
    call Clean_Performance_Evaluation()
    if (ranks == 0) print *, "=== Circle Theta Test Completed ==="

end subroutine Run_SASAM_Circle_Test


!==============================================================
!  (6) Circle Qnum Test 
!==============================================================
subroutine Run_SASAM_Circle_QnumTest(base_output_dir)
    implicit none
    character(len=*), intent(in) :: base_output_dir

    integer :: ranks, num_procs, q, k, s, Qnum, ierr
    real(dp) :: Limit, arl, cur_amp, cur_rad
    type(ARL_Stats_type) :: res_OC, res_Exp, res_Exploit
    
    character(len=1024) :: output_dir, file_limit, file_res
    character(len=32) :: str_t1, str_t2, str_Qnum
    
    real(dp) :: t1 = 0.20d0, t2 = 0.70d0
    real(dp), parameter :: shift_value(3) = (/ 0.5d0, 1.0d0, 2.0d0 /)
    real(dp), parameter :: shift_rad(3)   = (/ 0.02d0, 0.05d0, 0.10d0 /)
    
    real(dp), allocatable :: Nodesset(:,:), NeighborDis(:,:), Kernel_NeighborDis(:,:)
    integer, allocatable :: NeighborId(:,:), Kernel_NeighborId(:,:)

    call Init_Experiment_Common(ranks, num_procs)
    output_dir = trim(base_output_dir) // '/SASAM/Circle-Q'
    call ensure_directory(output_dir)
    if (ranks == 0) call Print_Banner("SASAM: Circle Qnum Test", output_dir)

    call Prepare_Simulation_Data(ranks, cache_path('Prepared_settings-1'), &
            Nodesset=Nodesset, Kernel_NeighborDis=Kernel_NeighborDis, Kernel_NeighborId=Kernel_NeighborId)
    call Prepare_Spatial_Neighbors(ranks, cache_path('Prepared_settings-1'), &
            NeighborDis=NeighborDis, NeighborId=NeighborId, Nodesset=Nodesset)
    call Init_Performance_Evaluation(NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId, Nodesset)

    write(str_t1, '(F4.1)') t1; str_t1 = trim(adjustl(str_t1))
    write(str_t2, '(F4.2)') t2; str_t2 = trim(adjustl(str_t2))

    do q = 1, 20
        Qnum = 5 * q
        call set_generalparams(num_samplingnodes_in=Qnum)
        write(str_Qnum, '(I3)') Qnum; str_Qnum = trim(adjustl(str_Qnum))
        
        call Setup_SASAM_Config(ranks, output_dir, "QTest", &
                                0.02d0, 1.0d0, t1, t2, 1, &
                                NOISE_GAUSSIAN, BACK_NONE, TYPE_CIRCLE, &
                                tag_custom="Config_t1-"//trim(str_t1)//"_t2-"//trim(str_t2)//"_Q-"//trim(str_Qnum)//".txt")

        file_limit = trim(output_dir) // '/Limit_t1-'//trim(str_t1)//'_t2-'//trim(str_t2)//'_Q-'//trim(str_Qnum)//'.txt'
        call Execute_Limit_Search_Workflow_SASAM(ranks, NOISE_GAUSSIAN, BACK_NONE, &
                                                 3.0d0, 8.0d0, Limit, arl, file_limit, .false.)
        call MPI_Barrier(MPI_COMM_WORLD, ierr)

        file_res = trim(output_dir) // '/OC_t1-'//trim(str_t1)//'_t2-'//trim(str_t2)//'_Q-'//trim(str_Qnum)//'.txt'
        if (ranks == 0) call Write_OCARL_Header(file_res, "Amp Radius")

        do k = 1, SIZE(shift_value)
            cur_amp = shift_value(k) 
            do s = 1, SIZE(shift_rad)
                cur_rad = shift_rad(s)
                call SASAM_OCARL(NOISE_GAUSSIAN, BACK_NONE, TYPE_CIRCLE, &
                                 Limit=Limit, shiftradius=cur_rad, shiftvalue=cur_amp, &
                                 res_OC_ARL=res_OC, res_Exploration_ARL=res_Exp, res_Exploitation_Num=res_Exploit)
                if (ranks == 0) call Write_OCARL_Row(file_res, (/ cur_amp, cur_rad /), res_OC, res_Exp, res_Exploit)
            end do
        end do
        call MPI_Barrier(MPI_COMM_WORLD, ierr) 
    end do
    
    call Clean_Performance_Evaluation()
    if (ranks == 0) print *, "=== Circle Qnum Test Completed ==="

end subroutine Run_SASAM_Circle_QnumTest


!==============================================================
!  SASAM-specific helper subroutines
!==============================================================

! --- SASAM configuration setup ---
subroutine Setup_SASAM_Config(rank, log_dir, tag, epah, miu, t1, t2, topr, &
                              n_type, b_idx, a_type, tag_custom)
    integer, intent(in) :: rank, topr, n_type, b_idx, a_type
    real(dp), intent(in) :: epah, miu, t1, t2
    character(len=*), intent(in) :: log_dir, tag
    character(len=*), intent(in), optional :: tag_custom
    
    character(len=1024) :: fname
    integer :: stat
    
    call set_SASAMparams(Epah_SA_in=epah, miu_min_SA_in=miu, &
                         theta_1_SA_in=t1, theta_2_SA_in=t2, Topr_SA_in=topr)
                         
    if (rank == 0) then
        if (present(tag_custom)) then
            fname = trim(log_dir) // '/' // trim(tag_custom)
        else
            fname = trim(log_dir) // '/Config_' // trim(tag) // '.txt'
        end if
        
        open(unit=50, file=fname, status='replace', action='write', iostat=stat)
        if (stat == 0) then
            write(50, '(A)') "### SASAM CONFIG ###"
            write(50, '(A,F8.4)') " Epah: ", epah
            write(50, '(A,F8.4)') " Miu:  ", miu
            write(50, '(A,F8.4)') " T1:   ", t1
            write(50, '(A,F8.4)') " T2:   ", t2
            call get_SASAMparams(50)
            write(50, *)
            write(50, '(A,I4)') " Noise: ", n_type
            write(50, '(A,I4)') " Back:  ", b_idx
            write(50, '(A,I4)') " Ano:   ", a_type
            close(50)
        end if
    end if
end subroutine Setup_SASAM_Config

! --- SASAM Specific Helper: Limit Search Workflow ---
subroutine Execute_Limit_Search_Workflow_SASAM(rank, n_type, b_idx, bound_L, bound_R, &
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
        if (rank == 0) print *, ">> Estimating Limit Range (SASAM)..."
        limit_out = 50.0d0 ! Guess
        call SASAM_EstimateLimitRange(n_type, b_idx, limit_out, L, R)
    else
        L = bound_L
        R = bound_R
    end if

    if (rank == 0) open(10, file=filename, status='replace')
    
    call SASAM_limitSearch(n_type, b_idx, L, R, limit_out, arl_out, std_dummy, fid=10)
    
    if (rank == 0) then
        write(10, *) "Final Limit:", limit_out, " ARL:", arl_out
        close(10)
        print *, ">> Limit Found:", limit_out
    end if
end subroutine Execute_Limit_Search_Workflow_SASAM

end module Experiment_Runner_mod
