!==============================================================
!  Module: Experiment_Runner_mod
!  Purpose: Define CDS experiment configurations and execution flow.
!==============================================================
module Experiment_Runner_mod
    use ExperimentSupport_mod, only: Init_Experiment_Common, Print_Banner, &
        Write_OCARL_Header, Write_OCARL_Row
    use mpi
    use GlobalSettings_mod
    use CDS_mod              
    use utils_mod
    use SamplesGeneration_mod
    use performance_evaluation 
    use Data_Loader
    use RuntimeConfig_mod, only: cache_path, ensure_directory
    USE SVRGP_INT
    implicit none
    
    public

    ! Private helpers
    private :: Setup_CDS_Config, Execute_Limit_Search_Workflow_CDS, &
               Write_OCARL_Header, Write_OCARL_Row, Init_Experiment_Common, &
               Print_Banner

contains

!==============================================================
!  (1) Kernel CircleConv Test (CDS)
!==============================================================
subroutine Run_CDS_Kernel_CircleConv_Test(base_output_dir)
    implicit none
    character(len=*), intent(in) :: base_output_dir

    ! Locals
    integer :: ranks, num_procs, i, ierr, k, s
    real(dp) :: Limit, arl, cur_amp
    real(dp) :: miu_min
    type(ARL_Stats_type) :: res_OC, res_Exp, res_Exploit
    
    ! Config
    integer :: noise_type, back_index, anomaly_type, shift_pointnum
    character(len=1024) :: output_dir, file_limit, file_res_rad, file_res_amp,file_ICp

    ! Spatial data
    real(dp), allocatable :: Nodesset(:,:)
    
    
    
    ! CDS Data
    real(dp), allocatable :: Cov_diag(:)
    real(dp), allocatable :: Cov_Neighborval(:,:)
    real(dp), allocatable :: Cov_dense_global(:,:)
    integer, allocatable :: Cov_NeighborId(:,:)
    real(dp):: falserate

    real(dp), parameter :: shift_size(6)    = (/5.0d0, 10.0d0, 20.0d0,30.0d0, 50.0d0, 100.0d0 /)
    integer, parameter :: shift_num(3) = (/ 1, 2, 3 /)
    
    ! --- 1. Init ---
    call Init_Experiment_Common(ranks, num_procs)
    output_dir = trim(base_output_dir) // '/CDS/Kernel_CircleConv-6'
    call ensure_directory(output_dir)
    if (ranks == 0) call Print_Banner("CDS: Kernel + CircleConv", output_dir)

    ! --- 2. Setup Params ---
    noise_type = NOISE_GAUSSIAN; back_index = BACK_KERNEL_RANDOM; anomaly_type = TYPE_RANDOM_POINTS_CONV
    
    call set_generalparams(num_x_in=20, num_y_in=20, num_samplingnodes_in=10)
    call compute_conv_stddev(miu_min)
    call Setup_CDS_Config(ranks, output_dir, "KernelTest", &
                          miu_min, 1, 0.99d0, &
                          noise_type, back_index, anomaly_type)

    ! --- 3. Load Data ---
    ! CDS needs covariance data, so we use the extended Prepare_Simulation_Data call
    call Prepare_Simulation_Data(ranks, cache_path('Prepared_settings-3'), &
        Nodesset, &
        Kernel_Cov_NeighborId=Cov_NeighborId, Kernel_Cov_NeighborVal=Cov_Neighborval, &
        Kernel_Cov_diag=Cov_diag, Cov_dense_global=Cov_dense_global)

    call Init_Performance_Evaluation(Nodesset, &
        Cov_diag=Cov_diag, &
        Cov_Neighborval=Cov_Neighborval, &
        Cov_NeighborId=Cov_NeighborId, Cov_dense_global=Cov_dense_global)
    ! --- 4. Limit Search ---
    call Set_SampleParams(noise_type_in=noise_type, back_index_in=back_index)
    file_limit = trim(output_dir) // '/Limit_KernelBack_Gauss-3.txt'
    

    call Execute_Limit_Search_Workflow_CDS(ranks, noise_type, back_index, &
                                         118.0d0, 130.0d0, Limit, arl, file_limit )

    file_res_rad = trim(output_dir) // '/OC_Arl-3.txt'
    if (ranks == 0) call Write_OCARL_Header(file_res_rad, "Amp num_point")

    do s =  SIZE(shift_num), 1, -1
        shift_pointnum = shift_num(s)
        do k = SIZE(shift_size),1,-1
            cur_amp = shift_size(k)
            call CDS_OCARL(noise_type, back_index, anomaly_type, Limit=Limit, shiftvalue=cur_amp, shift_pointnum=shift_pointnum,&
                                res_OC_ARL=res_OC, res_Exploration_ARL=res_Exp, res_Exploitation_Num=res_Exploit)
            if (ranks == 0) call Write_OCARL_Row(file_res_rad, (/ cur_amp, 1.0d0*shift_pointnum /), res_OC, res_Exp, res_Exploit)
        end do
    end do

    ! --- 7. Cleanup ---
    call Clean_Performance_Evaluation()
    if (allocated(Nodesset)) deallocate(Nodesset)
    if (allocated(Cov_diag))        deallocate(Cov_diag)
    if (allocated(Cov_Neighborval)) deallocate(Cov_Neighborval)
    if (allocated(Cov_NeighborId))   deallocate(Cov_NeighborId)
    
    if (ranks == 0) print *, "=== CDS Kernel Test Completed ==="

end subroutine Run_CDS_Kernel_CircleConv_Test


!==============================================================
!  PRIVATE HELPER SUBROUTINES
!==============================================================

! --- Reuse Common Helpers ---
! --- Write Headers/Rows ---
! --- CDS Specific Helper: Config Setup ---
subroutine Setup_CDS_Config(rank, log_dir, tag, miu, topr, alpha, &
                            n_type, b_idx, a_type)
    integer, intent(in) :: rank, topr, n_type, b_idx, a_type
    real(dp), intent(in) :: miu, alpha
    character(len=*), intent(in) :: log_dir, tag
    
    character(len=1024) :: fname
    integer :: stat
    
    call set_CDSparams(miu_min_CDS_in=miu, Topr_CDS_in=topr, alpha_CDS_in=alpha)
                         
    if (rank == 0) then
        fname = trim(log_dir) // '/Config_' // trim(tag) // '.txt'
        open(unit=50, file=fname, status='replace', action='write', iostat=stat)
        if (stat == 0) then
            write(50, '(A)') "### CDS CONFIG ###"
            write(50, '(A,F8.4)') " Miu Min: ", miu
            write(50, '(A,F8.4)') " Alpha:   ", alpha
            write(50, '(A,I4)')   " Top R:   ", topr
            call get_CDSparams(50)
            write(50, *)
            write(50, '(A,I4)') " Noise: ", n_type
            write(50, '(A,I4)') " Back:  ", b_idx
            write(50, '(A,I4)') " Ano:   ", a_type
            close(50)
        end if
    end if
end subroutine Setup_CDS_Config

! --- CDS Specific Helper: Limit Search Workflow ---
subroutine Execute_Limit_Search_Workflow_CDS(rank, n_type, b_idx, bound_L, bound_R, &
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
        if (rank == 0) print *, ">> Estimating Limit Range (CDS)..."
        call cds_estimatelimitrange(n_type, b_idx, L, R)
    else
        L = bound_L
        R = bound_R
    end if

    if (rank == 0) open(10, file=filename, status='replace')
    
    call CDS_limitSearch(n_type, b_idx, L, R, limit_out, arl_out, std_dummy, fid=10)
    
    if (rank == 0) then
        write(10, *) "Final Limit:", limit_out, " ARL:", arl_out
        close(10)
        print *, ">> Limit Found:", limit_out
    end if
end subroutine Execute_Limit_Search_Workflow_CDS

end module Experiment_Runner_mod
