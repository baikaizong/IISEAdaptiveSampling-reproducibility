!==============================================================
!  Module: Experiment_Runner_mod
!  Purpose: Define TSSPR experiment configurations and execution flow.
!==============================================================
module Experiment_Runner_mod
    use ExperimentSupport_mod, only: Init_Experiment_Common, Print_Banner, &
        Write_OCARL_Header, Write_OCARL_Row
    use mpi
    use GlobalSettings_mod
    use TSSPR_mod            
    use utils_mod
    use SamplesGeneration_mod
    use performance_evaluation 
    use Data_Loader 
    use RuntimeConfig_mod, only: cache_path, ensure_directory
    implicit none
    
    public

    ! Private helpers
    private :: Setup_TSSPR_Config, Execute_Limit_Search_Workflow_TSSPR, &
               Write_OCARL_Header, Write_OCARL_Row, Init_Experiment_Common, &
               Print_Banner

contains

!==============================================================
!  (1) B-Spline Test
!==============================================================
subroutine Run_TSSPR_Bspline_Test(base_output_dir)
    implicit none
    character(len=*), intent(in) :: base_output_dir
    
    ! Locals
    integer :: ranks, num_procs, k, ierr
    real(dp) :: Limit, arl, cur_amp
    type(ARL_Stats_type) :: res_OC, res_Exp, res_Exploit
    
    ! Config
    integer :: noise_type, back_index, anomaly_type
    character(len=1024) :: output_dir, file_limit, file_res

    ! Spatial data
    real(dp), allocatable :: Nodesset(:,:)
    type(Bbasis_sparse_type) :: B0_all, B1_all

    ! --- 1. Init ---
    call Init_Experiment_Common(ranks, num_procs)
    output_dir = trim(base_output_dir) // '/TSSPR/Bspline'
    call ensure_directory(output_dir)
    if (ranks == 0) call Print_Banner("TSSPR: B-Spline Test", output_dir)

    ! --- 2. Setup Params ---
    noise_type = NOISE_GAUSSIAN; back_index = BACK_BSPLINE; anomaly_type = TYPE_BSPLINE
    
    call set_generalparams()
    call Setup_TSSPR_Config(ranks, output_dir, "Standard", &
                            1.0d0, 1, &
                            noise_type, back_index, anomaly_type)

    ! --- 3. Load Data ---
    call Prepare_Simulation_Data(ranks, cache_path('Prepared_settings'), &
         Nodesset=Nodesset, B0_out=B0_all, B1_out=B1_all)
    call Init_Performance_Evaluation(Nodesset=Nodesset, B0_all=B0_all, B1_all=B1_all)

    ! --- 4. Limit Search ---
    call Set_SampleParams(noise_type_in=noise_type, back_index_in=back_index)
    file_limit = trim(output_dir) // '/Limit_BsplineBack_Gauss.txt'
    
    call Execute_Limit_Search_Workflow_TSSPR(ranks, noise_type, back_index, &
                                             0.0d0, 0.0d0, Limit, arl, file_limit, .true.) 
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    ! --- 5. OCARL Experiment ---
    file_res = trim(output_dir) // '/OC_BsplineAno_VaryAmp.txt'
    if (ranks == 0) call Write_OCARL_Header(file_res, "Amp")

    do k = 1, 20
        cur_amp = 0.5d0 * dble(k) 
        
        call TSSPR_OCARL(noise_type, back_index, anomaly_type, &
                         Limit=Limit, amplitude=cur_amp, &
                         res_OC_ARL=res_OC, res_Exploration_ARL=res_Exp, res_Exploitation_Num=res_Exploit)
        
        if (ranks == 0) call Write_OCARL_Row(file_res, (/ cur_amp /), res_OC, res_Exp, res_Exploit)
    end do

    ! --- 6. Cleanup ---
    call Clean_Performance_Evaluation()
    call free_Bbasis_sparse(B0_all); call free_Bbasis_sparse(B1_all)
    if (allocated(Nodesset)) deallocate(Nodesset)
    
    if (ranks == 0) print *, "=== TSSPR B-Spline Completed ==="

end subroutine Run_TSSPR_Bspline_Test


!==============================================================
!  (2) Kernel CircleConv Test
!==============================================================
subroutine Run_TSSPR_Kernel_CircleConv_Test(base_output_dir)
    implicit none
    character(len=*), intent(in) :: base_output_dir

    integer :: ranks, num_procs, i, ierr
    real(dp) :: Limit, arl, cur_amp, cur_rad
    type(ARL_Stats_type) :: res_OC, res_Exp, res_Exploit
    
    integer :: noise_type, back_index, anomaly_type
    character(len=1024) :: output_dir, file_limit, file_res_rad, file_res_amp
    
    real(dp), allocatable :: Nodesset(:,:)
    type(Bbasis_sparse_type) :: B0_all, B1_all

    call Init_Experiment_Common(ranks, num_procs)
    output_dir = trim(base_output_dir) // '/TSSPR/Kernel_CircleConv'
    call ensure_directory(output_dir)
    if (ranks == 0) call Print_Banner("TSSPR: Kernel + CircleConv", output_dir)

    ! --- Setup ---
    noise_type = NOISE_GAUSSIAN; back_index = BACK_KERNEL_RANDOM; anomaly_type = TYPE_CIRCLE_CONV
    call set_generalparams()
    call Setup_TSSPR_Config(ranks, output_dir, "KernelTest", &
                            1.0d0, 1, &
                            noise_type, back_index, anomaly_type)

    ! --- Load Data ---
    call Prepare_Simulation_Data(ranks, cache_path('Prepared_settings'), &
         Nodesset=Nodesset, B0_out=B0_all, B1_out=B1_all)
    call Init_Performance_Evaluation(Nodesset=Nodesset, B0_all=B0_all, B1_all=B1_all)

    ! --- Limit Search ---
    call Set_SampleParams(noise_type_in=noise_type, back_index_in=back_index)
    file_limit = trim(output_dir) // '/Limit_KernelBack_Gauss.txt'
    
    call Execute_Limit_Search_Workflow_TSSPR(ranks, noise_type, back_index, &
                                             100.0d0, 120.0d0, Limit, arl, file_limit, .false.)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    ! --- Exp A: Vary Radius ---
    file_res_rad = trim(output_dir) // '/OC_VaryRadius_FixShift2.0.txt'
    if (ranks == 0) call Write_OCARL_Header(file_res_rad, "Amp Radius")

    do i = 1, 20
        cur_rad = 0.01d0 * dble(i)
        cur_amp = 2.0d0
        call TSSPR_OCARL(noise_type, back_index, anomaly_type, Limit=Limit, shiftradius=cur_rad, amplitude=cur_amp, &
                         res_OC_ARL=res_OC, res_Exploration_ARL=res_Exp, res_Exploitation_Num=res_Exploit)
        if (ranks == 0) call Write_OCARL_Row(file_res_rad, (/ cur_amp, cur_rad /), res_OC, res_Exp, res_Exploit)
    end do

    ! --- Exp B: Vary Amp ---
    file_res_amp = trim(output_dir) // '/OC_VaryShift_FixRadius0.05.txt'
    if (ranks == 0) call Write_OCARL_Header(file_res_amp, "Amp Radius")

    do i = 1, 20
        cur_rad = 0.05d0
        cur_amp = 0.2d0 * dble(i) 
        call TSSPR_OCARL(noise_type, back_index, anomaly_type, Limit=Limit, shiftradius=cur_rad, amplitude=cur_amp, &
                         res_OC_ARL=res_OC, res_Exploration_ARL=res_Exp, res_Exploitation_Num=res_Exploit)
        if (ranks == 0) call Write_OCARL_Row(file_res_amp, (/ cur_amp, cur_rad /), res_OC, res_Exp, res_Exploit)
    end do

    call Clean_Performance_Evaluation()
    call free_Bbasis_sparse(B0_all); call free_Bbasis_sparse(B1_all)
    if (allocated(Nodesset)) deallocate(Nodesset)
    
    if (ranks == 0) print *, "=== TSSPR Kernel Test Completed ==="

end subroutine Run_TSSPR_Kernel_CircleConv_Test


!==============================================================
!  PRIVATE HELPER SUBROUTINES
!==============================================================

! --- Reuse Common Helpers ---
! --- Write Headers/Rows ---
! --- TSSPR Specific Helper: Config Setup ---
subroutine Setup_TSSPR_Config(rank, log_dir, tag, miu, topr, &
                              n_type, b_idx, a_type)
    integer, intent(in) :: rank, topr, n_type, b_idx, a_type
    real(dp), intent(in) :: miu
    character(len=*), intent(in) :: log_dir, tag
    
    character(len=1024) :: fname
    integer :: stat
    
    call set_TSSPRparams(miu_min_TS_in=miu, Topr_TS_in=topr)
                         
    if (rank == 0) then
        fname = trim(log_dir) // '/Config_' // trim(tag) // '.txt'
        open(unit=50, file=fname, status='replace', action='write', iostat=stat)
        if (stat == 0) then
            write(50, '(A)') "### TSSPR CONFIG ###"
            write(50, '(A,F8.4)') " Miu Min: ", miu
            write(50, '(A,I4)')   " Top R:   ", topr
            call get_TSSPRparams(50)
            write(50, *)
            write(50, '(A,I4)') " Noise: ", n_type
            write(50, '(A,I4)') " Back:  ", b_idx
            write(50, '(A,I4)') " Ano:   ", a_type
            close(50)
        end if
    end if
end subroutine Setup_TSSPR_Config

! --- TSSPR Specific Helper: Limit Search Workflow ---
subroutine Execute_Limit_Search_Workflow_TSSPR(rank, n_type, b_idx, bound_L, bound_R, &
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
        if (rank == 0) print *, ">> Estimating Limit Range (TSSPR)..."
        call TSSPR_EstimateLimitRange(n_type, b_idx, L, R)
    else
        L = bound_L
        R = bound_R
    end if

    if (rank == 0) open(10, file=filename, status='replace')
    
    call TSSPR_limitSearch(n_type, b_idx, L, R, limit_out, arl_out, std_dummy, fid=10)
    
    if (rank == 0) then
        write(10, *) "Final Limit:", limit_out, " ARL:", arl_out
        close(10)
        print *, ">> Limit Found:", limit_out
    end if
end subroutine Execute_Limit_Search_Workflow_TSSPR

end module Experiment_Runner_mod
