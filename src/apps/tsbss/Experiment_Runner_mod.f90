!==============================================================
!  Module: Experiment_Runner_mod
!  Purpose: Define TSBSS experiment configurations and execution flow.
!==============================================================
module Experiment_Runner_mod
    use ExperimentSupport_mod, only: Init_Experiment_Common, Print_Banner, &
        Write_OCARL_Header, Write_OCARL_Row
    use mpi
    use GlobalSettings_mod
    use TSBSS_mod            
    use utils_mod
    use SamplesGeneration_mod
    use performance_evaluation 
    use Data_Loader 
    use RuntimeConfig_mod, only: cache_path, ensure_directory
    implicit none
    
    public

    ! Private helpers
    private :: Setup_TSBSS_Config, Execute_Limit_Search_Workflow_TSBSS, &
               Write_OCARL_Header, Write_OCARL_Row, Init_Experiment_Common, &
               Print_Banner

contains

!==============================================================
!  (1) B-Spline Test (TSBSS)
!==============================================================
subroutine Run_TSBSS_Bspline_Test(base_output_dir)
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
    type(Bbasis_sparse_type) :: B0_all, B1_all
    
    real(dp), parameter :: shift_size(2)    = (/8.0d0,10.0d0 /)

    ! --- 1. Init ---
    call Init_Experiment_Common(ranks, num_procs)
    output_dir = trim(base_output_dir) // '/TSBSS/BsplineNoBack-3'
    call ensure_directory(output_dir)
    if (ranks == 0) call Print_Banner("TSBSS: B-Spline Test", output_dir)

    ! --- 2. Setup Params ---
    noise_type = NOISE_GAUSSIAN; back_index =  BACK_NONE; anomaly_type = TYPE_BSPLINE
    
    call set_generalparams(num_samplingnodes_in=30)
    call Setup_TSBSS_Config(ranks, output_dir, "Bspline", &
                            0.9d0, 0.01d0, 5.0d0, 1.0d0, 3.0d0, &
                            noise_type, back_index, anomaly_type)

    ! --- 3. Load Data ---
    call Prepare_Simulation_Data(ranks, cache_path('Prepared_settings-1'), &
        B0_out=B0_all, B1_out=B1_all)
    call Init_Performance_Evaluation(B0_all=B0_all, B1_all=B1_all)

    ! --- 4. In-control charting limit ---
    Limit = 8.06767069183581d0
    ! --- 5. OCARL Experiment ---
    file_res = trim(output_dir) // '/OC_Bspline-add.txt'
    if (ranks == 0) call Write_OCARL_Header(file_res, "Amp")

    do k = 1, SIZE(shift_size)
        cur_amp = shift_size(k)  
        call TSBSS_OCARL(noise_type, back_index, anomaly_type, &
                         Limit=Limit, amplitude=cur_amp, &
                         res_OC_ARL=res_OC, res_Exploration_ARL=res_Exp, res_Exploitation_Num=res_Exploit)
        
        if (ranks == 0) call Write_OCARL_Row(file_res, (/ cur_amp /), res_OC, res_Exp, res_Exploit)
    end do

    ! --- 6. Cleanup ---
    call Clean_Performance_Evaluation()
    call free_Bbasis_sparse(B0_all); call free_Bbasis_sparse(B1_all)
    
    if (ranks == 0) print *, "=== TSBSS B-Spline Completed ==="

end subroutine Run_TSBSS_Bspline_Test


!==============================================================
!  PRIVATE HELPER SUBROUTINES
!==============================================================

! --- Reuse Common Helpers ---
! --- Write Headers/Rows ---
! --- TSBSS Specific Helper: Config Setup ---
subroutine Setup_TSBSS_Config(rank, log_dir, tag, lam0, w0, s0, se, sb, &
                              n_type, b_idx, a_type)
    integer, intent(in) :: rank, n_type, b_idx, a_type
    real(dp), intent(in) :: lam0, w0, s0, se, sb
    character(len=*), intent(in) :: log_dir, tag
    
    character(len=1024) :: fname
    integer :: stat
    
    call set_TSBSSparams(lambda0_TS_in=lam0, w0_TS_in=w0, &
                         sigma0_TS_in=s0, sigmae_TS_in=se, &
                         sigmab_TS_in=sb)
                         
    if (rank == 0) then
        fname = trim(log_dir) // '/Config_' // trim(tag) // '.txt'
        open(unit=50, file=fname, status='replace', action='write', iostat=stat)
        if (stat == 0) then
            write(50, '(A)') "### TSBSS CONFIG ###"
            write(50, '(A,F8.4)') " Lambda0:", lam0
            write(50, '(A,F8.4)') " W0:     ", w0
            write(50, '(A,F8.4)') " Sigma0: ", s0
            call get_TSBSSparams(50)
            write(50, *)
            write(50, '(A,I4)') " Noise: ", n_type
            write(50, '(A,I4)') " Back:  ", b_idx
            write(50, '(A,I4)') " Ano:   ", a_type
            close(50)
        end if
    end if
end subroutine Setup_TSBSS_Config

! --- TSBSS Specific Helper: Limit Search Workflow ---
subroutine Execute_Limit_Search_Workflow_TSBSS(rank, n_type, b_idx, bound_L, bound_R, &
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
        if (rank == 0) print *, ">> Estimating Limit Range (TSBSS)..."
        call TSBSS_EstimateLimitRange(n_type, b_idx, L, R)
    else
        L = bound_L
        R = bound_R
    end if

    if (rank == 0) open(10, file=filename, status='replace')
    
    call TSBSS_limitSearch(n_type, b_idx, L, R, limit_out, arl_out, std_dummy, fid=10)
    
    if (rank == 0) then
        write(10, *) "Final Limit:", limit_out, " ARL:", arl_out
        close(10)
        print *, ">> Limit Found:", limit_out
    end if
end subroutine Execute_Limit_Search_Workflow_TSBSS

end module Experiment_Runner_mod
