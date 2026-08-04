!==============================================================
!  Module: performance_evaluation
!  Purpose: Evaluate CDS detection limits and run-length performance.
!==============================================================
module performance_evaluation
    use GlobalSettings_mod
    use SparseMatrix_mod
    use SamplesGeneration_mod
    use CDS_mod              
    use utils_mod
    use mpi
    use RNSET_INT
    use RuntimeConfig_mod, only: get_clock_seed
    use RNUND_INT
    use SVRGP_INT
    implicit none
    public 
    
    logical, save :: perf_initialized = .false.

    contains
    
!==============================================================
!  Subroutine: Init_Performance_Evaluation
!  Inputs include the covariance structures required by CDS.
!==============================================================
subroutine Init_Performance_Evaluation( Nodesset,CovMat, &
    Cov_diag, Cov_Neighborval, Cov_NeighborId, Cov_dense_global) 
    
    implicit none

    
    ! Sample Gen Inputs
    real(dp), intent(in), optional :: Nodesset(:,:)
    type(Bbasis_sparse_type), intent(in), optional ::  CovMat

    ! CDS Specific Inputs (Optional if running other methods, but required for CDS)
    real(dp), intent(in), optional :: Cov_diag(:)
    real(dp), intent(in), optional :: Cov_Neighborval(:,:)
    real(dp), intent(in), optional :: Cov_dense_global(:,:)
    integer, intent(in), optional  :: Cov_NeighborId(:,:)

    ! 1. Init Sample Generation
    call Init_SamplesGeneration(Nodesset,  CovMat)

    ! 2. Init CDS Data (Check presence)
    if (present(Cov_diag) .and. present(Cov_Neighborval) .and. present(Cov_NeighborId) .and. present(Cov_dense_global)) then
        call Init_CDS_Data(Cov_diag_in=Cov_diag, &
                           Cov_Neighborval_in=Cov_Neighborval, &
                           Cov_NeighborId_in=Cov_NeighborId, &
                           Cov_dense_global_in=Cov_dense_global)
    end if

    
    perf_initialized = .true.

end subroutine Init_Performance_Evaluation

!==============================================================
!  Subroutine: Clean_Performance_Evaluation
!==============================================================
subroutine Clean_Performance_Evaluation()
    implicit none
    call Clean_CDS()
    call Clean_SamplesGeneration()
    perf_initialized = .false.
end subroutine Clean_Performance_Evaluation

!==============================================================
!  Subroutine: CDS_EstimateLimitRange
!  Purpose: Fast estimation of limit range (Single Pass)
!==============================================================
subroutine CDS_EstimateLimitRange(noise_type, back_index, &
                                    LimitL, LimitR)
    implicit none

    ! Inputs
    integer, intent(in) :: noise_type, back_index

    ! Outputs
    real(dp), intent(out) :: LimitL, LimitR

    ! MPI & Locals
    integer :: ierr, num_procs, my_rank, status(MPI_STATUS_SIZE)
    integer :: task, next_task, active_workers, worker_rank
    integer :: recv_count
    logical :: is_master
    integer, parameter :: WORK_TAG = 1, STOP_TAG = 2, RESULT_TAG = 3

    integer :: i, j, seed, time_seed
    real(dp) :: local_max_stat, sim_max_stat
    real(dp), allocatable :: all_max_stats(:)
    
    integer :: sampling_index(num_samplingnodes)
    real(dp) :: OnlineSample(num_samplingnodes)
    real(dp) :: charting_statistic
    
    integer :: sampling_mask(num_allnodes), nodes_count, one_vec(1)

    integer :: curr_percent, prev_percent
    integer, parameter :: bar_width = 40
    integer :: filled_len
    character(len=bar_width) :: bar_str

    if (.not. perf_initialized) stop "Error: Perf Eval not initialized."

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
    is_master = (my_rank == 0)

    if (is_master) then
        allocate(all_max_stats(IC_TestRuns))
        print *, "==== Auto-Estimating Limit Range (CDS) ===="
        prev_percent = -1 
    end if

    ! --- Parallel Execution ---
    if (is_master) then
        next_task = 1
        recv_count = 0
        active_workers = num_procs - 1
        
        do worker_rank = 1, num_procs - 1
            if (next_task <= IC_TestRuns) then
                call MPI_SEND(next_task, 1, MPI_INTEGER, worker_rank, WORK_TAG, MPI_COMM_WORLD, ierr)
                next_task = next_task + 1
            else
                call MPI_SEND(0, 1, MPI_INTEGER, worker_rank, STOP_TAG, MPI_COMM_WORLD, ierr)
                active_workers = active_workers - 1
            end if
        end do

        do while (active_workers > 0)
            call MPI_RECV(local_max_stat, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, RESULT_TAG, &
                        MPI_COMM_WORLD, status, ierr)
            worker_rank = status(MPI_SOURCE)
            
            recv_count = recv_count + 1
            if (recv_count <= IC_TestRuns) then
                all_max_stats(recv_count) = local_max_stat
                
                curr_percent = int(100.0_dp * dble(recv_count) / dble(IC_TestRuns))
                
    
                if (curr_percent > prev_percent) then
                    filled_len = int((real(curr_percent) / 100.0) * bar_width)
                    
                    bar_str = ""
                    if (filled_len > 0) bar_str(1:filled_len) = repeat('=', filled_len)
                    if (filled_len < bar_width) bar_str(filled_len+1:bar_width) = repeat('.', bar_width-filled_len)
                    

                    write(*, '(A1, A, A, A, I3, A)', advance='no') &
                        char(13), " Progress: [", bar_str, "] ", curr_percent, "%"
                        
                    prev_percent = curr_percent
                end if
            end if
            
            if (next_task <= IC_TestRuns) then
                call MPI_SEND(next_task, 1, MPI_INTEGER, worker_rank, WORK_TAG, MPI_COMM_WORLD, ierr)
                next_task = next_task + 1
            else
                call MPI_SEND(0, 1, MPI_INTEGER, worker_rank, STOP_TAG, MPI_COMM_WORLD, ierr)
                active_workers = active_workers - 1
            end if
        end do

        print *, "" 
        
    else
        ! Worker Logic
        do
            call MPI_RECV(task, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
            if (status(MPI_TAG) == STOP_TAG) exit
            
            call get_clock_seed(time_seed)
            seed = mod(time_seed + task*313 + my_rank*997, 2147483647) + 1
            call RNSET(seed)

            ! 1. Random Sampling Layout
            sampling_mask = 0; nodes_count = 0
            do while (nodes_count < num_samplingnodes)
                call RNUND(num_allnodes, one_vec)
                if (sampling_mask(one_vec(1)) == 0) then
                    nodes_count = nodes_count + 1
                    sampling_index(nodes_count) = one_vec(1)
                    sampling_mask(one_vec(1)) = 1
                end if
            end do

            ! 2. Reset CDS & Set Generation Params
            call Reset_CDS_State()
            call Set_SampleParams(ICorOC_in=IC, noise_type_in=noise_type, &
                                  back_index_in=back_index)
            do j = 1, Time_window
                call GenerateOnlineSample(sampling_index, OnlineSample)
                call CDS(OnlineSample, sampling_index, charting_statistic)
            end do            
            sim_max_stat = -1.0d20
            ! 3. Run Window
            do j = 1, IcARL
                call GenerateOnlineSample(sampling_index, OnlineSample)
                call CDS(OnlineSample, sampling_index, charting_statistic)
                if (charting_statistic > sim_max_stat) sim_max_stat = charting_statistic
            end do
            print*, " Worker ", my_rank, " Run ", task, " Max Stat: ", sim_max_stat
            call MPI_SEND(sim_max_stat, 1, MPI_DOUBLE_PRECISION, 0, RESULT_TAG, MPI_COMM_WORLD, ierr)
        end do
    end if
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    ! Master Calculation
    if (is_master) then
        ! Sort to find quantiles
        call partial_quickselect(values=all_max_stats, num_all=IC_TestRuns, &
                                 top_k=up_gap, &
                                 order_index=1, values_sub=all_max_stats)
        LimitL = all_max_stats(down_gap)
        LimitR = all_max_stats(up_gap) 

        print '(A,F10.4,A,F10.4,A)', " CDS Estimated Range: [", LimitL, ",", LimitR, "]"
        deallocate(all_max_stats)
    end if
    
    call MPI_BCAST(LimitL, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(LimitR, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

end subroutine CDS_EstimateLimitRange

!==============================================================
!  Subroutine: CDS_limitSearch
!  Cumulative simulation CPU time is recorded.
!==============================================================
subroutine CDS_limitSearch(noise_type, back_index, &
                           BLimitL, BLimitR, Limit, arl, std, fid)
    implicit none

    !---------------- Inputs ----------------
    integer, intent(in) :: noise_type, back_index
    integer, intent(in), optional :: fid
    real(dp), intent(in) :: BLimitL, BLimitR
    
    !---------------- Outputs ----------------
    real(dp), intent(out) :: Limit, arl, std

    !---------------- MPI Variables -----------------
    integer :: ierr, num_procs, my_rank, status(MPI_STATUS_SIZE)
    integer :: task, next_task, active_workers, worker_rank
    integer :: recv_count
    logical :: is_master, is_converged
    integer, parameter :: WORK_TAG = 1, STOP_TAG = 2, RESULT_TAG = 3

    !---------------- Locals -----------------
    integer :: seed, time_seed, i, j, print_step
    real(dp) :: running_sum, local_result, sim_result, global_arl
    real(dp) :: LimitL, LimitR, charting_statistic
    real(dp), allocatable :: all_results(:)
    
    integer :: sampling_index(num_samplingnodes)
    real(dp) :: OnlineSample(num_samplingnodes)
    integer :: sampling_mask(num_allnodes), nodes_count, one_vec(1)

    !---------------- Timing variables ----------------
    real(dp) :: t_sim_start, t_sim_end, local_sim_time
    real(dp) :: total_accumulated_time ! Sum of all simulation times
    
    ! Packet for sending [ARL_Result, Time_Taken]
    real(dp) :: send_packet(2), recv_packet(2)

    if (.not. perf_initialized) stop "Error: Perf Eval not initialized."

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
    is_master = (my_rank == 0)

    if (is_master) then
        LimitL = BLimitL
        LimitR = BLimitR
        global_arl = 0.0d0
        allocate(all_results(simu))
        print *, "==== Begin CDS_limitSearch (Cumulative Timing) ===="
        if (present(fid)) then
             write(fid, '(A)') "Iter_Limit,  Result_ARL,  Total_Calc_Time(s)"
        end if
    end if

    ! Binary Search Loop
    do         
      if (is_master) then        
            if (abs(global_arl - IcArl) < Icstd) then
                is_converged = .true.
            else
                is_converged = .false.
                Limit = 0.5d0 * (LimitL + LimitR)
                if (.not. present(fid)) then
                    print '(A,F12.6,A,F12.6,A,F12.6)', " Search: [", LimitL, ",", LimitR, "] Try:", Limit
                end if
            end if
        end if
        
        call MPI_BCAST(is_converged, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        if (is_converged) exit
        call MPI_BCAST(Limit, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        if (is_master) then
            next_task = 1
            recv_count = 0
            running_sum = 0.0d0   
            total_accumulated_time = 0.0d0 
            active_workers = num_procs - 1
            print_step = max(10, int(simu * 0.1)) 
            
            do worker_rank = 1, num_procs - 1
                if (next_task <= simu) then
                    call MPI_SEND(next_task, 1, MPI_INTEGER, worker_rank, WORK_TAG, MPI_COMM_WORLD, ierr)
                    next_task = next_task + 1
                else
                    call MPI_SEND(0, 1, MPI_INTEGER, worker_rank, STOP_TAG, MPI_COMM_WORLD, ierr)
                    active_workers = active_workers - 1
                end if
            end do

            do while (active_workers > 0)
                call MPI_RECV(recv_packet, 2, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, RESULT_TAG, &
                            MPI_COMM_WORLD, status, ierr)
                worker_rank = status(MPI_SOURCE)
                
                recv_count = recv_count + 1
                if (recv_count <= simu) then
                    all_results(recv_count) = recv_packet(1)
                    running_sum = running_sum + recv_packet(1)
                    total_accumulated_time = total_accumulated_time + recv_packet(2)
                end if
                
                if (mod(recv_count, print_step) == 0) then
                    print '(A,I0,A,I0,A,F10.4)', " Progress: ", recv_count, "/", simu, " Avg:", running_sum/dble(recv_count)
                end if

                if (next_task <= simu) then
                    call MPI_SEND(next_task, 1, MPI_INTEGER, worker_rank, WORK_TAG, MPI_COMM_WORLD, ierr)
                    next_task = next_task + 1
                else
                    call MPI_SEND(0, 1, MPI_INTEGER, worker_rank, STOP_TAG, MPI_COMM_WORLD, ierr)
                    active_workers = active_workers - 1
                end if
            end do

            global_arl = running_sum / dble(simu)
            print '(A,F10.4,A,F10.4,A,F10.4,A)', " >>>> Iteration Done. Result ARL =", global_arl, &
                  " Total Calc Time: ", total_accumulated_time, "s"

            if (present(fid)) then
                write(fid, '(F14.6, 2X, F14.6, 2X, F14.4)') Limit, global_arl, total_accumulated_time
            end if

            if (global_arl < IcArl) then
                LimitL = Limit
            else
                LimitR = Limit
            end if

        else
            ! Worker Logic
            do
                call MPI_RECV(task, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
                if (status(MPI_TAG) == STOP_TAG) exit
                
                call get_clock_seed(time_seed)
                seed = mod(time_seed + task*7919 + my_rank*99991, 2147483647) + 1
                call RNSET(seed)

                ! Init Sampling
                sampling_mask = 0; nodes_count = 0
                do while (nodes_count < num_samplingnodes)
                    call RNUND(num_allnodes, one_vec)
                    if (sampling_mask(one_vec(1)) == 0) then
                        nodes_count = nodes_count + 1
                        sampling_index(nodes_count) = one_vec(1)
                        sampling_mask(one_vec(1)) = 1
                    end if
                end do

                call Reset_CDS_State()
                call Set_SampleParams(ICorOC_in=IC, noise_type_in=noise_type, &
                                      back_index_in=back_index)

                ! 2. Warm-up
                do j = 1, Time_window
                    call GenerateOnlineSample(sampling_index, OnlineSample)
                    call CDS(OnlineSample, sampling_index, charting_statistic)
                end do

                ! 3. Start Timer
                t_sim_start = MPI_Wtime()

                ! 4. Monitoring
                sim_result = 0.0d0
                do while (charting_statistic < Limit)
                    call GenerateOnlineSample(sampling_index, OnlineSample)
                    call CDS(OnlineSample, sampling_index, charting_statistic)
                    sim_result = sim_result + 1.0d0
                    if (sim_result >= Maxarl) exit
                end do
                
                ! 5. Stop Timer
                t_sim_end = MPI_Wtime()
                local_sim_time = t_sim_end - t_sim_start
                
                ! 6. Send Packet
                send_packet(1) = sim_result
                send_packet(2) = local_sim_time
                call MPI_SEND(send_packet, 2, MPI_DOUBLE_PRECISION, 0, RESULT_TAG, MPI_COMM_WORLD, ierr)
            end do
        end if
        
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    end do 

    if (is_master) then
        arl = global_arl
        std = 0.0d0
        do i = 1, simu
            std = std + (all_results(i) - arl)**2
        end do
        std = sqrt(std / dble(max(1, simu-1))) / sqrt(dble(simu))
        deallocate(all_results)
    end if

    call MPI_BCAST(arl, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(std, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

end subroutine CDS_limitSearch
                           
!==============================================================
!  Subroutine: CDS_OCARL
!  Returns out-of-control, exploration, and exploitation statistics.
!==============================================================
subroutine CDS_OCARL(noise_type, back_index, anomaly_type, &
                 theta, shiftradius, shiftvalue, shiftarea, shiftellipticity, shiftradius_zero, &
                 delta_radius, delta_time, amplitude, shift_pointnum, & 
                 Limit, &
                 res_OC_ARL, res_Exploration_ARL, res_Exploitation_Num)
                 
    implicit none
    
    ! --- Inputs ---
    integer, intent(in) :: noise_type, back_index, anomaly_type
    real(dp), intent(in) :: Limit
    
    real(dp), intent(in), optional :: theta
    real(dp), intent(in), optional :: shiftradius, shiftvalue
    real(dp), intent(in), optional :: shiftarea, shiftellipticity
    real(dp), intent(in), optional :: shiftradius_zero, delta_radius, delta_time
    real(dp), intent(in), optional :: amplitude
    integer, intent(in),  optional :: shift_pointnum
    

    ! --- Outputs (Structured) ---
    type(ARL_Stats_type), intent(out) :: res_OC_ARL
    type(ARL_Stats_type), intent(out), optional :: res_Exploration_ARL
    type(ARL_Stats_type), intent(out), optional :: res_Exploitation_Num

    ! --- Locals ---
    integer :: ierr, num_procs, my_rank, status(MPI_STATUS_SIZE)
    integer :: task, next_task, active_workers, worker_rank
    integer :: i, j, seed, time_seed, recv_count, print_step
    logical :: is_master
    integer, parameter :: WORK_TAG = 1, STOP_TAG = 2, RESULT_TAG = 3
    real(dp) :: idx_vec(2)
    
    ! Data Packet: [OC_ARL, Exp_ARL, Exploit_Num]
    real(dp) :: local_packet(3), recv_packet(3)
    real(dp) :: running_sum_OC
    
    ! Arrays for Master
    real(dp), allocatable :: all_OC(:), all_Exp(:), all_Exploit(:)
    real(dp), allocatable :: random_coords(:,:)
    
    real(dp) :: charting_statistic, sim_time
    integer :: sampling_index(num_samplingnodes)
    real(dp) :: OnlineSample(num_samplingnodes)
    integer :: one_vec(1)
    integer :: sampling_mask(num_allnodes), nodes_count
    
    ! Logic variables
    integer :: current_anom_count
    logical :: exploration_found
    real(dp) :: exploration_time, total_exploited_nodes
    
    type(AnomalyParams_type) :: local_ano_params

    if (.not. perf_initialized) stop "Error: Perf Eval not initialized."

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
    is_master = (my_rank == 0)

    if (is_master) then
        allocate(all_OC(simu))
        allocate(all_Exp(simu))
        allocate(all_Exploit(simu))
        all_OC = 0.0d0; all_Exp = 0.0d0; all_Exploit = 0.0d0
        print *, "==== Begin CDS_OCARL_Extended ===="
    end if
    ! --- Fill Static Parameters ---
    local_ano_params%type_id = anomaly_type
    if (present(shiftradius))      local_ano_params%radius = shiftradius
    if (present(shiftvalue))       local_ano_params%value = shiftvalue
    if (present(shiftarea))      local_ano_params%area = shiftarea
    if (present(shiftellipticity))  local_ano_params%ellipticity = shiftellipticity            
    if (present(theta))            local_ano_params%theta = theta
    if (present(shiftradius_zero)) local_ano_params%radius = shiftradius_zero 
    if (present(delta_radius))     local_ano_params%delta_r = delta_radius
    if (present(delta_time))       local_ano_params%delta_t = delta_time
    if (present(amplitude))        local_ano_params%value = amplitude
    if (present(shift_pointnum))        local_ano_params%num_points  = shift_pointnum
    ! --- Master Distribution ---
    if (is_master) then
        next_task = 1
        recv_count = 0
        running_sum_OC = 0.0d0
        active_workers = num_procs - 1
        print_step = max(10, int(simu * 0.1)) 

        do worker_rank = 1, num_procs - 1
            if (next_task <= simu) then
                call MPI_SEND(next_task, 1, MPI_INTEGER, worker_rank, WORK_TAG, MPI_COMM_WORLD, ierr)
                next_task = next_task + 1
            else
                call MPI_SEND(0, 1, MPI_INTEGER, worker_rank, STOP_TAG, MPI_COMM_WORLD, ierr)
                active_workers = active_workers - 1
            end if
        end do

        do while (active_workers > 0)
            call MPI_RECV(recv_packet, 3, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, RESULT_TAG, &
                        MPI_COMM_WORLD, status, ierr)
            worker_rank = status(MPI_SOURCE)
            
            recv_count = recv_count + 1
            if (recv_count <= simu) then
                all_OC(recv_count)      = recv_packet(1)
                all_Exp(recv_count)     = recv_packet(2)
                all_Exploit(recv_count) = recv_packet(3)
                running_sum_OC = running_sum_OC + recv_packet(1)
            end if
            
            if (mod(recv_count, print_step) == 0) then
                print '(A,I0,A,I0,A,F10.4)', " CDS OCARL Progress: ", recv_count, "/", simu, &
                    " Avg:", running_sum_OC/dble(recv_count)
            end if
            
            if (next_task <= simu) then
                call MPI_SEND(next_task, 1, MPI_INTEGER, worker_rank, WORK_TAG, MPI_COMM_WORLD, ierr)
                next_task = next_task + 1
            else
                call MPI_SEND(0, 1, MPI_INTEGER, worker_rank, STOP_TAG, MPI_COMM_WORLD, ierr)
                active_workers = active_workers - 1
            end if
        end do

    else
        ! --- Worker Logic ---
        do
            call MPI_RECV(task, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
            if (status(MPI_TAG) == STOP_TAG) exit
            
            call get_clock_seed(time_seed)
            seed = mod(time_seed + task*7919 + my_rank*99991, 2147483647) + 1
            call RNSET(seed)

            ! 1. Init Sampling
            sampling_mask = 0; nodes_count = 0
            do while (nodes_count < num_samplingnodes)
                call RNUND(num_allnodes, one_vec)
                if (sampling_mask(one_vec(1)) == 0) then
                    nodes_count = nodes_count + 1
                    sampling_index(nodes_count) = one_vec(1)
                    sampling_mask(one_vec(1)) = 1
                end if
            end do

            ! 2. Reset
            call Reset_CDS_State()

            ! 3. Warm-up
            call Set_SampleParams(ICorOC_in=IC, noise_type_in=noise_type, &
                                  back_index_in=back_index)

            do j = 1, Time_window
                call GenerateOnlineSample(sampling_index, OnlineSample)
                call CDS(OnlineSample, sampling_index, charting_statistic)
            end do

            select case (anomaly_type)
            case (TYPE_CIRCLE, TYPE_ST, TYPE_ELLIPSE, TYPE_CRESCENT)
                CALL DRNUN (2, idx_vec)
                local_ano_params%center_idx = idx_vec 
            case (TYPE_BSPLINE)
                call RNUND(max(1,knot1_square_TS), one_vec)
                local_ano_params%bspline_idx = min(knot1_square_TS, one_vec(1))
            case (TYPE_RANDOM_POINTS_CONV)
                local_ano_params%num_points = shift_pointnum
                if (allocated(local_ano_params%points_coords)) deallocate(local_ano_params%points_coords)
                allocate(local_ano_params%points_coords(2, shift_pointnum))
                call generate_unique_random_points(shift_pointnum, random_coords)
                local_ano_params%points_coords=random_coords
            end select
        
            ! 5. Run OCARL (Monitoring)
            call Set_SampleParams(ICorOC_in=OC)
            ! 5. Monitoring
            sim_time = 0.0d0
            exploration_found = .false.
            exploration_time  = 0.0d0
            total_exploited_nodes = 0.0d0

            do while (charting_statistic < Limit)
                local_ano_params%time_idx = int(sim_time)
                call Set_SampleParams(ICorOC_in=OC, AnoParams_in=local_ano_params)
                
                ! Capture anomaly count
                call GenerateOnlineSample(sampling_index, OnlineSample, count_out=current_anom_count)
                
                ! Stats Logic
                if ((.not. exploration_found) .and. (current_anom_count > 0)) then
                    exploration_found = .true.
                    exploration_time = sim_time + 1.0d0 
                end if
                total_exploited_nodes = total_exploited_nodes + dble(current_anom_count)
                
                ! Algo Call
                call CDS(OnlineSample, sampling_index, charting_statistic)
                
                sim_time = sim_time + 1.0d0
                if (sim_time >= Maxarl) exit
            end do
            
            if (.not. exploration_found) exploration_time = sim_time
            
            local_packet(1) = sim_time
            local_packet(2) = exploration_time
            local_packet(3) = total_exploited_nodes
            
            call MPI_SEND(local_packet, 3, MPI_DOUBLE_PRECISION, 0, RESULT_TAG, MPI_COMM_WORLD, ierr)
        end do
    end if

    if (is_master) then
        call Compute_ARL_Stats(all_OC, simu, res_OC_ARL)
        if (present(res_Exploration_ARL)) call Compute_ARL_Stats(all_Exp, simu, res_Exploration_ARL)
        if (present(res_Exploitation_Num)) call Compute_ARL_Stats(all_Exploit, simu, res_Exploitation_Num)
        deallocate(all_OC, all_Exp, all_Exploit)
    end if

end subroutine CDS_OCARL
                 
        !======================================================================
! Helper Subroutine: Compute_ARL_Stats
! Computes Mean, Std, Min, Max, Q1, Median, Q3 from an array
!======================================================================
subroutine Compute_ARL_Stats(data_arr, n, stats)
    implicit none
    real(dp), intent(inout) :: data_arr(:) ! inout for sorting
    integer, intent(in) :: n
    type(ARL_Stats_type), intent(out) :: stats
    integer :: idx_arr(n), i
    
    ! Mean
    stats%mean = sum(data_arr) / dble(n)
    
    ! Std Error
    stats%std_err = 0.0d0
    do i = 1, n
        stats%std_err = stats%std_err + (data_arr(i) - stats%mean)**2.0d0
    end do
    stats%std_err = sqrt(stats%std_err / dble(max(1,n-1))) / sqrt(dble(n))
    
    ! Sort for Quantiles & Min/Max
    idx_arr = [(i, i=1,n)]
    call SVRGP(data_arr, data_arr, idx_arr) ! IMSL Sort or similar
    
    stats%min_val = data_arr(1)
    stats%max_val = data_arr(n)
    
    ! Quantiles (Simple interpolation or nearest rank)
    stats%q1     = data_arr(max(1, int(0.25d0 * n)))
    stats%median = data_arr(max(1, int(0.50d0 * n)))
    stats%q3     = data_arr(max(1, int(0.75d0 * n)))
    
end subroutine Compute_ARL_Stats         

!==============================================================
!  Subroutine: CDS_ICPerformance
!  Runs on rank 0 to preserve process state and drift.
!==============================================================
subroutine CDS_ICPerformance(noise_type, back_index, Limit, falserate)
    use, intrinsic :: iso_fortran_env, only: output_unit, real64, int64
    implicit none

    ! --- Parameters ---
    integer, parameter :: dp = real64

    ! --- Input/Output Arguments ---
    integer, intent(in) :: noise_type, back_index
    real(dp), intent(in) :: Limit
    real(dp), intent(out) :: falserate

    ! --- MPI Variables ---
    integer :: ierr, num_procs, my_rank
    ! Note: Worker/Tag variables removed as they are not needed for serial execution

    ! --- Simulation Variables ---
    integer :: j, seed, time_seed, nodes_count, print_step
    integer :: one_vec(1)
    integer, allocatable :: sampling_index(:), sampling_mask(:)
    real(dp), allocatable :: OnlineSample(:)
    real(dp) :: charting_statistic, loc_back_bw
    
    ! --- Counters (Use 64-bit for long runs) ---
    integer(int64) :: global_falsecount_long
    integer(int64) :: total_steps_long

    ! --- Check Initialization ---
    if (.not. perf_initialized) stop "Error: Perf Eval not initialized."

    ! --- MPI Setup ---
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

    ! --- Master Execution Logic (Serial) ---
    if (my_rank == 0) then
        print *, "==== Begin CDS_ICPerformance (Serial/Long-Run Mode) ===="
        print *, "Total Steps Target:", IC_TestRuns
        print *, "Control Limit:", Limit

        ! 1. Memory Allocation
        allocate(sampling_index(num_samplingnodes))
        allocate(sampling_mask(num_allnodes))
        allocate(OnlineSample(num_allnodes))

        ! 2. Parameter Setup
        ! Ensure bandwidth or other kernel params are set
        loc_back_bw = kernel_bandwidth 
        ! Note: Assuming Set_SampleParams sets up necessary module variables for CDS
        call Set_SampleParams(ICorOC_in=IC, noise_type_in=noise_type, &
                              back_index_in=back_index)

        ! 3. Random Seed Initialization
        call get_clock_seed(time_seed)
        seed = mod(time_seed + 54321, 2147483647) + 1
        call RNSET(seed)

        ! 4. Sensor/Node Selection (Fixed for the entire run)
        sampling_mask = 0
        nodes_count = 0
        do while (nodes_count < num_samplingnodes)
            call RNUND(num_allnodes, one_vec)
            if (sampling_mask(one_vec(1)) == 0) then
                nodes_count = nodes_count + 1
                sampling_index(nodes_count) = one_vec(1)
                sampling_mask(one_vec(1)) = 1
            end if
        end do

        ! 5. Simulation Start
        ! A. Reset CDS Algorithm State (Once at the beginning)
        call Reset_CDS_State()

        ! B. Warm-up (Fill window / initialize history)
        do j = 1, Time_window
            call GenerateOnlineSample(sampling_index, OnlineSample)
            call CDS(OnlineSample, sampling_index, charting_statistic)
        end do

        ! C. Continuous Monitoring Loop (The Long Run)
        global_falsecount_long = 0
        print_step = max(1000, int(IC_runs * 0.1))

        do j = 1, IC_runs
            call GenerateOnlineSample(sampling_index, OnlineSample)
            call CDS(OnlineSample, sampling_index, charting_statistic)
            
            ! Check for False Alarm
            if (charting_statistic > Limit) then
                global_falsecount_long = global_falsecount_long + 1
            end if

            ! Optional: Progress Log
            if (mod(j, print_step) == 0) then
                print '(A, I0, A, I0, A, F6.2, A)', " Step: ", j, " / ", IC_runs, " (", &
                      (real(j, dp)/real(IC_runs, dp)*100.0_dp), "%)"
            end if
        end do

        ! 6. Calculate Results
        total_steps_long = int(IC_runs, int64)
        falserate = real(global_falsecount_long, dp) / real(total_steps_long, dp)

        print '(A, I0)', " Total False Alarms: ", global_falsecount_long
        print '(A, ES12.5)', " Final IC FAR (CDS): ", falserate

        ! Cleanup
        deallocate(sampling_index, sampling_mask, OnlineSample)
    end if

    ! --- Synchronization ---
    ! Broadcast the result to all other ranks so they can exit the subroutine consistently
    call MPI_BCAST(falserate, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

end subroutine CDS_ICPerformance

end module performance_evaluation
