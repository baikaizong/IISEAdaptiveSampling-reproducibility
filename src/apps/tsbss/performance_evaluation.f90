!==============================================================
!  Module: performance_evaluation
!  Purpose:
!      Evaluate TSBSS detection limits and run-length performance.
!      Supports Limit Search, IC Performance, and OCARL.
!==============================================================
module performance_evaluation
    use GlobalSettings_mod
    use SparseMatrix_mod
    use SamplesGeneration_mod
    use TSBSS_mod            
    use utils_mod
    use mpi
    use RNSET_INT
    use RuntimeConfig_mod, only: get_clock_seed
    use RNUND_INT
    use SVRGP_INT
    implicit none
    public 
    
    ! Internal State Flags
    logical, save :: perf_initialized = .false.

    contains
    
!==============================================================
!  Subroutine: Init_Performance_Evaluation
!  Purpose: Initialize TSBSS and SampleGeneration static data.
!==============================================================
subroutine Init_Performance_Evaluation(B0_all, B1_all, CovMat)
    
    implicit none
    type(Bbasis_sparse_type), intent(in), optional :: B0_all, B1_all, CovMat

    ! 1. Initialize TSBSS Static Data (Basis Matrices)
    ! Note: TSBSS relies heavily on B0 and B1
    if (present(B1_all)) then
        call Init_TSBSS_Data(B0_in=B0_all, B1_in=B1_all)
    else
        stop "Error: B1_all is required for TSBSS initialization."
    end if

    ! 2. Initialize SamplesGeneration Static Data
    call Init_SamplesGeneration(B0_in=B0_all, B1_in=B1_all, Cov_in=CovMat)

    perf_initialized = .true.

end subroutine Init_Performance_Evaluation

!==============================================================
!  Subroutine: Clean_Performance_Evaluation
!==============================================================
subroutine Clean_Performance_Evaluation()
    implicit none
    call Clean_TSBSS()
    call Clean_SamplesGeneration()
    perf_initialized = .false.
end subroutine Clean_Performance_Evaluation

!==============================================================
!  Subroutine: TSBSS_limitSearch
!  Cumulative simulation CPU time is recorded.
!==============================================================
subroutine TSBSS_limitSearch(noise_type, back_index, &
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
    real(dp) :: running_sum
    real(dp) :: LimitL, LimitR, charting_statistic
    real(dp) :: sim_result, global_arl
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

    ! Configure IC Params
    call Set_SampleParams(ICorOC_in=IC, noise_type_in=noise_type, &
                            back_index_in=back_index)

    if (is_master) then
        LimitL = BLimitL
        LimitR = BLimitR
        global_arl = 0.0d0
        allocate(all_results(simu))
        print *, "==== Begin TSBSS_limitSearch (Cumulative Timing) ===="
        
        if (present(fid)) then
             ! Header: Limit, ARL, Total_CPU_Time
             write(fid, '(A)') "Iter_Limit,  Result_ARL,  Total_Calc_Time(s)"
        end if
    end if

    !================ Outer binary search loop =================
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

        !================ Dynamic Task Allocation =================
        if (is_master) then
            next_task = 1
            recv_count = 0
            running_sum = 0.0d0   
            total_accumulated_time = 0.0d0 ! Reset timer accumulator
            active_workers = num_procs - 1
            print_step = max(10, int(simu * 0.1)) 
            
            ! Send initial tasks
            do worker_rank = 1, num_procs - 1
                if (next_task <= simu) then
                    call MPI_SEND(next_task, 1, MPI_INTEGER, worker_rank, WORK_TAG, MPI_COMM_WORLD, ierr)
                    next_task = next_task + 1
                else
                    call MPI_SEND(0, 1, MPI_INTEGER, worker_rank, STOP_TAG, MPI_COMM_WORLD, ierr)
                    active_workers = active_workers - 1
                end if
            end do

            ! Collect Results
            do while (active_workers > 0)
                ! Receive PACKET (Size 2)
                call MPI_RECV(recv_packet, 2, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, RESULT_TAG, &
                            MPI_COMM_WORLD, status, ierr)
                worker_rank = status(MPI_SOURCE)
                
                recv_count = recv_count + 1
                if (recv_count <= simu) then
                    ! Packet(1) is ARL
                    all_results(recv_count) = recv_packet(1)
                    running_sum = running_sum + recv_packet(1)
                    
                    ! Packet(2) is Time
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
            
            ! Write to file: Limit, ARL, and Time
            if (present(fid)) then
                write(fid, '(F14.6, 2X, F14.6, 2X, F14.4)') Limit, global_arl, total_accumulated_time
            end if
            
            if (global_arl < IcArl) then
                LimitL = Limit
            else
                LimitR = Limit
            end if

        else
            ! ==================== WORKER PROCESS ====================
            do
                call MPI_RECV(task, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
                if (status(MPI_TAG) == STOP_TAG) exit
                
                call get_clock_seed(time_seed)
                seed = mod(time_seed + task*7919 + my_rank*99991, 2147483647) + 1
                call RNSET(seed)

                ! 1. Initialize Sampling
                sampling_mask = 0; nodes_count = 0
                do while (nodes_count < num_samplingnodes)
                    call RNUND(num_allnodes, one_vec)
                    if (sampling_mask(one_vec(1)) == 0) then
                        nodes_count = nodes_count + 1
                        sampling_index(nodes_count) = one_vec(1)
                        sampling_mask(one_vec(1)) = 1
                    end if
                end do

                ! 2. Reset TSBSS State
                call Reset_TSBSS_State()

                ! 4. Warm-up
                do j = 1, Time_window
                    call GenerateOnlineSample(sampling_index, OnlineSample)
                    call TSBSS(back_index, OnlineSample, sampling_index, charting_statistic)
                end do

                ! 5. Start Timer
                t_sim_start = MPI_Wtime()

                ! 6. Monitoring
                sim_result = 0.0d0
                do while (charting_statistic < Limit)
                    call GenerateOnlineSample(sampling_index, OnlineSample)
                    call TSBSS(back_index, OnlineSample, sampling_index, charting_statistic)
                    sim_result = sim_result + 1.0d0
                    if (sim_result >= Maxarl) exit
                end do
                
                ! 7. Stop Timer
                t_sim_end = MPI_Wtime()
                local_sim_time = t_sim_end - t_sim_start

                ! 8. Pack and Send
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

end subroutine TSBSS_limitSearch

!==============================================================
!  Subroutine: TSBSS_OCARL
!  Returns out-of-control, exploration, and exploitation statistics.
!==============================================================
subroutine TSBSS_OCARL(noise_type, back_index, anomaly_type, &
                 theta, shiftradius, shiftvalue, shiftarea, shiftellipticity, shiftradius_zero, &
                 delta_radius, delta_time, amplitude, &
                 Limit, &
                 res_OC_ARL, res_Exploration_ARL, res_Exploitation_Num)
                 
    implicit none
    
    ! --- Inputs ---
    integer, intent(in) :: noise_type, back_index, anomaly_type
    real(dp), intent(in) :: Limit
    
    real(dp), intent(in), optional :: theta, shiftradius, shiftvalue
    real(dp), intent(in), optional :: shiftarea, shiftellipticity
    real(dp), intent(in), optional :: shiftradius_zero, delta_radius, delta_time
    real(dp), intent(in), optional :: amplitude

    ! --- Outputs ---
    type(ARL_Stats_type), intent(out) :: res_OC_ARL
    type(ARL_Stats_type), intent(out), optional :: res_Exploration_ARL
    type(ARL_Stats_type), intent(out), optional :: res_Exploitation_Num

    ! --- Locals ---
    integer :: ierr, num_procs, my_rank, status(MPI_STATUS_SIZE)
    integer :: task, next_task, active_workers, worker_rank
    integer :: i, j, seed, time_seed, recv_count, print_step
    logical :: is_master
    integer, parameter :: WORK_TAG = 1, STOP_TAG = 2, RESULT_TAG = 3
    
    real(dp) :: local_packet(3), recv_packet(3), running_sum_OC
    real(dp), allocatable :: all_OC(:), all_Exp(:), all_Exploit(:)
    
    real(dp) :: charting_statistic, sim_time
    integer :: sampling_index(num_samplingnodes)
    real(dp) :: OnlineSample(num_samplingnodes)
    integer :: one_vec(1)
    integer :: sampling_mask(num_allnodes), nodes_count
    real(dp) :: Idx_vec(2)
    
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
        print *, "==== Begin TSBSS_OCARL_Extended ===="
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
    call Set_SampleParams( noise_type_in=noise_type,back_index_in=back_index)    

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
                print '(A,I0,A,I0,A,F10.4)', " TSBSS OCARL Progress: ", recv_count, "/", simu, &
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

            ! 2. Reset Stats
            call Reset_TSBSS_State()

            ! 3. Warm-up
            call Set_SampleParams(ICorOC_in=IC)
            do j = 1, Time_window
                call GenerateOnlineSample(sampling_index, OnlineSample)
                call TSBSS(back_index, OnlineSample, sampling_index, charting_statistic)
            end do

            ! 4. Randomize Location
            select case (anomaly_type)
            case (TYPE_CIRCLE, TYPE_ST, TYPE_CIRCLE_CONV, TYPE_ELLIPSE, TYPE_CRESCENT)
                CALL DRNUN (2, idx_vec)
                local_ano_params%center_idx = idx_vec 
            case (TYPE_BSPLINE)
                call RNUND(max(1,knot1_square_TS), one_vec)
                local_ano_params%bspline_idx = min(knot1_square_TS, one_vec(1))
            end select

            ! 5. Monitoring
            call Set_SampleParams(ICorOC_in=OC) 
            
            sim_time = 0.0d0
            exploration_found = .false.
            exploration_time  = 0.0d0
            total_exploited_nodes = 0.0d0

            do while (charting_statistic < Limit)
                local_ano_params%time_idx = int(sim_time)
                call Set_SampleParams(ICorOC_in=OC, AnoParams_in=local_ano_params)
                
                call GenerateOnlineSample(sampling_index, OnlineSample, count_out=current_anom_count)
                
                ! Stats
                if ((.not. exploration_found) .and. (current_anom_count > 0)) then
                    exploration_found = .true.
                    exploration_time = sim_time + 1.0d0 
                end if
                total_exploited_nodes = total_exploited_nodes + dble(current_anom_count)
                
                call TSBSS(back_index, OnlineSample, sampling_index, charting_statistic)
                
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

end subroutine TSBSS_OCARL
                 
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
!  Subroutine: TSBSS_ICPerformance
!  Evaluate in-control performance for TSBSS.
!==============================================================
subroutine TSBSS_ICPerformance(noise_type, back_index, Limit, falserate)
    implicit none

    integer, intent(in) :: noise_type, back_index
    real(dp), intent(in) :: Limit
    real(dp), intent(out) :: falserate

    integer :: ierr, num_procs, my_rank, status(MPI_STATUS_SIZE)
    integer :: run, next_run, active_workers, worker_rank
    logical :: is_master
    integer, parameter :: WORK_TAG = 1, STOP_TAG = 2, RESULT_TAG = 3
    integer :: false_result(2)
    
    integer :: i, j, seed, time_seed, local_falsecount, completed_runs, print_step 
    real(dp) :: charting_statistic
    integer  :: sampling_index(num_samplingnodes)
    real(dp) :: OnlineSample(num_samplingnodes)
    integer :: global_falsecount, one_vec(1)
    integer :: sampling_mask(num_allnodes), nodes_count


    if (.not. perf_initialized) stop "Error: Perf Eval not initialized."

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
    is_master = (my_rank == 0)



    if (is_master) then
        print *, "==== Begin TSBSS_ICPerformance_dynamic ===="
        global_falsecount = 0
    end if
    call Set_SampleParams(ICorOC_in=IC, noise_type_in=noise_type, &
                            back_index_in=back_index)

    if (is_master) then
        next_run = 1
        active_workers = num_procs - 1
        completed_runs = 0         
        print_step = max(100, int(IC_runs * 0.05)) 
        
        do worker_rank = 1, num_procs - 1
            if (next_run <= IC_runs) then
                call MPI_SEND(next_run, 1, MPI_INTEGER, worker_rank, WORK_TAG, MPI_COMM_WORLD, ierr)
                next_run = next_run + 1
            else
                call MPI_SEND(0, 1, MPI_INTEGER, worker_rank, STOP_TAG, MPI_COMM_WORLD, ierr)
                active_workers = active_workers - 1
            end if
        end do

        do while (active_workers > 0)
            call MPI_RECV(false_result, 2, MPI_INTEGER, MPI_ANY_SOURCE, RESULT_TAG, &
                        MPI_COMM_WORLD, status, ierr)
            worker_rank = status(MPI_SOURCE)
            
            global_falsecount = global_falsecount + false_result(1)
            completed_runs = completed_runs + 1 
            
            if (mod(completed_runs, print_step) == 0) then
                print '(A,I0,A,F8.5)', " IC Checked: ", completed_runs, " FAR:", dble(global_falsecount)/dble(completed_runs)
            end if
            
            if (next_run <= IC_runs) then
                call MPI_SEND(next_run, 1, MPI_INTEGER, worker_rank, WORK_TAG, MPI_COMM_WORLD, ierr)
                next_run = next_run + 1
            else
                call MPI_SEND(0, 1, MPI_INTEGER, worker_rank, STOP_TAG, MPI_COMM_WORLD, ierr)
                active_workers = active_workers - 1
            end if
        end do

        falserate = dble(global_falsecount) / dble(IC_runs)
        print *, "Final IC FAR =", falserate

    else
        do
            call MPI_RECV(run, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
            if (status(MPI_TAG) == STOP_TAG) exit
            
            call get_clock_seed(time_seed)
            seed = mod(time_seed + run*7919 + my_rank*99991, 2147483647) + 1
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
            call Reset_TSBSS_State()

            ! 4. Warm-up
            do j = 1, Time_window
                call GenerateOnlineSample(sampling_index, OnlineSample)
                call TSBSS(back_index, OnlineSample, sampling_index, charting_statistic)
            end do
            
            ! 5. Check
            call GenerateOnlineSample(sampling_index, OnlineSample)
            call TSBSS(back_index, OnlineSample, sampling_index, charting_statistic)
            
            local_falsecount = 0
            if (charting_statistic > Limit) local_falsecount = 1
            
            false_result(1) = local_falsecount
            false_result(2) = 1
            call MPI_SEND(false_result, 2, MPI_INTEGER, 0, RESULT_TAG, MPI_COMM_WORLD, ierr)
        end do
    end if

    call MPI_BCAST(falserate, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

end subroutine TSBSS_ICPerformance

!==============================================================
!  Subroutine: TSBSS_Estimate_Limit_Range
!  Purpose: 
!      Estimates Limit range for TSBSS.
!      Since TSBSS statistics do not depend on the limit itself,
!      only a SINGLE run is required (no iteration/convergence needed).
!==============================================================
subroutine TSBSS_EstimateLimitRange(noise_type, back_index, &
                                LimitL, LimitR)
    implicit none

    !---------------- Inputs ----------------
    integer, intent(in) :: noise_type, back_index

    !---------------- Outputs ----------------
    real(dp), intent(out) :: LimitL, LimitR

    !---------------- MPI Variables -----------------
    integer :: ierr, num_procs, my_rank, status(MPI_STATUS_SIZE)
    integer :: task, next_task, active_workers, worker_rank
    integer :: recv_count
    logical :: is_master
    integer, parameter :: WORK_TAG = 1, STOP_TAG = 2, RESULT_TAG = 3

    !---------------- Locals -----------------
    integer :: i, j, seed, time_seed
    real(dp) :: local_max_stat, sim_max_stat
    real(dp), allocatable :: all_max_stats(:)
    
    integer :: sampling_index(num_samplingnodes)
    real(dp) :: OnlineSample(num_samplingnodes)
    real(dp) :: charting_statistic
    
    integer :: sampling_mask(num_allnodes), nodes_count, one_vec(1)
    integer :: target_idx_L, target_idx_U
    real(dp) :: target_prob_L, target_prob_U

    ! Safety Check
    if (.not. perf_initialized) stop "Error: Perf Eval not initialized."

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
    is_master = (my_rank == 0)



    ! Initialize Master Buffers
    if (is_master) then
        allocate(all_max_stats(IC_TestRuns))
        print *, "==== Auto-Estimating Limit Range (TSBSS Single Pass) ===="
        print *, "Target Runs:", IC_TestRuns, " Target ARL:", IcArl
    end if

    call Set_SampleParams(ICorOC_in=IC, noise_type_in=noise_type, &
                            back_index_in=back_index)

    ! 2. Execute Parallel Simulations
    if (is_master) then
        next_task = 1
        recv_count = 0
        active_workers = num_procs - 1
        
        ! Dispatch
        do worker_rank = 1, num_procs - 1
            if (next_task <= IC_TestRuns) then
                call MPI_SEND(next_task, 1, MPI_INTEGER, worker_rank, WORK_TAG, MPI_COMM_WORLD, ierr)
                next_task = next_task + 1
            else
                call MPI_SEND(0, 1, MPI_INTEGER, worker_rank, STOP_TAG, MPI_COMM_WORLD, ierr)
                active_workers = active_workers - 1
            end if
        end do

        ! Collect
        do while (active_workers > 0)
            call MPI_RECV(local_max_stat, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, RESULT_TAG, &
                        MPI_COMM_WORLD, status, ierr)
            worker_rank = status(MPI_SOURCE)
            
            recv_count = recv_count + 1
            if (recv_count <= IC_TestRuns) then
                all_max_stats(recv_count) = local_max_stat
            end if
            
            if (next_task <= IC_TestRuns) then
                call MPI_SEND(next_task, 1, MPI_INTEGER, worker_rank, WORK_TAG, MPI_COMM_WORLD, ierr)
                next_task = next_task + 1
            else
                call MPI_SEND(0, 1, MPI_INTEGER, worker_rank, STOP_TAG, MPI_COMM_WORLD, ierr)
                active_workers = active_workers - 1
            end if
        end do
        
    else
        ! Worker Logic
        do
            call MPI_RECV(task, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
            if (status(MPI_TAG) == STOP_TAG) exit
            
            ! Random Seed
            call get_clock_seed(time_seed)
            seed = mod(time_seed + task*313 + my_rank*997, 2147483647) + 1
            call RNSET(seed)

            ! Init Sampling
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

            ! --- TSBSS Specific Reset ---
            call Reset_TSBSS_State()
            

            do j = 1, Time_window
                call GenerateOnlineSample(sampling_index, OnlineSample)
                call TSBSS(back_index, OnlineSample, sampling_index, charting_statistic)
            end do
            sim_max_stat = -1.0d20
            do j = 1, Time_window
                call GenerateOnlineSample(sampling_index, OnlineSample)
                call TSBSS(back_index, OnlineSample, sampling_index, charting_statistic)
                if (charting_statistic > sim_max_stat) sim_max_stat = charting_statistic
            end do
            
            call MPI_SEND(sim_max_stat, 1, MPI_DOUBLE_PRECISION, 0, RESULT_TAG, MPI_COMM_WORLD, ierr)
        end do
    end if
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    ! 3. Process Results (Master Only)
    if (is_master) then
        target_prob_L = (1.0d0*Time_window/(IcArl-30.0d0))
        target_idx_L = int(dble(IC_TestRuns) * target_prob_L)
        target_idx_L = max(1, min(IC_TestRuns, target_idx_L))
        
        target_prob_U = (1.0d0*Time_window/(IcArl+60.0d0))   
        target_idx_U = int(dble(IC_TestRuns) * target_prob_U)
        target_idx_U = max(1, min(IC_TestRuns, target_idx_U))
        
        ! Find LimitL & LimitR
        call partial_quickselect(values=all_max_stats, num_all=IC_TestRuns, &
                                 top_k=IC_TestRuns-target_idx_L+1, &
                                 order_index=-1, values_sub=all_max_stats)
        LimitL = all_max_stats(target_idx_L)
        LimitR = all_max_stats(target_idx_U) 

        print '(A,F10.4,A,F10.4)', " Estimated Range: [", LimitL, ",", LimitR, "]"
    end if
    
    ! Ensure all procs have final limits
    call MPI_BCAST(LimitL, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(LimitR, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    if (is_master) deallocate(all_max_stats)

end subroutine TSBSS_EstimateLimitRange
                                
end module performance_evaluation
