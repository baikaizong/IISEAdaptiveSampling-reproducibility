!==============================================================
!  Module: performance_evaluation
!  Purpose:
!      Evaluation routines supporting all background,
!      noise, and anomaly types defined in SamplesGeneration_mod.
!==============================================================
module performance_evaluation
    use GlobalSettings_mod
    use SparseMatrix_mod
    use SamplesGeneration_mod
    use mMSTD_mod
    use Exploration_mod
    use utils_mod
    use mpi
    use RNSET_INT
    use RuntimeConfig_mod, only: get_clock_seed
    use RNUND_INT
    use SVRGP_INT
    USE RNUN_INT
    implicit none
    public 
    
    ! Internal State Flags
    logical, save :: perf_initialized = .false.

    contains
    
!==============================================================
!  Subroutine: Init_Performance_Evaluation
!==============================================================
subroutine Init_Performance_Evaluation( &
    NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId, &
    Nodesset, B0_all, B1_all, CovMat)
    
    implicit none
    real(dp), intent(in) :: NeighborDis(:,:), Kernel_NeighborDis(:,:)
    integer, intent(in) :: NeighborId(:,:), Kernel_NeighborId(:,:)
    real(dp), intent(in), optional :: Nodesset(:,:)
    type(Bbasis_sparse_type), intent(in), optional :: B0_all, B1_all, CovMat

    call Init_mMSTD_Data(NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId)
    call Init_SamplesGeneration(Nodesset, B0_all, B1_all, CovMat)

    perf_initialized = .true.

end subroutine Init_Performance_Evaluation

!==============================================================
!  Subroutine: Clean_Performance_Evaluation
!==============================================================
subroutine Clean_Performance_Evaluation()
    implicit none
    call Clean_mMSTD()
    call Clean_SamplesGeneration()
    perf_initialized = .false.
end subroutine Clean_Performance_Evaluation

     
!==============================================================
!  Subroutine: mMSTD_limitSearch
!  Records cumulative simulation CPU time independent of parallel speedup.
!==============================================================
subroutine mMSTD_limitSearch(noise_type, back_index, &
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

    !---------------- Other Locals -----------------
    integer :: time_seed, seed, i, j, print_step
    real(dp) :: running_sum
    real(dp) :: LimitL, LimitR, charting_statistic
    real(dp) :: sim_result, global_arl
    real(dp), allocatable :: all_results(:)
    
    integer :: sampling_index(num_samplingnodes)
    real(dp) :: alpha_value, OnlineSample(num_samplingnodes)
    integer :: sampling_mask(num_allnodes), nodes_count, one_vec(1)

    !---------------- Timing variables ----------------
    real(dp) :: t_sim_start, t_sim_end, local_sim_time
    real(dp) :: total_accumulated_time ! Sum of all simulation times
    
    ! Packet for sending [ARL_Result, Time_Taken]
    real(dp) :: send_packet(2), recv_packet(2) 

    if (.not. perf_initialized) stop "Error: Performance Evaluation not initialized."

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
    is_master = (my_rank == 0)

    if (is_master) then
        LimitL = BLimitL
        LimitR = BLimitR
        global_arl = 0.0d0
        allocate(all_results(simu))
        print *, "==== Begin mMKSTD limit search (cumulative timing) ===="
        if (present(fid)) then
             ! Header: Limit, ARL, Total_CPU_Time
             write(fid, '(A)') "Iter_Limit,  Result_ARL,  Total_Calc_Time(s)"
        end if
    end if

    call Set_SampleParams(ICorOC_in=IC, noise_type_in=noise_type, &
                            back_index_in=back_index)

    !================ Outer binary search loop =================
    do         
        ! 1. Check Convergence
        if (is_master) then        
            if (abs(global_arl - IcArl) < Icstd) then
                is_converged = .true.
            else
                is_converged = .false.
                Limit = 0.5d0 * (LimitL + LimitR)
                alpha_value = alpha_mM
                if (.not. present(fid)) then
                    print '(A,F12.6,A,F12.6,A,F12.6)', " Search: [", LimitL, ",", LimitR, "] Try:", Limit
                end if
            end if
        end if
        
        call MPI_BCAST(is_converged, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        if (is_converged) exit
        call MPI_BCAST(Limit, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(alpha_value, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        !================ Dynamic Task Allocation =================
        if (is_master) then
            next_task = 1
            recv_count = 0
            running_sum = 0.0d0 
            total_accumulated_time = 0.0d0 ! Reset timer accumulator for this iteration
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

            ! Collect results
            do while (active_workers > 0)
                ! Receive PACKET (Size 2) instead of single value
                call MPI_RECV(recv_packet, 2, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, RESULT_TAG, &
                            MPI_COMM_WORLD, status, ierr)
                worker_rank = status(MPI_SOURCE)
                
                recv_count = recv_count + 1
                if (recv_count <= simu) then
                    ! Packet(1) is the Run Length (sim_result)
                    all_results(recv_count) = recv_packet(1)
                    running_sum = running_sum + recv_packet(1)
                    
                    ! Packet(2) is the Simulation Time for that specific run
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

            ! Write to file: Limit, ARL, and the SUM of all simulation times
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

                ! 1. Initialize
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

                call Reset_mMSTD_State()

                ! 2. Warm-up
                do j = 1, Time_window
                    call GenerateOnlineSample(sampling_index, OnlineSample)
                    call mMSTDII(alpha_value, OnlineSample, sampling_index, charting_statistic)
                end do

                ! 3. Start Timer for this simulation
                ! Note: Warm-up is EXCLUDED as per requirement. 
                t_sim_start = MPI_Wtime()

                ! 4. Monitoring Loop
                sim_result = 0.0d0
                do while (charting_statistic < Limit)
                    call GenerateOnlineSample(sampling_index, OnlineSample)
                    call mMSTDII(alpha_value, OnlineSample, sampling_index, charting_statistic)
                    sim_result = sim_result + 1.0d0
                    if (sim_result >= Maxarl) exit
                end do
                
                ! 5. Stop Timer
                t_sim_end = MPI_Wtime()
                local_sim_time = t_sim_end - t_sim_start

                ! 6. Pack and Send (Result + Time)
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

end subroutine mMSTD_limitSearch                                 
                                 
!==============================================================
!  Subroutine: mMSTD_OCARL
!  Returns out-of-control, exploration, and exploitation statistics.
!==============================================================
subroutine mMSTD_OCARL(noise_type, back_index, anomaly_type, &
                 theta, shiftradius, shiftvalue, shiftarea, shiftellipticity, shiftradius_zero, &
                 delta_radius, delta_time, amplitude, shift_pointnum, &
                 Limit, &
                 res_OC_ARL, res_Exploration_ARL, res_Exploitation_Num)
                 
    implicit none
    
    ! --- Inputs ---
    integer, intent(in) :: noise_type, back_index, anomaly_type
    real(dp), intent(in) :: Limit
    
    ! Optional Anomaly Params
    real(dp), intent(in), optional :: theta
    real(dp), intent(in), optional :: shiftradius, shiftvalue
    real(dp), intent(in), optional :: shiftarea, shiftellipticity
    real(dp), intent(in), optional :: shiftradius_zero, delta_radius, delta_time
    real(dp), intent(in), optional :: amplitude
    integer, intent(in),  optional :: shift_pointnum

    ! --- Outputs (Structured & Optional) ---
    type(ARL_Stats_type), intent(out) :: res_OC_ARL
    type(ARL_Stats_type), intent(out), optional :: res_Exploration_ARL
    type(ARL_Stats_type), intent(out), optional :: res_Exploitation_Num

    ! --- Locals ---
    integer :: ierr, num_procs, my_rank, status(MPI_STATUS_SIZE)
    integer :: task, next_task, active_workers, worker_rank
    integer :: i, j, seed, time_seed, recv_count, print_step
    logical :: is_master
    integer, parameter :: WORK_TAG = 1, STOP_TAG = 2, RESULT_TAG = 3
    
    ! Data Packet for MPI Communication
    ! 1: Run Length (OC ARL)
    ! 2: First Detection Time (Exploration)
    ! 3: Total Anomalous Nodes Detected (Exploitation)
    real(dp) :: local_packet(3), recv_packet(3)
    
    real(dp) :: running_sum_OC, running_sum_Exp, running_sum_Exploit
    real(dp) :: charting_statistic, sim_time, alpha_value
    
    ! Storage for all simulation results to calc quantiles
    real(dp), allocatable :: all_OC(:), all_Exp(:), all_Exploit(:)
    real(dp), allocatable :: random_coords(:,:)
    
    integer :: sampling_index(num_samplingnodes)
    real(dp) :: OnlineSample(num_samplingnodes)
    real(dp) :: idx_vec(2)
    integer :: one_vec(1), simu_index(simu)
    integer :: sampling_mask(num_allnodes), nodes_count
    
    ! New variables for logic
    integer :: current_anom_count
    logical :: exploration_found
    real(dp) :: exploration_time
    real(dp) :: total_exploited_nodes
    
    type(AnomalyParams_type) :: local_ano_params
    
    if (.not. perf_initialized) stop "Error: Perf Eval not initialized."

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
    is_master = (my_rank == 0)

    if (is_master) then
        allocate(all_OC(simu))
        allocate(all_Exp(simu))
        allocate(all_Exploit(simu))
        all_OC = 0.0d0
        all_Exp = 0.0d0
        all_Exploit = 0.0d0
        print *, "==== Begin mMKSTD out-of-control evaluation ===="
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
    
    call Set_SampleParams( noise_type_in=noise_type, back_index_in=back_index)
    
    if (is_master) alpha_value = alpha_mM 
    call MPI_BCAST(alpha_value, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    if (is_master) then
        next_task = 1
        recv_count = 0
        running_sum_OC = 0.0d0
        active_workers = num_procs - 1
        print_step = max(10, int(simu * 0.1)) 

        ! Initial Dispatch
        do worker_rank = 1, num_procs - 1
            if (next_task <= simu) then
                call MPI_SEND(next_task, 1, MPI_INTEGER, worker_rank, WORK_TAG, MPI_COMM_WORLD, ierr)
                next_task = next_task + 1
            else
                call MPI_SEND(0, 1, MPI_INTEGER, worker_rank, STOP_TAG, MPI_COMM_WORLD, ierr)
                active_workers = active_workers - 1
            end if
        end do

        ! Collection Loop
        do while (active_workers > 0)
            ! Receive Packet of size 3
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
                print '(A,I0,A,I0,A,F10.4)', " Progress: ", recv_count, "/", simu, " Avg ARL:", running_sum_OC/dble(recv_count)
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
        ! Worker Logic
        do
            call MPI_RECV(task, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
            if (status(MPI_TAG) == STOP_TAG) exit
            
            call get_clock_seed(time_seed)
            seed = mod(time_seed + task*7919 + my_rank*99991, 2147483647) + 1
            call RNSET(seed)

            ! 1. Init Sampling
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

            ! 2. Reset Stats
            call Reset_mMSTD_State()

            ! 3. Warm-up (IC)
            call Set_SampleParams(ICorOC_in=IC)
            do j = 1, Time_window
                call GenerateOnlineSample(sampling_index, OnlineSample) 
                call mMSTDII(alpha_value, OnlineSample, sampling_index, charting_statistic)
            end do

            ! 4. Randomize Location
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
            
            sim_time = 0.0d0
            exploration_found = .false.
            exploration_time  = 0.0d0 ! Default if never found (should be handled)
            total_exploited_nodes = 0.0d0
            
            do while (charting_statistic < Limit)
                ! Update Dynamic params
                local_ano_params%time_idx = int(sim_time)
                call Set_SampleParams(AnoParams_in=local_ano_params)
                
                ! Call GenerateOnlineSample with count_out
                call GenerateOnlineSample(sampling_index, OnlineSample, count_out=current_anom_count)
                
                ! --- Logic for Exploration ARL ---
                ! First time we see ANY anomalous node
                if ((.not. exploration_found) .and. (current_anom_count > 0)) then
                    exploration_found = .true.
                    ! Record time. Since sim_time starts at 0, the first step is 1. 
                    ! Adjust based on your time definition. Usually step index.
                    exploration_time = sim_time + 1.0d0 
                end if
                
                ! --- Logic for Exploitation Num ---
                ! Accumulate count of anomalous nodes found
                total_exploited_nodes = total_exploited_nodes + dble(current_anom_count)
                
                call mMSTDII(alpha_value, OnlineSample, sampling_index, charting_statistic)
                
                sim_time = sim_time + 1.0d0
                if (sim_time >= Maxarl) exit
            end do
            
            ! If exploration never happened (e.g. max_arl reached or weak anomaly),
            ! set exploration time to max_arl or sim_time.
            if (.not. exploration_found) exploration_time = sim_time
            
            ! Pack Results
            local_packet(1) = sim_time
            local_packet(2) = exploration_time
            local_packet(3) = total_exploited_nodes
            
            call MPI_SEND(local_packet, 3, MPI_DOUBLE_PRECISION, 0, RESULT_TAG, MPI_COMM_WORLD, ierr)
        end do
    end if

    ! --- Master: Calculate Statistics for Structure ---
    if (is_master) then
        ! 1. OC ARL Stats
        call Compute_ARL_Stats(all_OC, simu, res_OC_ARL)
        
        ! 2. Exploration ARL Stats (Optional)
        if (present(res_Exploration_ARL)) then
            call Compute_ARL_Stats(all_Exp, simu, res_Exploration_ARL)
        end if
        
        ! 3. Exploitation Num Stats (Optional)
        if (present(res_Exploitation_Num)) then
            call Compute_ARL_Stats(all_Exploit, simu, res_Exploitation_Num)
        end if
        
        deallocate(all_OC, all_Exp, all_Exploit)
    end if


end subroutine mMSTD_OCARL


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
!  Subroutine: mMSTD_ICPerformance
!  Supports background-bandwidth configuration.
!==============================================================
subroutine mMSTD_ICPerformance(noise_type, back_index, Limit, falserate)
    use, intrinsic :: iso_fortran_env, only: output_unit, real64
    implicit none

    ! --- Arguments ---
    integer, intent(in) :: noise_type, back_index
    real(dp), intent(in) :: Limit
    real(dp), intent(out) :: falserate

    ! --- MPI Variables ---
    integer :: ierr, num_procs, my_rank
    
    ! --- Simulation Variables (Only used by Master) ---
    integer :: j, nodes_count, seed, time_seed, print_step
    integer :: one_vec(1)
    integer, allocatable :: sampling_index(:), sampling_mask(:)
    real(dp), allocatable :: OnlineSample(:)
    real(dp) :: alpha_value, charting_statistic
    integer(8) :: global_falsecount_long ! 64-bit for large counts

    ! --- Check Init ---
    if (.not. perf_initialized) stop "Error: Perf Eval not initialized."

    ! --- MPI Setup ---
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

    ! --- Master Execution Logic ---
    if (my_rank == 0) then
        
        print *, "==== mMKSTD in-control evaluation: serial mode ===="
        print *, "Running on Master Node only. Total Steps:", IC_runs

        ! 1. Allocations
        allocate(sampling_index(num_samplingnodes))
        allocate(sampling_mask(num_allnodes))
        allocate(OnlineSample(num_allnodes))

        ! 2. Setup Params
        call Set_SampleParams(ICorOC_in=IC, noise_type_in=noise_type, back_index_in=back_index)
        alpha_value = alpha_mM
        
        ! 3. Seeding
        call get_clock_seed(time_seed)
        seed = mod(time_seed + 12345, 2147483647) + 1
        call RNSET(seed)

        ! 4. Sensor Selection
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

        ! 5. Simulation Loop
        call Reset_mMSTD_State()

        ! Warm-up
        do j = 1, Time_window
            call GenerateOnlineSample(sampling_index, OnlineSample)
            call mMSTD(alpha_value, OnlineSample, sampling_index, charting_statistic)
        end do

        ! Continuous Run
        global_falsecount_long = 0
        print_step = max(1000, int(IC_runs * 0.1))

        do j = 1, IC_runs
            call GenerateOnlineSample(sampling_index, OnlineSample)
            call mMSTD(alpha_value, OnlineSample, sampling_index, charting_statistic)
            
            if (charting_statistic > Limit) then
                global_falsecount_long = global_falsecount_long + 1
            end if

            ! Optional: Progress log for very long serial runs
            if (mod(j, print_step) == 0) then
                print '(A, I0, A, I0, A, F6.2, A)', " Step: ", j, " / ", IC_runs, " (", &
                      (real(j, dp)/real(IC_runs, dp)*100.0_dp), "%)"
            end if
        end do

        ! Calculate Rate
        falserate = real(global_falsecount_long, dp) / real(IC_runs, dp)
        
        print '(A, ES12.5)', " Final IC FAR (Serial): ", falserate

        deallocate(sampling_index, sampling_mask, OnlineSample)
    end if

    ! --- Synchronization ---
    ! All other ranks wait here until Master finishes and broadcasts the result
    call MPI_BCAST(falserate, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

end subroutine mMSTD_ICPerformance

!==============================================================
!  Subroutine: mMSTD_EstimateLimitRange
!  Purpose: 
!      Estimates Limit range based on the instantaneous distribution 
!      of statistics gathered from a long IC run.
!==============================================================
subroutine mMSTD_EstimateLimitRange(noise_type, back_index, &
                                    Limit_Init, LimitL, LimitR, &
                                    back_bandwidth)
    implicit none

    !---------------- Inputs ----------------
    integer, intent(in) :: noise_type, back_index
    real(dp), intent(inout) :: Limit_Init 
    real(dp), intent(in), optional :: back_bandwidth

    !---------------- Outputs ----------------
    real(dp), intent(out) :: LimitL, LimitR

    !---------------- MPI & Load Balance -----------------
    integer :: ierr, num_procs, my_rank
    integer :: runs_per_proc, remainder
    integer :: my_run_count
    
    ! Buffers for MPI_GATHERV
    integer, allocatable :: recv_counts(:), displs(:)
    real(dp), allocatable :: local_stats(:)
    real(dp), allocatable :: global_stats(:)

    !---------------- Simulation Locals -----------------
    integer :: i, j, iter, max_iter
    integer :: sampling_index(num_samplingnodes)
    real(dp) :: OnlineSample(num_samplingnodes)
    
    ! mMKSTD state variables
    real(dp) :: ProbabilityTerm(num_allnodes), KernelVarianceTerm(num_allnodes)
    real(dp) :: DistanceTerm(num_allnodes), stdTerm(num_allnodes)
    real(dp) :: charting_statistic
    
    integer :: sampling_mask(num_allnodes), nodes_count, one_vec(1)
    integer :: seed, time_seed, idx_L, idx_U
    real(dp) :: alpha_value
    logical :: range_accepted

    ! Quantile indices
    integer :: idx_lower, idx_upper
    real(dp) :: prob_lower, prob_upper

    ! Safety Check
    if (.not. perf_initialized) stop "Error: Perf Eval not initialized."

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

    call Set_SampleParams(ICorOC_in=IC, noise_type_in=noise_type, &
                            back_index_in=back_index)
    !------------------------------------------------------------
    ! 1. Static Task Allocation (Calculate my_run_count)
    !------------------------------------------------------------
    runs_per_proc = IC_TestRuns / num_procs
    remainder = mod(IC_TestRuns, num_procs)
    
    if (my_rank < remainder) then
        my_run_count = runs_per_proc + 1
    else
        my_run_count = runs_per_proc
    end if
    
    ! Allocate local buffer
    if (my_run_count > 0) allocate(local_stats(my_run_count))

    ! Master prepares for gathering
    if (my_rank == 0) then
        allocate(global_stats(IC_TestRuns))
        allocate(recv_counts(num_procs))
        allocate(displs(num_procs))
        
        ! Pre-calculate Gatherv parameters
        do i = 0, num_procs - 1
            if (i < remainder) then
                recv_counts(i+1) = runs_per_proc + 1
            else
                recv_counts(i+1) = runs_per_proc
            end if
            
            if (i == 0) then
                displs(i+1) = 0
            else
                displs(i+1) = displs(i) + recv_counts(i)
            end if
        end do
        
        print *, "==== Auto-Estimating Limit Range (Static & Continuous) ===="
        print *, "Total Samples:", IC_TestRuns, " Target ARL:", IcArl
    end if

    max_iter = 10 
    range_accepted = .false.

    !================ Iteration Loop =================
    do iter = 1, max_iter
        
        ! A. Set the sampling-weight parameter.
        alpha_value = alpha_mM
        
        call MPI_BCAST(alpha_value, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        ! B. Run Continuous Simulation (Local)
        if (my_run_count > 0) then
            call get_clock_seed(time_seed)
            seed = mod(time_seed + my_rank*997 + iter*1009, 2147483647) + 1
            call RNSET(seed)

            ! 1. Initialize Sampling Nodes
            sampling_mask = 0; nodes_count = 0
            do while (nodes_count < num_samplingnodes)
                call RNUND(num_allnodes, one_vec)
                if (sampling_mask(one_vec(1)) == 0) then
                    nodes_count = nodes_count + 1
                    sampling_index(nodes_count) = one_vec(1)
                    sampling_mask(one_vec(1)) = 1
                end if
            end do

            ! 2. Reset Algorithm State
            call Reset_mMSTD_State()

            ! 3. Warm-up Phase (Run for Time_window steps, ignore results)
            do j = 1, Time_window
                call GenerateOnlineSample(sampling_index, OnlineSample)
                call mMSTD(alpha_value, OnlineSample, sampling_index, charting_statistic)
            end do

            ! 4. Collection Phase (Run continuously for my_run_count steps)
            do j = 1, my_run_count
                call GenerateOnlineSample(sampling_index, OnlineSample)
                call mMSTD(alpha_value, OnlineSample, sampling_index, charting_statistic)
                local_stats(j) = charting_statistic
            end do
        end if
        
        ! C. Gather Results to Master
        ! Using MPI_GATHERV because load might differ by 1 sample
        call MPI_GATHERV(local_stats, my_run_count, MPI_DOUBLE_PRECISION, &
                         global_stats, recv_counts, displs, MPI_DOUBLE_PRECISION, &
                         0, MPI_COMM_WORLD, ierr)

        ! D. Process Results (Master Only)
        ! 3. Process Results (Master Only)
        if (my_rank == 0) then
            idx_L = ceiling(dble(IC_TestRuns) * (1.0d0 / (IcArl - down_gap)))
            idx_L = max(1, min(IC_TestRuns, idx_L))
            
            ! Index for LimitR, the cutoff for the tighter upper range.
            idx_U = ceiling(dble(IC_TestRuns) * (1.0d0 / (IcArl + up_gap)))
            idx_U = max(1, min(IC_TestRuns, idx_U))

            call partial_quickselect(values=global_stats, num_all=IC_TestRuns, &
                                     top_k=idx_L, &
                                     order_index=-1, & 
                                     values_sub=global_stats) 


            LimitR = global_stats(idx_U) ! Upper Limit (Smaller Index, Larger Value)
            LimitL = global_stats(idx_L) ! Lower Limit (Larger Index, Smaller Value)

            print '(A,I2,A,F10.4,A,F10.4,A,F10.4)', &
                  " Iter ", iter, ": Current=", Limit_Init, " -> Range [", LimitL, ",", LimitR, "]"
            ! 4. Check Convergence
            if (Limit_Init >= LimitL .and. Limit_Init <= LimitR) then
                range_accepted = .true.
                print *, " >> Range Accepted."
            else
                Limit_Init = 0.5d0 * (LimitL + LimitR)
                print *, " >> Updating Guess -> ", Limit_Init
            end if
        end if   
        call MPI_BCAST(range_accepted, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        if (range_accepted) exit
        call MPI_BCAST(Limit_Init, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    end do

    ! Ensure all procs have final limits
    call MPI_BCAST(LimitL, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(LimitR, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    if (allocated(local_stats)) deallocate(local_stats)
    if (my_rank == 0) then
        deallocate(global_stats)
        deallocate(recv_counts)
        deallocate(displs)
    end if
end subroutine mMSTD_EstimateLimitRange

subroutine Exploration_Evaluation(data_fid, config_fid)
    use, intrinsic :: iso_fortran_env, only: output_unit, real64
    implicit none


    ! --- Input Arguments ---
    integer, intent(in) :: data_fid   
    integer, intent(in) :: config_fid 

    ! --- MPI Variables ---
    integer :: my_rank, num_procs, ierr
    integer :: local_runs, remainder, total_elements
    integer, allocatable :: counts(:), displs(:)

    ! --- Simulation Variables ---
    integer :: i, j, nodes_count, time_seed, seed
    integer :: one_vec(1) 
    integer, allocatable :: sampling_index(:), sampling_mask(:)
    real(dp), allocatable :: OnlineSample(:)
    real(dp) :: alpha_value, charting_statistic, current_maxmin

    ! --- Result Storage ---
    real(dp), allocatable :: Local_DisVec(:,:), Global_DisVec(:,:) 

    ! --- Statistics Processing Variables ---
    real(dp), allocatable :: temp_slice(:)
    type(ARL_Stats_type)  :: current_stats 

    ! --- Progress Display Variables ---
    integer :: t_start, t_now, t_rate, t_ela, t_rem, t_min, t_sec
    integer :: prog_pct, last_pct, bar_len, filled_len
    character(len=20) :: bar_str
    character(len=1)  :: cr

    ! =================================================================
    ! 1. MPI Initialization & Setup
    ! =================================================================
    call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)

    local_runs = IC_TestRuns / num_procs
    remainder  = mod(IC_TestRuns, num_procs)
    if (my_rank < remainder) local_runs = local_runs + 1

    allocate(Local_DisVec(INT(Exploration_IcArl), local_runs))
    allocate(sampling_index(num_samplingnodes))
    allocate(sampling_mask(num_allnodes))
    allocate(OnlineSample(num_allnodes))

    if (my_rank == 0) then
        call system_clock(count_rate=t_rate)
        call system_clock(count=t_start)
        last_pct = -1; bar_len = 20; cr = char(13)
        print *, ">> Simulation Started..."
        write(config_fid, '(A, I4, A)') " MPI Run Started: ", num_procs, " cores."
    end if

    ! =================================================================
    ! 2. Local Simulation Loop
    ! =================================================================
    do i = 1, local_runs
        ! Initialize progress reporting.
        if (my_rank == 0) then
            prog_pct = int((real(i)/real(local_runs))*100.0)
            if (prog_pct > last_pct .or. i == local_runs) then
                call system_clock(count=t_now)
                t_ela = int(real(t_now - t_start) / real(t_rate))
                if (i > 5 .and. t_ela > 0) then
                    t_rem = int(real(t_ela) / real(i) * real(local_runs - i))
                else; t_rem = 0; end if
                t_min = t_rem / 60; t_sec = mod(t_rem, 60)
                filled_len = int(real(prog_pct)/100.0 * real(bar_len))
                if (filled_len < 0) filled_len = 0; if (filled_len > bar_len) filled_len = bar_len
                bar_str = repeat('#', filled_len) // repeat('.', bar_len - filled_len)
                write(output_unit, '(A, A, I3, A, A, A, A, I2.2, A, I2.2, A)', advance='no') &
                    cr, " Status:", prog_pct, "% [", trim(bar_str), "] ", "ETA: ", t_min, "m ", t_sec, "s"
                flush(output_unit)
                last_pct = prog_pct
            end if
        end if

        ! --- Simulation ---
        call get_clock_seed(time_seed)
        seed = mod(time_seed + my_rank*997 + i*1009, 2147483647) + 1
        call RNSET(seed)

        sampling_mask = 0; nodes_count = 0
        do while (nodes_count < num_samplingnodes)
            call RNUND(num_allnodes, one_vec)
            if (sampling_mask(one_vec(1)) == 0) then
                nodes_count = nodes_count + 1
                sampling_index(nodes_count) = one_vec(1)
                sampling_mask(one_vec(1)) = 1
            end if
        end do

        call Reset_mMSTD_State()
        call Set_SampleParams(ICorOC_in=IC, noise_type_in=NOISE_GAUSSIAN, back_index_in=BACK_NONE)

        do j = 1, Time_window
            call GenerateOnlineSample(sampling_index, OnlineSample)
            call mMSTD(alpha_value, OnlineSample, sampling_index, charting_statistic)
        end do

        call Reset_Exploration_State()
        do j = 1, Exploration_IcArl
            call GenerateOnlineSample(sampling_index, OnlineSample)
            call mMSTD(alpha_value, OnlineSample, sampling_index, charting_statistic)
            call ExplorationEvaluation(sampling_index, current_maxmin)
            Local_DisVec(j, i) = current_maxmin
        end do
    end do

    ! =================================================================
    ! 3. Gather & Write
    ! =================================================================
    if (my_rank == 0) then
        write(output_unit, *) ! Newline after progress bar
        print *, ">> Gathering data..."
        allocate(Global_DisVec(INT(Exploration_IcArl), IC_TestRuns))
        allocate(counts(num_procs), displs(num_procs))
    end if

    ! MPI Gather Setup
    total_elements = local_runs * Exploration_IcArl
    call MPI_Gather(total_elements, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    if (my_rank == 0) then
        displs(1) = 0
        do i = 2, num_procs
            displs(i) = displs(i-1) + counts(i-1)
        end do
    end if

    call MPI_Gatherv(Local_DisVec, total_elements, MPI_DOUBLE_PRECISION, &
                     Global_DisVec, counts, displs, MPI_DOUBLE_PRECISION, &
                     0, MPI_COMM_WORLD, ierr)

    if (my_rank == 0) then
        ! Write Data
        do i = 1, IC_TestRuns
            write(data_fid, *) Global_DisVec(:, i)
        end do
        
        ! Log completion
        write(config_fid, '(A)') " Data successfully gathered and written."
        print *, ">> Done."
        
        deallocate(Global_DisVec, counts, displs)
    end if



    deallocate(Local_DisVec, sampling_index, sampling_mask, OnlineSample)


    end subroutine Exploration_Evaluation

subroutine Random_Exploration_Evaluation(data_fid, config_fid)
    use, intrinsic :: iso_fortran_env, only: output_unit, real64
    implicit none


    ! --- Input Arguments ---
    integer, intent(in) :: data_fid   
    integer, intent(in) :: config_fid 

    ! --- MPI Variables ---
    integer :: my_rank, num_procs, ierr
    integer :: local_runs, remainder, total_elements
    integer, allocatable :: counts(:), displs(:)

    ! --- Simulation Variables ---
    integer :: i, j, nodes_count, time_seed, seed
    integer :: one_vec(1) 
    integer, allocatable :: sampling_index(:), sampling_mask(:)
    real(dp) ::  current_maxmin

    ! --- Result Storage ---
    real(dp), allocatable :: Local_DisVec(:,:), Global_DisVec(:,:) 

    ! --- Statistics Processing Variables ---
    real(dp), allocatable :: temp_slice(:)
    type(ARL_Stats_type)  :: current_stats 

    ! --- Progress Display Variables ---
    integer :: t_start, t_now, t_rate, t_ela, t_rem, t_min, t_sec
    integer :: prog_pct, last_pct, bar_len, filled_len
    character(len=20) :: bar_str
    character(len=1)  :: cr

    ! =================================================================
    ! 1. MPI Initialization & Setup
    ! =================================================================
    call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)

    local_runs = IC_TestRuns / num_procs
    remainder  = mod(IC_TestRuns, num_procs)
    if (my_rank < remainder) local_runs = local_runs + 1

    allocate(Local_DisVec(INT(Exploration_IcArl), local_runs))
    allocate(sampling_index(num_samplingnodes))
    allocate(sampling_mask(num_allnodes))

    if (my_rank == 0) then
        call system_clock(count_rate=t_rate)
        call system_clock(count=t_start)
        last_pct = -1; bar_len = 20; cr = char(13)
        print *, ">> Simulation Started..."
        write(config_fid, '(A, I4, A)') " MPI Run Started: ", num_procs, " cores."
    end if

    ! =================================================================
    ! 2. Local Simulation Loop
    ! =================================================================
    do i = 1, local_runs
        ! Initialize progress reporting.
        if (my_rank == 0) then
            prog_pct = int((real(i)/real(local_runs))*100.0)
            if (prog_pct > last_pct .or. i == local_runs) then
                call system_clock(count=t_now)
                t_ela = int(real(t_now - t_start) / real(t_rate))
                if (i > 5 .and. t_ela > 0) then
                    t_rem = int(real(t_ela) / real(i) * real(local_runs - i))
                else; t_rem = 0; end if
                t_min = t_rem / 60; t_sec = mod(t_rem, 60)
                filled_len = int(real(prog_pct)/100.0 * real(bar_len))
                if (filled_len < 0) filled_len = 0; if (filled_len > bar_len) filled_len = bar_len
                bar_str = repeat('#', filled_len) // repeat('.', bar_len - filled_len)
                write(output_unit, '(A, A, I3, A, A, A, A, I2.2, A, I2.2, A)', advance='no') &
                    cr, " Status:", prog_pct, "% [", trim(bar_str), "] ", "ETA: ", t_min, "m ", t_sec, "s"
                flush(output_unit)
                last_pct = prog_pct
            end if
        end if

        ! --- Simulation ---
        call get_clock_seed(time_seed)
        seed = mod(time_seed + my_rank*997 + i*1009, 2147483647) + 1
        call RNSET(seed)

        sampling_mask = 0; nodes_count = 0
        do while (nodes_count < num_samplingnodes)
            call RNUND(num_allnodes, one_vec)
            if (sampling_mask(one_vec(1)) == 0) then
                nodes_count = nodes_count + 1
                sampling_index(nodes_count) = one_vec(1)
                sampling_mask(one_vec(1)) = 1
            end if
        end do


        do j = 1, Time_window
            sampling_mask = 0; nodes_count = 0
            do while (nodes_count < num_samplingnodes)
                call RNUND(num_allnodes, one_vec)
                if (sampling_mask(one_vec(1)) == 0) then
                    nodes_count = nodes_count + 1
                    sampling_index(nodes_count) = one_vec(1)
                    sampling_mask(one_vec(1)) = 1
                end if
            end do
        end do

        call Reset_Exploration_State()
        do j = 1, Exploration_IcArl
            sampling_mask = 0; nodes_count = 0
            do while (nodes_count < num_samplingnodes)
                call RNUND(num_allnodes, one_vec)
                if (sampling_mask(one_vec(1)) == 0) then
                    nodes_count = nodes_count + 1
                    sampling_index(nodes_count) = one_vec(1)
                    sampling_mask(one_vec(1)) = 1
                end if
            end do
            call ExplorationEvaluation(sampling_index, current_maxmin)
            Local_DisVec(j, i) = current_maxmin
        end do
    end do

    ! =================================================================
    ! 3. Gather & Write
    ! =================================================================
    if (my_rank == 0) then
        write(output_unit, *) ! Newline after progress bar
        print *, ">> Gathering data..."
        allocate(Global_DisVec(INT(Exploration_IcArl), IC_TestRuns))
        allocate(counts(num_procs), displs(num_procs))
    end if

    ! MPI Gather Setup
    total_elements = local_runs * Exploration_IcArl
    call MPI_Gather(total_elements, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    if (my_rank == 0) then
        displs(1) = 0
        do i = 2, num_procs
            displs(i) = displs(i-1) + counts(i-1)
        end do
    end if

    call MPI_Gatherv(Local_DisVec, total_elements, MPI_DOUBLE_PRECISION, &
                     Global_DisVec, counts, displs, MPI_DOUBLE_PRECISION, &
                     0, MPI_COMM_WORLD, ierr)

    if (my_rank == 0) then
        ! Write Data
        do i = 1, IC_TestRuns
            write(data_fid, *) Global_DisVec(:, i)
        end do
        
        ! Log completion
        write(config_fid, '(A)') " Data successfully gathered and written."
        print *, ">> Done."
        
        deallocate(Global_DisVec, counts, displs)
    end if



    deallocate(Local_DisVec, sampling_index, sampling_mask)


end subroutine Random_Exploration_Evaluation


end module performance_evaluation
