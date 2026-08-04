!==============================================================
!  Module: performance_evaluation
!  Purpose: Evaluate NAS detection limits and run-length performance.
!==============================================================
module performance_evaluation
    use GlobalSettings_mod
    use SparseMatrix_mod
    use SamplesGeneration_mod
    use NAS_mod
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
!  Purpose: Initialize NAS and SampleGen static data.
!==============================================================
subroutine Init_Performance_Evaluation( &
    NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId, &
    Nodesset, B0_all, B1_all, CovMat)
    
    implicit none
    ! Generic Inputs (Passed to SampleGen)
    real(dp), intent(in) :: NeighborDis(:,:), Kernel_NeighborDis(:,:)
    integer, intent(in) :: NeighborId(:,:), Kernel_NeighborId(:,:)
    
    ! SampleGen Inputs
    real(dp), intent(in), optional :: Nodesset(:,:)
    type(Bbasis_sparse_type), intent(in), optional :: B0_all, B1_all, CovMat

    ! 1. Initialize Sample Generation
    call Init_SamplesGeneration(Nodesset, B0_all, B1_all, CovMat)

    ! 2. Initialize NAS Data
    ! Allocates Smin_one, Smin_two, Smax_one, Smax_two
    call Init_NAS_Data()

    perf_initialized = .true.

end subroutine Init_Performance_Evaluation

!==============================================================
!  Subroutine: Clean_Performance_Evaluation
!==============================================================
subroutine Clean_Performance_Evaluation()
    implicit none
    call Clean_NAS()          ! Release NAS state
    call Clean_SamplesGeneration()
    perf_initialized = .false.
end subroutine Clean_Performance_Evaluation

!==============================================================
!  Subroutine: NAS_EstimateLimitRange
!  Purpose: Estimate Limit bracket [L, R] for Binary Search.
!==============================================================
subroutine NAS_EstimateLimitRange(noise_type, back_index, &
                                      Limit_Init, LimitL, LimitR, &
                                      back_bandwidth)
    implicit none

    !---------------- Inputs ----------------
    integer, intent(in) :: noise_type, back_index
    real(dp), intent(inout) :: Limit_Init 
    real(dp), intent(in), optional :: back_bandwidth

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
    real(dp) :: loc_back_bw
    integer :: target_idx_L, target_idx_U

    if (.not. perf_initialized) stop "Error: Perf Eval not initialized."

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
    is_master = (my_rank == 0)

    loc_back_bw = 1.0d0
    if (present(back_bandwidth)) loc_back_bw = back_bandwidth

    if (is_master) then
        allocate(all_max_stats(IC_TestRuns))
        print *, "==== Auto-Estimating Limit Range (NAS) ===="
    end if
    
    ! 2. Execute Parallel Simulations
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
            
            call get_clock_seed(time_seed)
            seed = mod(time_seed + task*313 + my_rank*997, 2147483647) + 1
            call RNSET(seed)

            ! Init Random Sampling
            sampling_mask = 0; nodes_count = 0
            do while (nodes_count < num_samplingnodes)
                call RNUND(num_allnodes, one_vec)
                if (sampling_mask(one_vec(1)) == 0) then
                    nodes_count = nodes_count + 1
                    sampling_index(nodes_count) = one_vec(1)
                    sampling_mask(one_vec(1)) = 1
                end if
            end do

            ! Reset NAS state before the run.
            call Reset_NAS_State()
            
            ! Warm-up Phase
            call Set_SampleParams(ICorOC_in=IC, noise_type_in=noise_type, &
                                  back_index_in=back_index, &
                                  back_bandwidth_in=loc_back_bw)

            do j = 1, Time_window
                call GenerateOnlineSample(sampling_index, OnlineSample)
                call NAS(OnlineSample, sampling_index, charting_statistic)
            end do                

            sim_max_stat = -1.0d20
            
            ! Measurement Phase
            do j = 1, Time_window
                call GenerateOnlineSample(sampling_index, OnlineSample)
                call NAS(OnlineSample, sampling_index, charting_statistic)
                
                if (charting_statistic > sim_max_stat) sim_max_stat = charting_statistic
            end do
            
            call MPI_SEND(sim_max_stat, 1, MPI_DOUBLE_PRECISION, 0, RESULT_TAG, MPI_COMM_WORLD, ierr)
        end do
    end if
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    ! 3. Process Results (Master Only)
    if (is_master) then
        ! Calculate quantiles for brackets
        target_idx_L = int(dble(IC_TestRuns) * (1.0d0*Time_window/(IcArl-50.0d0)))
        target_idx_L = max(1, min(IC_TestRuns, target_idx_L))
        
        target_idx_U = int(dble(IC_TestRuns) * (1.0d0*Time_window/(IcArl+50.0d0)))
        target_idx_U = max(1, min(IC_TestRuns, target_idx_U))
        
        call partial_quickselect(values=all_max_stats, num_all=IC_TestRuns, &
                                 top_k=IC_TestRuns-target_idx_L+1, &
                                 order_index=-1, values_sub=all_max_stats)
        
        LimitL = all_max_stats(target_idx_L)
        LimitR = all_max_stats(target_idx_U) 

        print '(A,F10.4,A,F10.4,A)', " Estimated Range [", LimitL, ",", LimitR, "]"
        Limit_Init = 0.5d0 * (LimitL + LimitR) ! Update Init guess
        deallocate(all_max_stats)
    end if

    ! Broadcast results
    call MPI_BCAST(LimitL, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(LimitR, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(Limit_Init, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

end subroutine NAS_EstimateLimitRange

!==============================================================
!  Subroutine: NAS_limitSearch
!==============================================================
subroutine NAS_limitSearch(noise_type, back_index, &
                             BLimitL, BLimitR, Limit, arl, std, fid, &
                             back_bandwidth)
    implicit none

    integer, intent(in) :: noise_type, back_index
    integer, intent(in), optional :: fid
    real(dp), intent(in) :: BLimitL, BLimitR
    real(dp), intent(in), optional :: back_bandwidth
    real(dp), intent(out) :: Limit, arl, std

    integer :: ierr, num_procs, my_rank, status(MPI_STATUS_SIZE)
    integer :: task, next_task, active_workers, worker_rank
    integer :: recv_count
    logical :: is_master, is_converged
    integer, parameter :: WORK_TAG = 1, STOP_TAG = 2, RESULT_TAG = 3

    integer :: seed, time_seed, i, j, print_step
    real(dp) :: running_sum, local_result, sim_result, global_arl
    real(dp) :: LimitL, LimitR, charting_statistic
    real(dp), allocatable :: all_results(:)
    
    integer :: sampling_index(num_samplingnodes)
    real(dp) :: OnlineSample(num_samplingnodes)
    integer :: sampling_mask(num_allnodes), nodes_count, one_vec(1)
    real(dp) :: loc_back_bw

    if (.not. perf_initialized) stop "Error: Perf Eval not initialized."

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
    is_master = (my_rank == 0)

    loc_back_bw = 1.0d0; if (present(back_bandwidth)) loc_back_bw = back_bandwidth

    if (is_master) then
        LimitL = BLimitL
        LimitR = BLimitR
        global_arl = 0.0d0
        allocate(all_results(simu))
        print *, "==== Begin NAS_limitSearch ===="
    end if

    call MPI_BCAST(loc_back_bw, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Binary Search Loop
    do         
      if (is_master) then        
            if (abs(global_arl - IcArl) < Icstd) then
                is_converged = .true.
            else
                is_converged = .false.
                Limit = 0.5d0 * (LimitL + LimitR)
                if (present(fid)) then
                    write(fid,*) 'NAS Search: [', LimitL, LimitR, '] -> Testing:', Limit
                else
                    print '(A,F12.6,A,F12.6,A,F12.6)', " NAS Search: [", LimitL, ",", LimitR, "] Try:", Limit
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
                call MPI_RECV(local_result, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, RESULT_TAG, &
                            MPI_COMM_WORLD, status, ierr)
                worker_rank = status(MPI_SOURCE)
                
                recv_count = recv_count + 1
                if (recv_count <= simu) then
                    all_results(recv_count) = local_result
                    running_sum = running_sum + local_result
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
            if (present(fid)) write(fid,*) 'Result ARL=', global_arl

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

                ! Reset NAS state before the run.
                call Reset_NAS_State()
                
                call Set_SampleParams(ICorOC_in=IC, noise_type_in=noise_type, &
                                      back_index_in=back_index, &
                                      back_bandwidth_in=loc_back_bw)

                ! Warm-up
                do j = 1, Time_window
                    call GenerateOnlineSample(sampling_index, OnlineSample)
                    call NAS(OnlineSample, sampling_index, charting_statistic)
                end do

                ! Monitoring ARL Check
                sim_result = 0.0d0
                do while (charting_statistic < Limit)
                    call GenerateOnlineSample(sampling_index, OnlineSample)
                    call NAS(OnlineSample, sampling_index, charting_statistic)
                    
                    sim_result = sim_result + 1.0d0
                    if (sim_result >= Maxarl) exit
                end do
                
                call MPI_SEND(sim_result, 1, MPI_DOUBLE_PRECISION, 0, RESULT_TAG, MPI_COMM_WORLD, ierr)
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

end subroutine NAS_limitSearch

!==============================================================
!  Subroutine: NAS_OCARL
!==============================================================
subroutine NAS_OCARL(noise_type, back_index, anomaly_type, &
                 theta, shiftradius, shiftvalue, shiftradius_zero, &
                 delta_radius, delta_time, amplitude, &
                 anomaly_bandwidth, back_bandwidth, & 
                 Limit, arl, std, max_arl, q1_arl, q2_arl, q3_arl)
    implicit none
    
    integer, intent(in) :: noise_type, back_index, anomaly_type
    real(dp), intent(in) :: Limit
    
    real(dp), intent(in), optional :: theta, shiftradius, shiftvalue
    real(dp), intent(in), optional :: shiftradius_zero, delta_radius, delta_time
    real(dp), intent(in), optional :: amplitude, anomaly_bandwidth, back_bandwidth

    real(dp), intent(out) :: arl, std, max_arl, q1_arl, q2_arl, q3_arl

    integer :: ierr, num_procs, my_rank, status(MPI_STATUS_SIZE)
    integer :: task, next_task, active_workers, worker_rank
    integer :: i, j, seed, time_seed, recv_count, print_step
    logical :: is_master
    integer, parameter :: WORK_TAG = 1, STOP_TAG = 2, RESULT_TAG = 3
    
    real(dp) :: running_sum, local_result, charting_statistic, sim_result
    real(dp), allocatable :: all_results(:)
    
    integer :: sampling_index(num_samplingnodes)
    real(dp) :: OnlineSample(num_samplingnodes)
    integer :: one_vec(1), simu_index(simu)
    integer :: sampling_mask(num_allnodes), nodes_count
    
    type(AnomalyParams_type) :: local_ano_params
    real(dp) :: loc_back_bw, loc_ano_bw
    ! The anomaly generator currently uses the shared kernel_bandwidth.
    ! loc_ano_bw is retained by the evaluation interface but is not applied independently.

    if (.not. perf_initialized) stop "Error: Perf Eval not initialized."

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
    is_master = (my_rank == 0)

    loc_back_bw = 1.0d0; if (present(back_bandwidth)) loc_back_bw = back_bandwidth
    loc_ano_bw  = 1.0d0; if (present(anomaly_bandwidth)) loc_ano_bw = anomaly_bandwidth

    if (is_master) then
        allocate(all_results(simu))
        all_results = 0.0d0
        print *, "==== Begin NAS_OCARL_dynamic ===="
    end if
    
    call MPI_BCAST(loc_back_bw, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(loc_ano_bw,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    if (is_master) then
        next_task = 1
        recv_count = 0
        running_sum = 0.0d0
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
            call MPI_RECV(local_result, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, RESULT_TAG, &
                        MPI_COMM_WORLD, status, ierr)
            worker_rank = status(MPI_SOURCE)
            
            recv_count = recv_count + 1
            if (recv_count <= simu) then
                all_results(recv_count) = local_result
                running_sum = running_sum + local_result
            end if
            
            if (mod(recv_count, print_step) == 0) then
                print '(A,I0,A,I0,A,F10.4)', " OCARL Progress: ", recv_count, "/", simu, " Avg:", running_sum/dble(recv_count)
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

            ! 2. Reset
            call Reset_NAS_State()

            ! 3. Warm-up
            call Set_SampleParams(ICorOC_in=IC, noise_type_in=noise_type, &
                                  back_index_in=back_index, &
                                  back_bandwidth_in=loc_back_bw)

            do j = 1, Time_window
                call GenerateOnlineSample(sampling_index, OnlineSample)
                call NAS(OnlineSample, sampling_index, charting_statistic)
            end do

            ! 4. Anomaly Params
            local_ano_params%type_id = anomaly_type
            
            select case (anomaly_type)
            case (TYPE_CIRCLE, TYPE_ST, TYPE_CIRCLE_CONV, TYPE_SPATIAL)
                call RNUND(num_allnodes, one_vec)
                local_ano_params%center_idx = one_vec(1)
            case (TYPE_BSPLINE)
                call RNUND(max(1,knot1_square_TS), one_vec)
                local_ano_params%bspline_idx = min(knot1_square_TS, one_vec(1))
            end select

            if (present(shiftradius))      local_ano_params%radius = shiftradius
            if (present(shiftvalue))       local_ano_params%value = shiftvalue
            if (present(theta))            local_ano_params%theta = theta
            if (present(shiftradius_zero)) local_ano_params%radius = shiftradius_zero
            if (present(delta_radius))     local_ano_params%delta_r = delta_radius
            if (present(delta_time))       local_ano_params%delta_t = delta_time
            if (present(amplitude))        local_ano_params%value = amplitude
            

            ! 5. Monitoring
            sim_result = 0.0d0
            do while (charting_statistic < Limit)
                local_ano_params%time_idx = int(sim_result)
                call Set_SampleParams(ICorOC_in=OC, AnoParams_in=local_ano_params)
                
                call GenerateOnlineSample(sampling_index, OnlineSample)
                call NAS(OnlineSample, sampling_index, charting_statistic)
                
                sim_result = sim_result + 1.0d0
                if (sim_result >= Maxarl) exit
            end do
            
            call MPI_SEND(sim_result, 1, MPI_DOUBLE_PRECISION, 0, RESULT_TAG, MPI_COMM_WORLD, ierr)
        end do
    end if

    if (is_master) then
        arl = sum(all_results) / dble(simu)
        simu_index = [(i, i=1,simu)]
        call SVRGP(all_results, all_results, simu_index)
        q1_arl  = all_results(max(1,int(0.25d0*simu)))
        q2_arl  = all_results(max(1,int(0.50d0*simu)))
        q3_arl  = all_results(max(1,int(0.75d0*simu)))
        max_arl = all_results(simu)

        std = 0.0d0
        do i=1,simu
            std = std + (all_results(i)-arl)**2.0d0
        end do
        std = sqrt(std / dble(max(1,simu-1))) / sqrt(dble(simu))
        deallocate(all_results)
    end if

    call MPI_BCAST(arl, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(std, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(q1_arl, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(q2_arl, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(q3_arl, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(max_arl, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

end subroutine NAS_OCARL

!==============================================================
!  Subroutine: NAS_ICPerformance
!==============================================================
subroutine NAS_ICPerformance(noise_type, back_index, Limit, falserate, back_bandwidth)
    implicit none

    integer, intent(in) :: noise_type, back_index
    real(dp), intent(in) :: Limit
    real(dp), intent(out) :: falserate
    real(dp), intent(in), optional :: back_bandwidth

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
    real(dp) :: loc_back_bw

    if (.not. perf_initialized) stop "Error: Perf Eval not initialized."

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
    is_master = (my_rank == 0)

    loc_back_bw = 1.0d0; if (present(back_bandwidth)) loc_back_bw = back_bandwidth

    if (is_master) then
        print *, "==== Begin NAS_ICPerformance_dynamic ===="
        global_falsecount = 0
    end if

    call MPI_BCAST(loc_back_bw, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

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
        ! Worker Logic
        do
            call MPI_RECV(run, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
            if (status(MPI_TAG) == STOP_TAG) exit
            
            call get_clock_seed(time_seed)
            seed = mod(time_seed + run*7919 + my_rank*99991, 2147483647) + 1
            call RNSET(seed)

            ! Random Init
            sampling_mask = 0; nodes_count = 0
            do while (nodes_count < num_samplingnodes)
                call RNUND(num_allnodes, one_vec)
                if (sampling_mask(one_vec(1)) == 0) then
                    nodes_count = nodes_count + 1
                    sampling_index(nodes_count) = one_vec(1)
                    sampling_mask(one_vec(1)) = 1
                end if
            end do

            ! Reset NAS state before the run.
            call Reset_NAS_State()
            
            call Set_SampleParams(ICorOC_in=IC, noise_type_in=noise_type, &
                                  back_index_in=back_index, &
                                  back_bandwidth_in=loc_back_bw)

            ! Warmup
            do j = 1, Time_window
                call GenerateOnlineSample(sampling_index, OnlineSample)
                call NAS(OnlineSample, sampling_index, charting_statistic)
            end do
            
            ! Single check for FAR
            call GenerateOnlineSample(sampling_index, OnlineSample)
            call NAS(OnlineSample, sampling_index, charting_statistic)
            
            local_falsecount = 0
            if (charting_statistic > Limit) local_falsecount = 1
            
            false_result(1) = local_falsecount
            false_result(2) = 1
            call MPI_SEND(false_result, 2, MPI_INTEGER, 0, RESULT_TAG, MPI_COMM_WORLD, ierr)
        end do
    end if

    call MPI_BCAST(falserate, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

end subroutine NAS_ICPerformance

end module performance_evaluation
