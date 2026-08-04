module mMSTD_mod
    use GlobalSettings_mod
    use RNNOF_INT
    use RNUND_INT
    use RNBIN_INT
    use SVRGP_INT
    use utils_mod
    use RNNOR_INT
    
    implicit none
    public 
    save    

    !---------------------- mMKSTD parameters ----------------------
    real(dp) :: Epah_mM = 0.02d0, alpha_mM = 4.0d0, pho_mM = 0.0001d0
    real(dp) :: sigmaprio_mM = 0.1d0, lambda_mM = 0.90d0
    integer  :: Topr_mM = 1
    real(dp) :: lambdasq_mM = 0.81
    
    ! These defaults can be set through set_mMSTDparams.
    integer :: near_Threshold_mM = 30,num_kelnear_mM = 15
    real(dp) :: Kernel_threshold_mM = 0.0d0           

    logical  :: initialized_mMSTD = .false.
    logical  :: data_loaded_mMSTD = .false.

    !---------------------- INTERNAL STATIC DATA --------------------
    ! Store geometric neighbors once to avoid passing them every time
    real(dp), allocatable :: mM_NeighborDis(:,:)
    integer,  allocatable :: mM_NeighborId(:,:)
    real(dp), allocatable :: mM_Kernel_NeighborDis(:,:)
    integer,  allocatable :: mM_Kernel_NeighborId(:,:)

    !---------------------- INTERNAL DYNAMIC STATE ------------------
    ! Store probability/distance terms internally to maintain state
    real(dp), allocatable :: mM_ProbabilityTerm(:)
    real(dp), allocatable :: mM_KernelVarianceTerm(:)
    real(dp), allocatable :: mM_DistanceTerm(:)
    real(dp), allocatable :: Maxmin_DistanceTerm(:)
    real(dp), allocatable :: mM_stdTerm(:)

contains    
    
!==================================================================
!  Subroutine: set_mMSTDparams
!  Purpose: Initialize or update mMKSTD-specific scalar parameters.
!==================================================================
subroutine set_mMSTDparams( &
    Epah_mM_in, alpha_mM_in, pho_mM_in, sigmaprio_mM_in, lambda_mM_in, &
    Topr_mM_in, near_Threshold_mM_in, Kernel_threshold_mM_in, &
    num_kelnear_mM_in)

    implicit none
    !------------------ Optional Inputs ------------------
    real(dp), intent(in), optional :: Epah_mM_in, alpha_mM_in
    real(dp), intent(in), optional :: pho_mM_in, sigmaprio_mM_in, lambda_mM_in
    integer, intent(in), optional :: Topr_mM_in, near_Threshold_mM_in
    integer, intent(in), optional :: num_kelnear_mM_in
    real(dp), intent(in), optional :: Kernel_threshold_mM_in

    !------------------ Parameter Assignments ------------------
    if (present(Epah_mM_in)) Epah_mM = Epah_mM_in
    if (present(alpha_mM_in)) alpha_mM = alpha_mM_in
    if (present(pho_mM_in)) pho_mM = pho_mM_in
    if (present(sigmaprio_mM_in)) sigmaprio_mM = sigmaprio_mM_in
    if (present(lambda_mM_in)) lambda_mM = lambda_mM_in
    if (present(Topr_mM_in)) Topr_mM = Topr_mM_in
    if (present(near_Threshold_mM_in)) near_Threshold_mM = near_Threshold_mM_in
    if (present(Kernel_threshold_mM_in)) Kernel_threshold_mM = Kernel_threshold_mM_in
    if (present(num_kelnear_mM_in)) num_kelnear_mM = num_kelnear_mM_in

    !------------------ Derived values ------------------
    lambdasq_mM = lambda_mM ** 2.0d0
    initialized_mMSTD = .true.

end subroutine set_mMSTDparams

!==================================================================
!  Subroutine: Init_mMSTD_Data
!  Purpose: 
!     1. Allocate and copy geometric matrices (NeighborDis, etc.) ONE TIME.
!     2. Allocate internal state vectors (ProbabilityTerm, etc.).
!     3. Reset state variables to initial values.
!  Note: This must be called before the first call to mMSTD.
!==================================================================
subroutine Init_mMSTD_Data(NeighborDis_in, NeighborId_in, &
                           Kernel_NeighborDis_in, Kernel_NeighborId_in)
    implicit none
    real(dp), intent(in) :: NeighborDis_in(:,:)
    integer, intent(in)  :: NeighborId_in(:,:)
    real(dp), intent(in) :: Kernel_NeighborDis_in(:,:)
    integer, intent(in)  :: Kernel_NeighborId_in(:,:)

    ! --- 1. Allocate and Copy Static Geometry Data ---
    if (allocated(mM_NeighborDis)) deallocate(mM_NeighborDis)
    if (allocated(mM_NeighborId))  deallocate(mM_NeighborId)
    if (allocated(mM_Kernel_NeighborDis)) deallocate(mM_Kernel_NeighborDis)
    if (allocated(mM_Kernel_NeighborId))  deallocate(mM_Kernel_NeighborId)

    ! Allocate based on global node count and input dimensions
    allocate(mM_NeighborDis(num_allnodes, num_disnearspatial))
    allocate(mM_NeighborId(num_allnodes, num_disnearspatial))
    allocate(mM_Kernel_NeighborDis(num_allnodes, num_kelnear_mM))
    allocate(mM_Kernel_NeighborId(num_allnodes, num_kelnear_mM))

    mM_NeighborDis = NeighborDis_in
    mM_NeighborId  = NeighborId_in
    mM_Kernel_NeighborDis = Kernel_NeighborDis_in
    mM_Kernel_NeighborId  = Kernel_NeighborId_in

    ! --- 2. Allocate Dynamic State Variables ---
    if (allocated(mM_ProbabilityTerm)) deallocate(mM_ProbabilityTerm)
    if (allocated(mM_KernelVarianceTerm)) deallocate(mM_KernelVarianceTerm)
    if (allocated(mM_DistanceTerm)) deallocate(mM_DistanceTerm)
    if (allocated(mM_stdTerm)) deallocate(mM_stdTerm)

    allocate(mM_ProbabilityTerm(num_allnodes))
    allocate(mM_KernelVarianceTerm(num_allnodes))
    allocate(mM_DistanceTerm(num_allnodes))
    allocate(mM_stdTerm(num_allnodes))

    ! --- 3. Initialize State ---
    call Reset_mMSTD_State()

    data_loaded_mMSTD = .true.

end subroutine Init_mMSTD_Data

!==================================================================
!  Subroutine: Reset_mMSTD_State
!  Purpose: Resets the dynamic probability and distance terms to 
!           initial conditions without reloading the geometry.
!==================================================================
subroutine Reset_mMSTD_State()
    implicit none
    if (.not. allocated(mM_ProbabilityTerm)) return

    mM_ProbabilityTerm    = 0.0d0
    mM_KernelVarianceTerm = 0.0d0
    mM_stdTerm            = 0.0d0
    
    ! Initialize DistanceTerm to a very large number so that the first
    ! minimum distance check (temp_val <= DistanceTerm) succeeds.
    mM_DistanceTerm       = huge(1.0_dp) 
end subroutine Reset_mMSTD_State

!==================================================================
!  Subroutine: Clean_mMSTD
!  Purpose: Frees allocated memory.
!==================================================================
subroutine Clean_mMSTD()
    implicit none
    if (allocated(mM_NeighborDis)) deallocate(mM_NeighborDis)
    if (allocated(mM_NeighborId))  deallocate(mM_NeighborId)
    if (allocated(mM_Kernel_NeighborDis)) deallocate(mM_Kernel_NeighborDis)
    if (allocated(mM_Kernel_NeighborId))  deallocate(mM_Kernel_NeighborId)

    if (allocated(mM_ProbabilityTerm)) deallocate(mM_ProbabilityTerm)
    if (allocated(mM_KernelVarianceTerm)) deallocate(mM_KernelVarianceTerm)
    if (allocated(mM_DistanceTerm)) deallocate(mM_DistanceTerm)
    if (allocated(mM_stdTerm)) deallocate(mM_stdTerm)

    data_loaded_mMSTD = .false.
end subroutine Clean_mMSTD

!==================================================================
!  Subroutine: get_mMSTDparams
!  Purpose: Print all mMKSTD parameters and initialization status.
!==================================================================
subroutine get_mMSTDparams(fid)
    implicit none
    integer, intent(in) :: fid  ! File handle for output

    write(fid,'(A)') "==============================================="
    write(fid,'(A)') "          mMKSTD PARAMETER CONFIGURATION"
    write(fid,'(A)') "==============================================="
    write(fid,'(A,I8)') "num_samplingnodes   = ", num_samplingnodes
    write(fid,'(A,F12.6)') "Epah_mM              = ", Epah_mM
    write(fid,'(A,F12.6)') "alpha_mM             = ", alpha_mM
    write(fid,'(A,F12.6)') "pho_mM               = ", pho_mM
    write(fid,'(A,F12.6)') "sigmaprio_mM         = ", sigmaprio_mM
    write(fid,'(A,F12.6)') "lambda_mM            = ", lambda_mM
    write(fid,'(A,F12.6)') "lambdasq_mM          = ", lambdasq_mM
    write(fid,'(A,I8)')   "Topr_mM              = ", Topr_mM
    write(fid,'(A,I8)')   "near_Threshold_mM    = ", near_Threshold_mM
    write(fid,'(A,I8)')   "num_kelnear_mM       = ", num_kelnear_mM
    write(fid,'(A,F12.6)') "Kernel_threshold_mM  = ", Kernel_threshold_mM
    write(fid,'(A,L1)')   "initialized_mMKSTD   = ", initialized_mMSTD
    write(fid,'(A,L1)')   "data_loaded_mMKSTD   = ", data_loaded_mMSTD

    write(fid,'(A)') "==============================================="
    write(fid,'(A)') "End of mMKSTD PARAMETER listing"
    write(fid,'(A)') "==============================================="

end subroutine get_mMSTDparams

!==================================================================
!  Subroutine: mMSTD
!  Purpose: 
!     Performs Adaptive Sampling using internal state variables.
!     The interface accepts the current sample and updates module state.
!
!  Arguments:
!     alpha_value        : (In) Scaling factor
!     OnlineSample       : (In) The sample values observed at current nodes
!     sampling_index     : (InOut) 
!                          - Input: Indices where OnlineSample was taken
!                          - Output: Indices selected for NEXT step
!     charting_statistic : (Out) Resulting statistic
!==================================================================
subroutine mMSTD(alpha_value, OnlineSample, &
                 sampling_index, charting_statistic)

    implicit none
    !===================== Input parameters ======================
    real(dp), intent(in) :: alpha_value
    real(dp), intent(in) :: OnlineSample(num_samplingnodes)

    !===================== In/Out parameters =====================
    ! sampling_index is used to locate current samples for Step 1,
    ! and then updated to store new samples in Step 3.
    integer,  intent(inout) :: sampling_index(num_samplingnodes)
    real(dp), intent(out)   :: charting_statistic

    !===================== Local variables ========================
    integer :: i, j, pick_idx, curr_idx, neighbor_idx
    integer :: near_index

    real(dp) :: TemVector(num_allnodes)
    real(dp) :: MeanTerm(num_allnodes)
    real(dp) :: noise(num_allnodes)
    
    real(dp) :: statisticmax, statisticmin, temp_val

    ! Safety Check
    if (.not. data_loaded_mMSTD) then
        print *, "Error: mMKSTD data not initialized. Call Init_mMSTD_Data first."
        stop
    end if

    !======================================================================
    ! STEP 1. Kernel-weighted accumulation update
    !         Updates spatial probability and variance based on current samples.
    !         Uses internal arrays: mM_Kernel_NeighborId, mM_ProbabilityTerm, etc.
    !======================================================================
    do i = 1, num_samplingnodes
        curr_idx = sampling_index(i) 

        do j = 1, num_kelnear_mM
            neighbor_idx = mM_Kernel_NeighborId(curr_idx, j)
            
            ! Vectorized-like accumulation for scalar elements
            mM_ProbabilityTerm(neighbor_idx) = mM_ProbabilityTerm(neighbor_idx) + &
                                            mM_Kernel_NeighborDis(curr_idx, j) * OnlineSample(i)

            mM_KernelVarianceTerm(neighbor_idx) = mM_KernelVarianceTerm(neighbor_idx) + &
                                               mM_Kernel_NeighborDis(curr_idx, j)**2

            ! Early stopping based on kernel threshold
            if (mM_Kernel_NeighborDis(curr_idx, j) <= Kernel_threshold_mM) exit
        end do
    end do

    !======================================================================
    ! STEP 2. Compute two-sided charting statistic 
    !======================================================================
    ! Vectorized calculation of standardized term
    mM_stdTerm = mM_ProbabilityTerm / sqrt(mM_KernelVarianceTerm + sigmaprio_mM)

    ! Compute Upper Statistic
    TemVector = mM_stdTerm 
    call partial_quickselect(values=TemVector, num_all=num_allnodes, &
                             top_k=Topr_mM, order_index=-1, values_sub=TemVector(1:Topr_mM))
    statisticmax = sum(TemVector(1:Topr_mM))
    
    ! Compute Lower Statistic
    TemVector = mM_stdTerm 
    call partial_quickselect(values=TemVector, num_all=num_allnodes, &
                             top_k=Topr_mM, order_index=1, values_sub=TemVector(1:Topr_mM))
    statisticmin = sum(TemVector(1:Topr_mM))
    
    charting_statistic = max(statisticmax, abs(statisticmin))

    !======================================================================
    ! STEP 3. Update adaptive sampling metric and select nodes
    !         Updates exploration (Distance) and exploitation (Mean) terms
    !======================================================================

    ! 3.1 Global Vector Updates
    mM_ProbabilityTerm    = mM_ProbabilityTerm * lambda_mM
    mM_KernelVarianceTerm = mM_KernelVarianceTerm * lambdasq_mM

    ! Uniform growth of exploration term (Distance increases over time)
    mM_DistanceTerm = mM_DistanceTerm + pho_mM

    ! Re-calculate detection score for next sampling decision
    mM_stdTerm  = (mM_ProbabilityTerm**2) / (mM_KernelVarianceTerm + sigmaprio_mM)
    MeanTerm = exp(alpha_value * mM_stdTerm)

    ! Initial Sampling Metric
    TemVector = mM_DistanceTerm * MeanTerm

    ! 3.2 Add Jitter 
    call RNNOR(noise)
    TemVector = TemVector + eps_perturb * noise

    ! 3.3 Greedy Node Selection Loop
    do i = 1, num_samplingnodes

        ! Pick the node with the highest metric
        pick_idx = maxloc(TemVector, dim=1)
        sampling_index(i) = pick_idx

        ! Update local distance information for the chosen node's neighbors
        near_index = 0
        
        do j = 1, num_disnearspatial
            neighbor_idx = mM_NeighborId(pick_idx, j)
            temp_val     = mM_NeighborDis(pick_idx, j)

            ! If the neighbor is closer to the new sample than to any previous sample
            ! (Uses internal mM_DistanceTerm)
            if (temp_val <= mM_DistanceTerm(neighbor_idx)) then
                mM_DistanceTerm(neighbor_idx) = temp_val
                
                ! Update the metric for this neighbor immediately
                TemVector(neighbor_idx) = temp_val * MeanTerm(neighbor_idx)
                
                near_index = 0 ! Reset counter since we found a relevant neighbor
            else
                near_index = near_index + 1
            end if
            
            ! Exit after the required neighborhood has been processed.
            if (near_index > near_Threshold_mM) exit
        end do

        ! Mask the selected node so it is not picked again
        TemVector(pick_idx) = -huge(1.0_dp) 
    end do
end subroutine mMSTD
!======================================================================
!  Subroutine: mMSTDII (Hybrid: Greedy Warm Start + ILP Optimization)
!  Purpose: 
!     Performs Adaptive Sampling.
!     1. Runs Greedy Strategy first to establish a baseline coverage (z_max).
!     2. Uses Binary Search + ILP (Gurobi) to find if a strictly better 
!        coverage radius (z < z_max) is feasible.
!     3. If ILP finds a partial solution, fills remainder using Greedy logic.
!======================================================================
subroutine mMSTDII(alpha_value, OnlineSample, &
                   sampling_index, charting_statistic)

    use Gurobi_Interface_mod 
    use iso_c_binding
    implicit none

    !===================== Input parameters ======================
    real(dp), intent(in) :: alpha_value
    real(dp), intent(in) :: OnlineSample(num_samplingnodes)

    !===================== In/Out parameters =====================
    integer,  intent(inout) :: sampling_index(num_samplingnodes)
    real(dp), intent(out)   :: charting_statistic

    !===================== Local variables ========================
    integer :: i, j, curr_idx, neighbor_idx
    
    ! Stats variables
    real(dp) :: TemVector(num_allnodes)
    real(dp) :: statisticmax, statisticmin
    
    ! Optimization variables
    real(dp) :: WeightTerm(num_allnodes)   
    real(dp) :: MetricState(num_allnodes) 
    
    ! Binary Search & ILP variables
    real(dp) :: z_min, z_max, z_mid, req_dist, temp_val
    integer  :: iter, max_bin_iter
    integer  :: current_sol_count
    
    ! Gurobi related (MUST be c_int for C interoperability)
    integer(c_int), allocatable :: row_ptr(:), col_idx(:), sol_vec(:)
    integer(c_int) :: n_cand_c, n_targ_c, stat_c, nnz 
    
    ! CSR Construction Helpers
    integer :: u_cnt, k, cand_node
    integer, allocatable :: target_nodes(:) 
    integer :: max_nnz_est ! For safe allocation
    
    ! Solution storage
    integer :: best_sampling_index(num_samplingnodes)
    integer :: greedy_index(num_samplingnodes)
    logical :: solution_found

    ! Greedy Temp Variables (for Phase A)
    real(dp) :: MeanTerm(num_allnodes) ! Used as temp weight storage
    real(dp) :: Greedy_Dist(num_allnodes)
    real(dp) :: noise(num_allnodes)
    integer :: pick_idx, near_index

    ! Safety Check
    if (.not. data_loaded_mMSTD) then
        print *, "Error: mMKSTD data not initialized. Call Init_mMSTD_Data first."
        stop
    end if

    !======================================================================
    ! STEP 1 & 2: Update State & Calculate Statistic (Standard Logic)
    !======================================================================
    do i = 1, num_samplingnodes
        curr_idx = sampling_index(i) 
        do j = 1, num_kelnear_mM
            neighbor_idx = mM_Kernel_NeighborId(curr_idx, j)
            mM_ProbabilityTerm(neighbor_idx) = mM_ProbabilityTerm(neighbor_idx) + &
                                            mM_Kernel_NeighborDis(curr_idx, j) * OnlineSample(i)
            mM_KernelVarianceTerm(neighbor_idx) = mM_KernelVarianceTerm(neighbor_idx) + &
                                               mM_Kernel_NeighborDis(curr_idx, j)**2
            if (mM_Kernel_NeighborDis(curr_idx, j) <= Kernel_threshold_mM) exit
        end do
    end do

    mM_stdTerm = mM_ProbabilityTerm / sqrt(mM_KernelVarianceTerm + sigmaprio_mM)
    
    ! Calc Charting Statistic
    TemVector = mM_stdTerm 
    call partial_quickselect(values=TemVector, num_all=num_allnodes, &
                             top_k=Topr_mM, order_index=-1, values_sub=TemVector(1:Topr_mM))
    statisticmax = sum(TemVector(1:Topr_mM))
    
    TemVector = mM_stdTerm 
    call partial_quickselect(values=TemVector, num_all=num_allnodes, &
                             top_k=Topr_mM, order_index=1, values_sub=TemVector(1:Topr_mM))
    statisticmin = sum(TemVector(1:Topr_mM))
    
    charting_statistic = max(statisticmax, abs(statisticmin))

    !======================================================================
    ! STEP 3. Optimal Sampling (HYBRID STRATEGY)
    !======================================================================

    ! 3.1 Global State Aging
    mM_ProbabilityTerm    = mM_ProbabilityTerm * lambda_mM
    mM_KernelVarianceTerm = mM_KernelVarianceTerm * lambdasq_mM
    mM_DistanceTerm       = mM_DistanceTerm + pho_mM 

    ! 3.2 Prepare Weights and Base Metric
    mM_stdTerm  = (mM_ProbabilityTerm**2) / (mM_KernelVarianceTerm + sigmaprio_mM)
    WeightTerm  = exp(alpha_value * mM_stdTerm)
    MetricState = WeightTerm * mM_DistanceTerm 

    !----------------------------------------------------------------------
    ! PHASE A: Run Greedy Algorithm First (Warm Start)
    !----------------------------------------------------------------------
    ! We simulate the greedy selection to find a baseline z_max.
    ! This ensures we have a fallback solution and tightens the search range.
    
    Greedy_Dist = mM_DistanceTerm 
    MeanTerm = WeightTerm ! Reuse array to store weights for greedy calc
    
    ! Add a small perturbation to break symmetric scores.
    TemVector = MetricState
    call RNNOR(noise)
    TemVector = TemVector + eps_perturb * noise
    
    do i = 1, num_samplingnodes
        pick_idx = maxloc(TemVector, dim=1)
        greedy_index(i) = pick_idx
        
        ! Greedy Update of local distances
        near_index = 0
        do j = 1, num_disnearspatial
            neighbor_idx = mM_NeighborId(pick_idx, j)
            temp_val     = mM_NeighborDis(pick_idx, j)
            
            if (temp_val <= Greedy_Dist(neighbor_idx)) then
                Greedy_Dist(neighbor_idx) = temp_val
                ! Update metric for next pick
                TemVector(neighbor_idx) = temp_val * MeanTerm(neighbor_idx) 
                near_index = 0
            else
                near_index = near_index + 1
            end if
            if (near_index > near_Threshold_mM) exit
        end do
        TemVector(pick_idx) = -huge(1.0_dp)
    end do
    
    ! Calculate the Max Weighted Distance resulting from Greedy
    ! z_greedy = max( Weight(r) * MinDist(r, GreedySet) )
    z_max = 0.0d0
    do i = 1, num_allnodes
        temp_val = WeightTerm(i) * Greedy_Dist(i)
        if (temp_val > z_max) z_max = temp_val
    end do
    
    ! Initialize "best" with Greedy result. 
    ! If ILP fails or crashes, we return this.
    best_sampling_index = greedy_index
    
    !----------------------------------------------------------------------
    ! PHASE B: Binary Search + ILP (Try to improve z_max)
    !----------------------------------------------------------------------
    
    z_min = 0.0d0
    ! z_max is already set to greedy result
    
    max_bin_iter = 10 
    
    ! Allocation (Moved outside loop)
    allocate(target_nodes(num_allnodes))
    allocate(sol_vec(num_allnodes))
    allocate(row_ptr(num_allnodes + 1)) 
    
    ! Allocate space for self-coverage and all neighboring nodes.
    ! Using int64 logic for safety in allocation size calculation
    max_nnz_est = num_allnodes * (num_disnearspatial + 1)
    allocate(col_idx(max_nnz_est))

    do iter = 1, max_bin_iter
        z_mid = 0.5d0 * (z_min + z_max)
        
        ! 1. Identify Uncovered Set U(z_mid)
        u_cnt = 0
        do i = 1, num_allnodes
            if (MetricState(i) > z_mid) then
                u_cnt = u_cnt + 1
                target_nodes(u_cnt) = i
            end if
        end do
        
        ! If z is so high that no nodes need covering, try smaller z
        if (u_cnt == 0) then
            z_max = z_mid
            cycle
        end if

        ! 2. Construct CSR Matrix
        nnz = 0
        row_ptr(1) = 1
        
        do k = 1, u_cnt
            i = target_nodes(k)
            req_dist = z_mid / WeightTerm(i)
            
            ! --- Self cover ---
            nnz = nnz + 1
            ! Safety Check
            if (nnz > max_nnz_est) then
                 print *, "Error: COL_IDX overflow. Increase size."
                 stop
            end if
            col_idx(nnz) = i 
            
            ! --- Neighbor cover ---
            do j = 1, num_disnearspatial
                if (mM_NeighborDis(i, j) > req_dist) exit 
                
                cand_node = mM_NeighborId(i, j)
                nnz = nnz + 1
                ! Safety Check
                if (nnz > max_nnz_est) then
                    print *, "Error: COL_IDX overflow (neighbors)."
                    stop
                end if
                col_idx(nnz) = cand_node
            end do
            
            row_ptr(k+1) = nnz + 1
        end do
        
        ! 3. Call Gurobi Wrapper
        n_cand_c = int(num_allnodes, c_int)
        n_targ_c = int(u_cnt, c_int)
        
        call solve_set_cover_wrapper(n_cand_c, n_targ_c, row_ptr, col_idx, sol_vec, stat_c)
        
        ! 4. Check Result
        if (stat_c == 0) then
            ! Gurobi solved it optimally
            current_sol_count = sum(sol_vec)
            
            if (current_sol_count <= num_samplingnodes) then
                ! FEASIBLE: We found a set of sensors <= budget that satisfies z_mid
                ! This solution is valid, store it.
                
                block
                    integer :: idx, count
                    count = 0
                    do idx = 1, num_allnodes
                        if (sol_vec(idx) == 1) then
                            count = count + 1
                            if (count <= num_samplingnodes) then
                                best_sampling_index(count) = idx
                            end if
                        end if
                    end do
                    
                    ! Smart Fill: If budget remains, use Greedy on residual
                    if (count < num_samplingnodes) then
                       call FillGreedySamples(best_sampling_index, count, num_samplingnodes, &
                                              WeightTerm, mM_DistanceTerm)
                    end if
                end block
                
                ! Current z_mid is feasible, try to find an even smaller (better) z
                z_max = z_mid
            else
                ! Budget exceeded, need larger z
                z_min = z_mid
            end if
        else
            ! Gurobi failed (Time limit / Infeasible model / License error)
            ! Assume z_mid is too strict
            z_min = z_mid
        end if
    end do

    ! Cleanup
    deallocate(row_ptr, col_idx, sol_vec, target_nodes)

    ! 3.4 Update Output and Internal State
    sampling_index = best_sampling_index

    ! 3.5 Update the distance term for the next time step
    !     Updates mM_DistanceTerm based on the FINAL decision.
    do i = 1, num_samplingnodes
        curr_idx = sampling_index(i)
        mM_DistanceTerm(curr_idx) = 0.0d0
        
        do j = 1, num_disnearspatial
            neighbor_idx = mM_NeighborId(curr_idx, j)
            temp_val     = mM_NeighborDis(curr_idx, j)
            
            if (temp_val < mM_DistanceTerm(neighbor_idx)) then
                mM_DistanceTerm(neighbor_idx) = temp_val
            end if
            if (near_Threshold_mM > 0 .and. j >= near_Threshold_mM) exit
        end do
    end do

contains

    !--------------------------------------------------------------
    ! Improved Filler: Uses Greedy strategy to fill remaining slots
    ! This is better than random filling as it reduces secondary risks.
    !--------------------------------------------------------------
    subroutine FillGreedySamples(arr, current_cnt, max_cnt, w_term, dist_term_in)
        integer, intent(inout) :: arr(:)
        integer, intent(in) :: current_cnt, max_cnt
        real(dp), intent(in) :: w_term(:), dist_term_in(:)
        
        real(dp) :: local_dist(num_allnodes)
        real(dp) :: current_metric(num_allnodes)
        integer :: k, f_i, f_j, pick, n_idx
        logical, allocatable :: used(:)
        
        if (current_cnt >= max_cnt) return
        
        allocate(used(num_allnodes))
        used = .false.
        
        ! 1. Initialize local distance map based on already selected ILP nodes
        local_dist = dist_term_in
        do k = 1, current_cnt
            used(arr(k)) = .true.
            pick = arr(k)
            local_dist(pick) = 0.0d0
            do f_j = 1, num_disnearspatial
                n_idx = mM_NeighborId(pick, f_j)
                if (mM_NeighborDis(pick, f_j) < local_dist(n_idx)) then
                    local_dist(n_idx) = mM_NeighborDis(pick, f_j)
                end if
            end do
        end do
        
        ! 2. Greedy Loop for remainder
        k = current_cnt
        do while (k < max_cnt)
            ! Metric = Weight * Current_Distance
            current_metric = w_term * local_dist
            
            ! Mask already used
            do f_i = 1, num_allnodes
                if (used(f_i)) current_metric(f_i) = -1.0d0
            end do
            
            ! Pick max weighted distance
            pick = maxloc(current_metric, dim=1)
            
            k = k + 1
            arr(k) = pick
            used(pick) = .true.
            
            ! Update distances
            local_dist(pick) = 0.0d0
            do f_j = 1, num_disnearspatial
                n_idx = mM_NeighborId(pick, f_j)
                if (mM_NeighborDis(pick, f_j) < local_dist(n_idx)) then
                    local_dist(n_idx) = mM_NeighborDis(pick, f_j)
                end if
            end do
        end do
        
        deallocate(used)
    end subroutine FillGreedySamples

end subroutine mMSTDII                 
end module mMSTD_mod
