module CDS_mod
    use GlobalSettings_mod
    use RNNOF_INT
    use RNUND_INT
    use RNBIN_INT
    use SVRGP_INT
    use linear_operators
    use utils_mod
    use SparseMatrix_mod
    use ANORIN_INT
    use RNNOR_INT
    implicit none
    public 
    save    

    !---------------------- CDS Parameters --------------------- 
    real(dp) :: miu_min_CDS=1.0d0
    real(dp) :: cuparame_CDS=0.5d0 
    real(dp) :: alpha_CDS=0.95d0, Phi_CDS=0.0627d0
    integer  :: Topr_CDS=1

    logical  :: initialized_CDS = .false.
    logical  :: data_loaded_CDS = .false.

    !---------------------- INTERNAL STATIC DATA --------------------
    ! Covariance Structure (Read-only after init)
    real(dp), allocatable :: CDS_Cov_diag(:)             ! (num_allnodes)
    real(dp), allocatable :: CDS_Cov_Neighborval(:,:)    ! (num_allnodes, num_varnear_CDS)
    integer,  allocatable :: CDS_Cov_NeighborId(:,:)      ! (num_allnodes, num_varnear_CDS)
    real(dp), allocatable :: Cov_dense_global (:, :)

    !---------------------- INTERNAL DYNAMIC STATE ------------------
    ! CUSUM Statistics (Updated iteratively)
    real(dp), allocatable :: CDS_cusumdouble_statistic(:,:) ! (num_allnodes, 2)
    real(dp), allocatable :: CDS_cusum_statistic(:)         ! (num_allnodes)

    contains

!==================================================================
!  Subroutine: set_CDSparams
!  Purpose: Initialize or update CDS-specific scalar parameters
!==================================================================
subroutine set_CDSparams(miu_min_CDS_in, Topr_CDS_in, alpha_CDS_in, num_varnear_CDS_in)
    implicit none
    real(dp), intent(in), optional :: miu_min_CDS_in, alpha_CDS_in
    integer, intent(in), optional :: Topr_CDS_in, num_varnear_CDS_in

    if (present(miu_min_CDS_in)) miu_min_CDS = miu_min_CDS_in
    if (present(Topr_CDS_in))    Topr_CDS    = Topr_CDS_in
    if (present(alpha_CDS_in))   alpha_CDS   = alpha_CDS_in
    if (present(num_varnear_CDS_in)) num_varnear_CDS = num_varnear_CDS_in
    
    cuparame_CDS = (miu_min_CDS**2.0d0) / 2.0d0
    Phi_CDS = ANORIN(1.0d0 - alpha_CDS/2.0d0)
    
    initialized_CDS = .true.
end subroutine set_CDSparams

!==================================================================
!  Subroutine: Init_CDS_Data
!  Purpose: 
!     1. Allocate and copy Covariance Data ONE TIME.
!     2. Allocate internal CUSUM state arrays.
!     3. Reset state to initial values.
!==================================================================
subroutine Init_CDS_Data(Cov_diag_in, Cov_Neighborval_in, Cov_NeighborId_in, Cov_dense_global_in)
    implicit none
    real(dp), intent(in) :: Cov_diag_in(:)
    real(dp), intent(in) :: Cov_Neighborval_in(:,:)
    integer, intent(in)  :: Cov_NeighborId_in(:,:)
    real(dp), intent(in) :: Cov_dense_global_in(:,:)
    
    ! Safety Check
    if (num_varnear_CDS == 0) then
        num_varnear_CDS = size(Cov_Neighborval_in, 2)
        print *, "Warning: num_varnear_CDS auto-set to", num_varnear_CDS
    end if

    ! 1. Allocate & Copy Static Data
    if (allocated(CDS_Cov_diag)) deallocate(CDS_Cov_diag)
    allocate(CDS_Cov_diag(num_allnodes))
    CDS_Cov_diag = Cov_diag_in 

    if (allocated(CDS_Cov_Neighborval)) deallocate(CDS_Cov_Neighborval)
    allocate(CDS_Cov_Neighborval(num_allnodes, num_varnear_CDS))
    CDS_Cov_Neighborval = Cov_Neighborval_in 

    if (allocated(CDS_Cov_NeighborId)) deallocate(CDS_Cov_NeighborId)
    allocate(CDS_Cov_NeighborId(num_allnodes, num_varnear_CDS))
    CDS_Cov_NeighborId = Cov_NeighborId_in
    
    if (allocated(Cov_dense_global)) deallocate(Cov_dense_global)
    allocate(Cov_dense_global(num_allnodes, num_allnodes))
    Cov_dense_global = Cov_dense_global_in

    ! 2. Allocate Dynamic State
    if (allocated(CDS_cusumdouble_statistic)) deallocate(CDS_cusumdouble_statistic)
    allocate(CDS_cusumdouble_statistic(num_allnodes, 2))

    if (allocated(CDS_cusum_statistic)) deallocate(CDS_cusum_statistic)
    allocate(CDS_cusum_statistic(num_allnodes))

    ! 3. Initialize State
    call Reset_CDS_State()

    data_loaded_CDS = .true.

end subroutine Init_CDS_Data

!==================================================================
!  Subroutine: Reset_CDS_State
!  Purpose: Reset CUSUM statistics to zero
!==================================================================
subroutine Reset_CDS_State()
    implicit none
    if (allocated(CDS_cusumdouble_statistic)) CDS_cusumdouble_statistic = 0.0d0
    if (allocated(CDS_cusum_statistic))       CDS_cusum_statistic       = 0.0d0
end subroutine Reset_CDS_State

!==================================================================
!  Subroutine: Clean_CDS
!  Purpose: Free all internal memory
!==================================================================
subroutine Clean_CDS()
    implicit none
    if (allocated(CDS_Cov_diag))          deallocate(CDS_Cov_diag)
    if (allocated(CDS_Cov_Neighborval))   deallocate(CDS_Cov_Neighborval)
    if (allocated(CDS_Cov_NeighborId))     deallocate(CDS_Cov_NeighborId)
    if (allocated(CDS_cusumdouble_statistic)) deallocate(CDS_cusumdouble_statistic)
    if (allocated(CDS_cusum_statistic))       deallocate(CDS_cusum_statistic)
    if (allocated(Cov_dense_global))      deallocate(Cov_dense_global)
    data_loaded_CDS = .false.
end subroutine Clean_CDS

!==================================================================
!  Subroutine: get_CDSparams
!  Purpose: Print status
!==================================================================
subroutine get_CDSparams(fid)
    implicit none
    integer, intent(in) :: fid

    write(fid,'(A)') "==============================================="
    write(fid,'(A)') "             CDS PARAMETER CONFIGURATION"
    write(fid,'(A)') "==============================================="
    write(fid,'(A,I8)') "num_samplingnodes   = ", num_samplingnodes
    write(fid,'(A,F12.6)') "miu_min_CDS         = ", miu_min_CDS
    write(fid,'(A,I8)') "Topr_CDS            = ", Topr_CDS
    write(fid,'(A,F12.6)') "alpha_CDS           = ", alpha_CDS
    write(fid,'(A,F12.6)') "Phi_CDS           = ", Phi_CDS
    write(fid,'(A,L1)')   "data_loaded_CDS     = ", data_loaded_CDS
    write(fid,'(A)') "==============================================="
end subroutine get_CDSparams

!==================================================================
!  Subroutine: CDS
!  Purpose: Update chart statistics and select the next sampled nodes.
!  Sampling ties are resolved with a random selection.
!==================================================================
subroutine CDS(OnlineSample, sampling_index, charting_statistic)
    implicit none

    ! --- Inputs ---
    real(dp), intent(in) :: OnlineSample(num_samplingnodes)
    
    ! --- InOuts / Outputs ---
    integer , intent(inout) :: sampling_index(num_samplingnodes)
    real(dp), intent(out)   :: charting_statistic

    ! --- Local variables ---
    integer :: i, j, k, max_id, idx, ii, jj
    integer :: sampling_iperm(num_allnodes)
    
    ! Buffers for conditional inference
    integer  :: sub_Covnum(num_allnodes)
    integer  :: sub_Covindex(num_allnodes, num_samplingnodes)
    real(dp) :: sub_samples(num_allnodes, num_samplingnodes)
    real(dp) :: sub_CovVector(num_allnodes, num_samplingnodes) 

    real(dp) :: Tsquare_add(num_samplingnodes)
    real(dp) :: Tsquare_add_m(num_allnodes)
    real(dp) :: cond_mean, cond_var, lower_bound, upper_bound
    
    real(dp) :: temp_vector(num_allnodes)

    integer :: candidates(num_allnodes) 
    integer :: cand_count               
    integer :: rand_picker(1)           
    real(dp) :: current_max_val    
    
    integer :: obs_list(num_samplingnodes)
    real(dp) :: obs_vals(num_samplingnodes)
    real(dp) :: pre_Vec_iS(num_samplingnodes)
    real(dp) :: Mat_SS(num_samplingnodes, num_samplingnodes)
    real(dp) :: Global_Sample_Map(num_allnodes) 
    
    logical  :: is_selected(num_allnodes)
    real(dp) :: Mat_SS_step(num_samplingnodes, num_samplingnodes)
    real(dp) :: pre_Vec_iS_step(num_samplingnodes)
    real(dp) :: obs_vals_step(num_samplingnodes)

    if (.not. data_loaded_CDS) then
        print *, "Error: CDS data not initialized."
        stop
    end if

    ! --- Initialization ---
    sub_Covnum = 0
    sampling_iperm = 0
    Tsquare_add = 0.0d0
    Tsquare_add_m = 0.0d0

    ! =================================================================
    ! STEP A: Update CUSUM for Sampled Nodes
    ! =================================================================
    sampling_iperm = 0
    Global_Sample_Map = 0.0d0
    
    do i = 1, num_samplingnodes
        idx = sampling_index(i)
        sampling_iperm(idx) = 1
        Global_Sample_Map(idx) = OnlineSample(i)
        
        CDS_cusumdouble_statistic(idx,1) = max(CDS_cusumdouble_statistic(idx,1) + &
                                           (miu_min_CDS * OnlineSample(i) - cuparame_CDS), 0.0d0)
        CDS_cusumdouble_statistic(idx,2) = max(CDS_cusumdouble_statistic(idx,2) + &
                                           (-miu_min_CDS * OnlineSample(i) - cuparame_CDS), 0.0d0)
        CDS_cusum_statistic(idx) = max(CDS_cusumdouble_statistic(idx,1), CDS_cusumdouble_statistic(idx,2))
    end do

    ! =================================================================
    ! STEP B: Inference for Non-sampled Nodes (Pull Method)
    ! =================================================================
    do i = 1, num_allnodes
        if (sampling_iperm(i) == 1) cycle 
        
        cond_mean = 0.0d0
        cond_var  = CDS_Cov_diag(i)
        k = 0
        
        do j = 1, num_varnear_CDS
            idx = CDS_Cov_NeighborId(i, j)
            if (idx > 0 .and. sampling_iperm(idx) == 1) then
                k = k + 1
                obs_list(k)   = idx
                obs_vals(k)   = Global_Sample_Map(idx)
                pre_Vec_iS(k) = CDS_Cov_NeighborVal(i, j) 
            end if
        end do

        if (k > 0) then
            do ii = 1, k
                do jj = 1, k
                    Mat_SS(ii, jj) = Cov_dense_global(obs_list(ii), obs_list(jj))
                end do
            end do
            
            call internal_conditional_solve( &
                target_var = CDS_Cov_diag(i), & 
                k          = k, &
                Mat_SS     = Mat_SS(1:k, 1:k), &
                pre_Vec_iS = pre_Vec_iS(1:k), & 
                obs_vals   = obs_vals(1:k), &
                out_mean   = cond_mean, &
                out_var    = cond_var &
            )
        end if
        

        !CDS_cusumdouble_statistic(i,1) = max(CDS_cusumdouble_statistic(i,1) + &
                                         !(miu_min_CDS * cond_mean - cuparame_CDS), 0.0d0)
        !CDS_cusumdouble_statistic(i,2) = max(CDS_cusumdouble_statistic(i,2) + &
                                         !(-miu_min_CDS * cond_mean - cuparame_CDS), 0.0d0)
        
         lower_bound = cond_mean - Phi_CDS * sqrt(max(cond_var, 0.0d0))
         upper_bound = cond_mean + Phi_CDS * sqrt(max(cond_var, 0.0d0))
         CDS_cusumdouble_statistic(i,1) = max(CDS_cusumdouble_statistic(i,1) + (miu_min_CDS * upper_bound - cuparame_CDS), 0.0d0)
         CDS_cusumdouble_statistic(i,2) = max(CDS_cusumdouble_statistic(i,2) + (-miu_min_CDS * lower_bound - cuparame_CDS), 0.0d0)

        CDS_cusum_statistic(i) = max(CDS_cusumdouble_statistic(i,1), CDS_cusumdouble_statistic(i,2))
    end do

    charting_statistic = MAXVAL(CDS_cusum_statistic)
! =================================================================
    ! STEP C: Calculate Initial T^2 Scores for all nodes
    ! =================================================================
    is_selected = .false.
    do i = 1, num_allnodes
        ! Assuming Kernel_Cov_diag contains the variance (inclusive of noise)
        Tsquare_add_m(i) = (CDS_cusum_statistic(i)**2) / Cov_dense_global(i,i)
    end do

    ! =================================================================
    ! STEP E: Select Next Sampling Nodes (Greedy Multivariate T^2 Maximization)
    ! =================================================================
    do i = 1, num_samplingnodes
        
        ! --- 1. Find the node with the maximum incremental T^2 score ---
        current_max_val = -huge(1.0_dp)
        cand_count = 0
        
        ! Single pass to find max and collect ties
        do k = 1, num_allnodes
            if (is_selected(k)) cycle ! Skip already picked nodes
            
            if (Tsquare_add_m(k) > current_max_val + 1.0d-12) then
                current_max_val = Tsquare_add_m(k)
                cand_count = 1
                candidates(1) = k
            else if (abs(Tsquare_add_m(k) - current_max_val) <= 1.0d-12) then
                cand_count = cand_count + 1
                candidates(cand_count) = k
            end if
        end do
        
        ! Random tie-breaking if multiple nodes have the exact same max score
        if (cand_count > 1) then
            call RNUND(cand_count, rand_picker) 
            max_id = candidates(rand_picker(1))
        else
            max_id = candidates(1)
        end if

        ! --- 2. Register the selected node ---
        sampling_index(i) = max_id
        is_selected(max_id) = .true.
        obs_vals_step(i) = CDS_cusum_statistic(max_id)
        
        ! --- 3. Build the Covariance Matrix of currently selected nodes (Mat_SS) ---
        ! This is the same for all unselected nodes, so we build it ONCE per step.
        do ii = 1, i
            do jj = 1, i
                Mat_SS_step(ii, jj) = Cov_dense_global(sampling_index(ii), sampling_index(jj))
            end do
        end do
        
        ! --- 4. Update the incremental T^2 score for ALL remaining candidates ---
        ! Math: Delta T^2 = (cusum - cond_mean)^2 / cond_var
        do k = 1, num_allnodes
            if (is_selected(k)) cycle
            
            ! Extract covariance vector between node 'k' and the 'i' selected nodes
            do ii = 1, i
                pre_Vec_iS_step(ii) = Cov_dense_global(k, sampling_index(ii))
            end do
            
            call internal_conditional_solve( &
                target_var   = Cov_dense_global(k,k), &
                k            = i, &
                Mat_SS       = Mat_SS_step(1:i, 1:i), &
                pre_Vec_iS   = pre_Vec_iS_step(1:i), &
                obs_vals     = obs_vals_step(1:i), &
                out_mean     = cond_mean, &
                out_var      = cond_var &
            )
                                   
            Tsquare_add_m(k) = ((CDS_cusum_statistic(k) - cond_mean)**2) / max(cond_var, 1.0d-16)
        end do
        
    end do


    contains

    subroutine internal_conditional_solve(target_var, k, Mat_SS, pre_Vec_iS, obs_vals, out_mean, out_var)
    implicit none
    integer, intent(in)   :: k
    real(dp), intent(in)  :: target_var         ! Target variable marginal variance
    real(dp), intent(in)  :: Mat_SS(k, k)       ! Covariance matrix of the k observations
    real(dp), intent(in)  :: pre_Vec_iS(k)      ! Covariances between target and k observations
    real(dp), intent(in)  :: obs_vals(k)        ! The k observation values
    real(dp), intent(out) :: out_mean, out_var  ! Output: conditional mean and variance

    real(dp) :: L(k, k)
    real(dp) :: y(k), w(k)
    integer  :: ii, jj

    ! 1. Cholesky Decomposition (L * L^T = Mat_SS)
    L = 0.0d0
    do ii = 1, k
        do jj = 1, ii - 1
            L(ii, jj) = (Mat_SS(ii, jj) - dot_product(L(ii, 1:jj-1), L(jj, 1:jj-1))) / L(jj, jj)
        end do
        ! Use max() to prevent negative square root due to floating point inaccuracies
        L(ii, ii) = sqrt(max(Mat_SS(ii, ii) - dot_product(L(ii, 1:ii-1), L(ii, 1:ii-1)), 1.0d-12))
    end do

    ! 2. Forward substitution (Solve L * y = pre_Vec_iS)
    do ii = 1, k
        y(ii) = (pre_Vec_iS(ii) - dot_product(L(ii, 1:ii-1), y(1:ii-1))) / L(ii, ii)
    end do

    ! 3. Compute Conditional Variance
    out_var = max(target_var - dot_product(y, y), 0.0d0)

    ! 4. Backward substitution (Solve L^T * w = y)
    do ii = k, 1, -1
        w(ii) = (y(ii) - dot_product(L(ii+1:k, ii), w(ii+1:k))) / L(ii, ii)
    end do

    ! 5. Compute Conditional Mean
    out_mean = dot_product(w, obs_vals)

end subroutine internal_conditional_solve

end subroutine CDS
end module CDS_mod
