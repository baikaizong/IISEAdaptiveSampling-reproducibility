module TSBSS_mod
    use GlobalSettings_mod
    use RNNOF_INT
    use RNUND_INT
    use RNBIN_INT
    use SVRGP_INT
    use RNNOR_INT
    use linear_operators
    use utils_mod
    use SparseMatrix_mod
    implicit none
    public 
    save
    
    !---------------------- TSBSS Parameters ------------------
    real(dp) :: lambda0_TS=0.9d0, w0_TS=0.01d0, sigma0_TS=5.0d0
    real(dp) :: sigmae_TS=1.0d0, sigmab_TS=3.0d0
    integer  :: niter_TS=20
    
    logical  :: initialized_TSBSS = .false.
    logical  :: data_loaded_TSBSS = .false.

    !---------------------- INTERNAL STATIC DATA --------------------
    ! Large sparse basis matrices (Assigned once, Read-only)
    type(Bbasis_sparse_type) :: TSBSS_B0_all
    type(Bbasis_sparse_type) :: TSBSS_B1_all

    !---------------------- INTERNAL DYNAMIC STATE ------------------
    ! Sufficient statistics updated iteratively
    real(dp), allocatable :: TSBSS_XB(:,:)   ! (knot1_square_TS, 1)
    real(dp), allocatable :: TSBSS_BB(:,:)   ! (knot1_square_TS, 1)
    real(dp), allocatable :: TSBSS_BBjk(:,:) ! (knot1_square_TS, knot1_square_TS)

    contains

!==================================================================
!  Subroutine: set_TSBSSparams
!  Purpose: Initialize or update TSBSS scalar parameters
!==================================================================
subroutine set_TSBSSparams( lambda0_TS_in, w0_TS_in, sigma0_TS_in, &
    sigmae_TS_in, sigmab_TS_in, niter_TS_in)
    implicit none
    real(dp), intent(in), optional :: lambda0_TS_in, w0_TS_in, sigma0_TS_in
    real(dp), intent(in), optional :: sigmae_TS_in, sigmab_TS_in
    integer, intent(in), optional :: niter_TS_in

    if (present(lambda0_TS_in)) lambda0_TS = lambda0_TS_in
    if (present(w0_TS_in))      w0_TS      = w0_TS_in
    if (present(sigma0_TS_in))  sigma0_TS  = sigma0_TS_in
    if (present(sigmae_TS_in))  sigmae_TS  = sigmae_TS_in
    if (present(sigmab_TS_in))  sigmab_TS  = sigmab_TS_in
    if (present(niter_TS_in))   niter_TS   = niter_TS_in

    initialized_TSBSS = .true.
end subroutine set_TSBSSparams

!==================================================================
!  Subroutine: Init_TSBSS_Data
!  Purpose: 
!     1. Copy large sparse matrices B0 and B1 internally ONE TIME.
!     2. Allocate internal state matrices (XB, BB, BBjk).
!     3. Reset state to initial conditions.
!==================================================================
subroutine Init_TSBSS_Data(B0_in, B1_in)
    implicit none
    type(Bbasis_sparse_type), intent(in), optional :: B0_in
    type(Bbasis_sparse_type), intent(in) :: B1_in
    call copy_Bbasis_sparse(B1_in, TSBSS_B1_all)

    if (present(B0_in)) then
        call copy_Bbasis_sparse(B0_in, TSBSS_B0_all)
    end if

    ! 2. Allocate Dynamic State
    if (allocated(TSBSS_XB)) deallocate(TSBSS_XB)
    if (allocated(TSBSS_BB)) deallocate(TSBSS_BB)
    if (allocated(TSBSS_BBjk)) deallocate(TSBSS_BBjk)

    allocate(TSBSS_XB(knot1_square_TS, 1))
    allocate(TSBSS_BB(knot1_square_TS, 1))
    allocate(TSBSS_BBjk(knot1_square_TS, knot1_square_TS))

    ! 3. Initialize
    call Reset_TSBSS_State()
    
    data_loaded_TSBSS = .true.

end subroutine Init_TSBSS_Data

!==================================================================
!  Subroutine: Reset_TSBSS_State
!  Purpose: Reset recursive statistics (XB, BB, BBjk) to zero/initial
!==================================================================
subroutine Reset_TSBSS_State()
    implicit none
    if (allocated(TSBSS_XB)) TSBSS_XB = 0.0d0
    if (allocated(TSBSS_BB)) TSBSS_BB = 0.0d0
    if (allocated(TSBSS_BBjk)) TSBSS_BBjk = 0.0d0
end subroutine Reset_TSBSS_State

!==================================================================
!  Subroutine: Clean_TSBSS
!  Purpose: Free all internal memory
!==================================================================
subroutine Clean_TSBSS()
    implicit none
    call free_Bbasis_sparse(TSBSS_B0_all)
    call free_Bbasis_sparse(TSBSS_B1_all)
    
    if (allocated(TSBSS_XB)) deallocate(TSBSS_XB)
    if (allocated(TSBSS_BB)) deallocate(TSBSS_BB)
    if (allocated(TSBSS_BBjk)) deallocate(TSBSS_BBjk)
    
    data_loaded_TSBSS = .false.
end subroutine Clean_TSBSS

!==================================================================
!  Subroutine: get_TSBSSparams
!  Purpose: Print parameters
!==================================================================
subroutine get_TSBSSparams(fid)
    implicit none
    integer, intent(in) :: fid  

    write(fid,'(A)') "==============================================="
    write(fid,'(A)') "           TSBSS PARAMETER CONFIGURATION"
    write(fid,'(A)') "==============================================="
    write(fid,'(A,I8)') "num_samplingnodes   = ", num_samplingnodes
    write(fid,'(A,F12.6)') "lambda0_TS          = ", lambda0_TS
    write(fid,'(A,F12.6)') "w0_TS               = ", w0_TS
    write(fid,'(A,F12.6)') "sigma0_TS           = ", sigma0_TS
    write(fid,'(A,F12.6)') "sigmae_TS           = ", sigmae_TS
    write(fid,'(A,F12.6)') "sigmab_TS           = ", sigmab_TS
    write(fid,'(A,I8)') "niter_TS            = ", niter_TS
    write(fid,'(A,L1)') "data_loaded_TSBSS   = ", data_loaded_TSBSS
    write(fid,'(A)') "==============================================="
end subroutine get_TSBSSparams  

!======================================================================
! Subroutine: TSBSS
! Description: Temporal-Spatial Bayesian Sparse Sampling update.
!              Updates posterior statistics based on new samples and 
!              determines the next best sampling locations.
!======================================================================
subroutine TSBSS(back_index, OnlineSample, &
                 sampling_index, charting_statistic)

    implicit none

    !-------------------- Input Arguments ----------------------------
    integer, intent(in)           :: back_index
    real(dp), intent(in)          :: OnlineSample(num_samplingnodes)
    
    !-------------------- Input / Output Arguments -------------------
    ! sampling_index is updated to point to the NEXT set of nodes
    integer, intent(inout)        :: sampling_index(num_samplingnodes)
    real(dp), intent(out)         :: charting_statistic

    !------------------------------------------------------------
    !  Local Variables
    !------------------------------------------------------------
    integer :: i, k
    integer :: one_vec(1)

    ! Automatic Arrays
    real(dp) :: alpha(knot1_square_TS), mu(knot1_square_TS), sigma_2(knot1_square_TS)
    real(dp) :: mu_tmp(knot1_square_TS), subsample(knot1_square_TS), XB_last(knot1_square_TS)
    
    real(dp) :: Tem_three(knot1_square_TS)     
    real(dp) :: Tem_four            
    
    real(dp) :: theta0(knot0_square_TS), temp_k0_vec(knot0_square_TS)
    
    real(dp) :: Enoise(num_allnodes)
    real(dp) :: subnoisesample(num_allnodes)
    real(dp) :: sampling_statistic(num_allnodes)
    real(dp) :: Tem_one_global(num_allnodes)
    real(dp) :: temp_vector(num_allnodes)
    
    real(dp) :: B1mu(num_samplingnodes), B0theta(num_samplingnodes)
    real(dp) :: terre(num_samplingnodes)
    
    real(dp), allocatable :: Tem_two(:,:)
    real(dp), allocatable :: Covb(:,:)

    ! Pre-computed constants
    real(dp) :: sigmae_sq, sigmae_inv_sq, sigma0_sq, inv_sigma0_sq
    real(dp) :: const_log_term, one_max_val

    ! Sparse temporary structures
    type(Bbasis_sparse_type) :: B0_subsparse, B1_subsparse

    if (.not. data_loaded_TSBSS) then
        print *, "Error: TSBSS data not initialized. Call Init_TSBSS_Data first."
        stop
    end if

    !============================================================
    ! PRE-COMPUTATION & INIT
    !============================================================
    sigmae_sq     = sigmae_TS**2.0d0
    sigmae_inv_sq = 1.0d0 / sigmae_sq
    sigma0_sq     = sigma0_TS**2.0d0
    inv_sigma0_sq = 1.0d0 / sigma0_sq
    
    const_log_term = log(w0_TS/(1.0d0-w0_TS)) + 0.5d0

    ! 1. Extract Sparse Submatrices from Internal Static Data
    if (back_index == BACK_BSPLINE) then
         call build_sub_sparse(TSBSS_B0_all, sampling_index, B0_subsparse)
    end if
    call build_sub_sparse(TSBSS_B1_all, sampling_index, B1_subsparse)

    !============================================================
    ! STEP 1: Update Sufficient Statistics (BB, BBjk)
    !============================================================
    if (.not. allocated(Tem_two)) allocate(Tem_two(knot1_square_TS, knot1_square_TS))
    
    ! Calculates B1^T * B1 (Dense output)
    call Sparse_MatMul_AtA(B1_subsparse, C_dense=Tem_two)
    
    ! Vectorized updates
    do i = 1, knot1_square_TS
        Tem_three(i) = Tem_two(i,i)
    end do

    TSBSS_BB(:,1) = TSBSS_BB(:,1) * lambda0_TS + Tem_three(:)
    TSBSS_BBjk    = TSBSS_BBjk * lambda0_TS + Tem_two

    !============================================================
    ! STEP 2: Initialize Posterior Parameters
    !============================================================
    sigma_2 = 1.0d0 / ((TSBSS_BB(:,1) * sigmae_inv_sq) + inv_sigma0_sq)
    mu      = 0.0d0
    alpha   = 0.5d0
    mu_tmp  = 0.0d0 

    !============================================================
    ! STEP 3: Bayesian Update Loop (Updates XB, mu, alpha)
    !============================================================
    if (back_index == BACK_BSPLINE) then
        !--------------------------------------------------------
        ! CASE 1: With Background
        !--------------------------------------------------------
        if (.not. allocated(Covb)) allocate(Covb(knot0_square_TS, knot0_square_TS))
        
        call Sparse_MatMul_AtA(B0_subsparse, C_dense=Covb)
        
        Covb = sigmae_inv_sq * Covb
        do i = 1, knot0_square_TS
            Covb(i,i) = Covb(i,i) + (sigmab_TS**(-2.0d0))
        end do
        Covb = .i. Covb 

        XB_last = TSBSS_XB(:,1)

        do k = 1, niter_TS
            call spmv_Bbasis_sparse(B1_subsparse, mu_tmp, B1mu)

            terre = OnlineSample - B1mu
            call spmv_Bbasis_sparse_transpose(B0_subsparse, terre, temp_k0_vec)
            
            theta0 = sigmae_inv_sq * (Covb .x. temp_k0_vec)

            call spmv_Bbasis_sparse(B0_subsparse, theta0, B0theta)
            terre = OnlineSample - B0theta

            call spmv_Bbasis_sparse_transpose(B1_subsparse, terre, Tem_three)
            TSBSS_XB(:,1) = XB_last * lambda0_TS + Tem_three

            do i = 1, knot1_square_TS
                one_max_val = dot_product(TSBSS_BBjk(:, i), mu_tmp)
                
                mu(i) = (TSBSS_XB(i,1) - (one_max_val - TSBSS_BBjk(i,i)*alpha(i)*mu(i))) * &
                        sigma_2(i) * sigmae_inv_sq

                Tem_four = const_log_term - &
                           (sigma_2(i) + mu(i)**2.0d0)/(2.0d0*sigma0_sq) + &
                           log(sqrt(sigma_2(i))/sigma0_TS) + &
                           (mu(i)**2.0d0)/sigma_2(i) - &
                           TSBSS_BB(i,1)*((mu(i)**2.0d0) + sigma_2(i))/(2.0d0*sigmae_sq)
                
                alpha(i) = 1.0d0 / (1.0d0 + exp(-Tem_four))
                
                if (alpha(i) /= alpha(i) .or. alpha(i) >= 0.9999999d0) alpha(i) = 0.9999999d0
                
                mu_tmp(i) = mu(i) * alpha(i)
            end do
        end do

    else
        !--------------------------------------------------------
        ! CASE 2: Without Background
        !--------------------------------------------------------
        XB_last = TSBSS_XB(:,1)
        
        call spmv_Bbasis_sparse_transpose(B1_subsparse, OnlineSample, Tem_three)
        TSBSS_XB(:,1) = XB_last * lambda0_TS + Tem_three
        
        do k = 1, niter_TS
            do i = 1, knot1_square_TS
                one_max_val = dot_product(TSBSS_BBjk(:, i), mu_tmp)
                
                mu(i) = (TSBSS_XB(i,1) - (one_max_val - TSBSS_BBjk(i,i)*alpha(i)*mu(i))) * &
                        sigma_2(i) * sigmae_inv_sq

                Tem_four = const_log_term - &
                           (sigma_2(i) + mu(i)**2.0d0)/(2.0d0*sigma0_sq) + &
                           log(sqrt(sigma_2(i))/sigma0_TS) + &
                           (mu(i)**2.0d0)/sigma_2(i) - &
                           TSBSS_BB(i,1)*((mu(i)**2.0d0) + sigma_2(i))/(2.0d0*sigmae_sq)

                alpha(i) = 1.0d0 / (1.0d0 + exp(-Tem_four))
                
                if (alpha(i) /= alpha(i) .or. alpha(i) >= 0.9999999d0) alpha(i) = 0.9999999d0
                
                mu_tmp(i) = mu(i) * alpha(i)
            end do
        end do
    end if

    !============================================================
    ! STEP 4: Compute Charting Statistic
    !============================================================
    if (back_index == BACK_BSPLINE) then
         call spmv_Bbasis_sparse_transpose(B1_subsparse, terre, Tem_three)
    else
         call spmv_Bbasis_sparse_transpose(B1_subsparse, OnlineSample, Tem_three)
    end if

    do i = 1, knot1_square_TS
        XB_last(i) = dot_product(Tem_two(:, i), mu_tmp) 
    end do
    
    charting_statistic = 2.0d0 * dot_product(Tem_three, mu_tmp) - &
                         dot_product(mu_tmp, XB_last)

    !============================================================
    ! STEP 5: Posterior Predictive Sampling (Global)
    !============================================================
    do i = 1, knot1_square_TS
        call RNBIN(1, real(alpha(i)), one_vec)
        if (one_vec(1) == 1) then
            subsample(i) = RNNOF() 
            subsample(i)=subsample(i)* sqrt(sigma_2(i)) + mu(i)
        else
            subsample(i) = 0.0d0
        end if
    end do

    call RNNOR(Enoise)
    
    ! Use Internal Global B1
    call spmv_Bbasis_sparse(TSBSS_B1_all, subsample, Tem_one_global)
    subnoisesample = Tem_one_global + (Enoise * sigma_noise)

    !============================================================
    ! STEP 6: Expected Improvement & Selection
    !============================================================
    sampling_statistic = 0.0d0
    
    call spmv_Bbasis_sparse(TSBSS_B1_all, mu_tmp, Tem_one_global)

    where (abs(Tem_one_global) > 1.0d-12)
        sampling_statistic = 2.0d0 * Tem_one_global * subnoisesample - &
                             Tem_one_global**2.0d0
    end where

    !============================================================
    ! STEP 7: Select Next Sampling Locations
    !============================================================
    call RNNOR(temp_vector)
    temp_vector = sampling_statistic + 1.0d-10 * temp_vector
    
    call partial_quickselect(values=temp_vector, num_all=num_allnodes, &
                             top_k=num_samplingnodes, order_index=-1, &
                             iperm_sub=sampling_index)

    !============================================================
    ! CLEANUP LOCAL RESOURCES
    !============================================================
    if (back_index == BACK_BSPLINE) call free_Bbasis_sparse(B0_subsparse)
    call free_Bbasis_sparse(B1_subsparse)

    if (allocated(Tem_two)) deallocate(Tem_two)
    if (allocated(Covb)) deallocate(Covb)

end subroutine TSBSS

end module TSBSS_mod