module TSSPR_mod
    use GlobalSettings_mod
    use RNNOF_INT
    use RNUND_INT
    use RNBIN_INT
    use RNUN_INT   ! IMSL uniform random-number interface
    use SVRGP_INT
    use linear_operators
    use utils_mod
    use SparseMatrix_mod
    use RNNOR_INT
    implicit none
    public 
    save
    
    !---------------------- TSSPR Parameters --------------------- 
    real(dp) :: miu_min_TS=1.0d0
    integer  :: Topr_TS=1 
    real(dp) :: cuparame_TS=0.5d0 ! (1.0^2)/2.0
    
    logical  :: initialized_TSSPR = .false.  
    logical  :: data_loaded_TSSPR = .false.

    !---------------------- INTERNAL DYNAMIC STATE ------------------
    ! L_statistic: Likelihood product (starts at 1.0)
    ! R_statistic: Accumulation statistic (starts at 0.0)
    real(dp), allocatable :: TSSPR_L_statistic(:,:) ! (num_allnodes, 2)
    real(dp), allocatable :: TSSPR_R_statistic(:,:) ! (num_allnodes, 2)
    
    contains
    
!==================================================================
!  Subroutine: set_TSSPRparams
!  Purpose: Initialize or update TSSPR-specific parameters
!==================================================================
subroutine set_TSSPRparams(miu_min_TS_in, Topr_TS_in)
    implicit none
    real(dp), intent(in), optional :: miu_min_TS_in
    integer, intent(in), optional :: Topr_TS_in

    if (present(miu_min_TS_in)) miu_min_TS = miu_min_TS_in
    if (present(Topr_TS_in))    Topr_TS    = Topr_TS_in

    cuparame_TS = (miu_min_TS**2.0d0) / 2.0d0
    initialized_TSSPR = .true.
end subroutine set_TSSPRparams

!==================================================================
!  Subroutine: Init_TSSPR_Data
!  Purpose: 
!     1. Allocate internal state arrays based on num_allnodes.
!     2. Reset state to initial values.
!==================================================================
subroutine Init_TSSPR_Data()
    implicit none
    
    ! 1. Allocate Dynamic State
    if (allocated(TSSPR_L_statistic)) deallocate(TSSPR_L_statistic)
    allocate(TSSPR_L_statistic(num_allnodes, 2))

    if (allocated(TSSPR_R_statistic)) deallocate(TSSPR_R_statistic)
    allocate(TSSPR_R_statistic(num_allnodes, 2))

    ! 2. Initialize State
    call Reset_TSSPR_State()

    data_loaded_TSSPR = .true.

end subroutine Init_TSSPR_Data

!==================================================================
!  Subroutine: Reset_TSSPR_State
!  Purpose: Reset statistics to initial values.
!           Note: L (Product) -> 1.0, R (Sum) -> 0.0
!==================================================================
subroutine Reset_TSSPR_State()
    implicit none
    if (allocated(TSSPR_R_statistic)) TSSPR_R_statistic = 0.0d0
    if (allocated(TSSPR_L_statistic)) TSSPR_L_statistic = 1.0d0
end subroutine Reset_TSSPR_State

!==================================================================
!  Subroutine: Clean_TSSPR
!  Purpose: Free all internal memory
!==================================================================
subroutine Clean_TSSPR()
    implicit none
    if (allocated(TSSPR_L_statistic)) deallocate(TSSPR_L_statistic)
    if (allocated(TSSPR_R_statistic)) deallocate(TSSPR_R_statistic)
    data_loaded_TSSPR = .false.
end subroutine Clean_TSSPR

!==================================================================
!  Subroutine: get_TSSPRparams
!  Purpose: Print parameter status
!==================================================================
subroutine get_TSSPRparams(fid)
    implicit none
    integer, intent(in) :: fid

    write(fid,'(A)') "==============================================="
    write(fid,'(A)') "           TSSPR PARAMETER CONFIGURATION"
    write(fid,'(A)') "==============================================="
    write(fid,'(A,I8)') "num_samplingnodes   = ", num_samplingnodes
    write(fid,'(A,F12.6)') "miu_min_TS          = ", miu_min_TS
    write(fid,'(A,I8)')   "Topr_TS             = ", Topr_TS
    write(fid,'(A,F12.6)') "cuparame_TS         = ", cuparame_TS
    write(fid,'(A,L1)')   "data_loaded_TSSPR   = ", data_loaded_TSSPR
    write(fid,'(A)') "==============================================="
end subroutine get_TSSPRparams   

!======================================================================
!  Subroutine: TSSPR
!  State variables are stored at module scope.
!======================================================================
subroutine TSSPR(OnlineSample, sampling_index, charting_statistic)
    implicit none

    !==============================================================
    ! Input / Output Declarations
    !==============================================================
    real(dp), intent(in)  :: OnlineSample(num_samplingnodes)
    
    ! sampling_index: Input = Current nodes, Output = Next nodes
    integer, intent(inout) :: sampling_index(num_samplingnodes) 
    real(dp), intent(out)   :: charting_statistic                 

    !==============================================================
    ! Local Variables
    !==============================================================
    real(dp) :: Sample_statistic(num_allnodes)
    real(dp) :: TS_statistic(num_allnodes)       
    real(dp) :: temp_vector(num_allnodes)        
    real(dp) :: Temvalue                        
    integer  :: i, idx

    if (.not. data_loaded_TSSPR) then
        print *, "Error: TSSPR data not initialized."
        stop
    end if

    !==============================================================
    ! Step 1: Update L_statistic and R_statistic for sampled nodes
    !         (Uses internal TSSPR_L_statistic, TSSPR_R_statistic)
    !==============================================================
    ! Global update for R (Shiryaev-Roberts style accumulation)
    TSSPR_R_statistic = TSSPR_R_statistic + 1.0d0
    
    do i = 1, num_samplingnodes
        idx = sampling_index(i)
        
        ! Positive direction update
        Temvalue = exp(miu_min_TS * OnlineSample(i) - cuparame_TS)
        TSSPR_R_statistic(idx,1) = TSSPR_R_statistic(idx,1) * Temvalue
        TSSPR_L_statistic(idx,1) = TSSPR_L_statistic(idx,1) * Temvalue

        ! Negative direction update
        Temvalue = exp(-miu_min_TS * OnlineSample(i) - cuparame_TS)
        TSSPR_R_statistic(idx,2) = TSSPR_R_statistic(idx,2) * Temvalue
        TSSPR_L_statistic(idx,2) = TSSPR_L_statistic(idx,2) * Temvalue
    end do

    !==============================================================
    ! Step 2: Compute charting statistic = sum of top-r R_statistic maxima
    !==============================================================
    temp_vector = maxval(TSSPR_R_statistic, dim=2)
    
    call partial_quickselect(values=temp_vector, num_all=num_allnodes, &
                             top_k=Topr_TS, order_index=-1, &
                             values_sub=temp_vector(1:Topr_TS))
                             
    charting_statistic = sum(temp_vector(1:Topr_TS))

    !==============================================================
    ! Step 3: Generate random uniform samples for stochastic update
    !==============================================================
    call RNUN(TS_statistic, num_allnodes)

    ! Combine left/right statistic branches for each node
    do i = 1, num_allnodes
        Sample_statistic(i) = max( &
            TSSPR_R_statistic(i,1) + TSSPR_L_statistic(i,1) * TS_statistic(i), &
            TSSPR_R_statistic(i,2) + TSSPR_L_statistic(i,2) * TS_statistic(i) )
    end do

    !==============================================================
    ! Step 4: Add perturbation noise and reselect top nodes
    !==============================================================
    call RNNOR(temp_vector)
    temp_vector = Sample_statistic + eps_perturb * temp_vector
    
    call partial_quickselect(values=temp_vector, num_all=num_allnodes, &
                             top_k=num_samplingnodes, order_index=-1, &
                             iperm_sub=sampling_index)

end subroutine TSSPR

end module TSSPR_mod
