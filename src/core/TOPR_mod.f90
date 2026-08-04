module TOPR_mod
    use GlobalSettings_mod
    use RNNOF_INT
    use RNUND_INT
    use RNBIN_INT
    use SVRGP_INT
    use linear_operators
    use utils_mod
    use SparseMatrix_mod
    use RNNOR_INT
    implicit none
    public 
    save   
    
    !---------------------- TOPR Parameters --------------------- 
    real(dp) :: miu_min_TO=1.0d0, delta_TO=0.00001d0
    integer  :: Topr_TO=1 
    real(dp) :: cuparame_TO=0.5d0 ! (1.0^2)/2.0
    
    logical  :: initialized_TOPR = .false.
    logical  :: data_loaded_TOPR = .false.

    !---------------------- INTERNAL DYNAMIC STATE ------------------
    ! CUSUM Statistics (Allocated once, Updated iteratively)
    real(dp), allocatable :: TOPR_cusumdouble_statistic(:,:) ! (num_allnodes, 2)
    real(dp), allocatable :: TOPR_cusum_statistic(:)         ! (num_allnodes)
    
    contains
    
!==================================================================
!  Subroutine: set_TOPRparams
!  Purpose: Initialize or update TOPR-specific parameters
!==================================================================
subroutine set_TOPRparams(miu_min_TO_in, delta_TO_in, Topr_TO_in)
    implicit none
    real(dp), intent(in), optional :: miu_min_TO_in, delta_TO_in
    integer, intent(in), optional :: Topr_TO_in

    if (present(miu_min_TO_in)) miu_min_TO = miu_min_TO_in
    if (present(delta_TO_in))   delta_TO   = delta_TO_in
    if (present(Topr_TO_in))    Topr_TO    = Topr_TO_in

    cuparame_TO = (miu_min_TO**2.0d0) / 2.0d0
    initialized_TOPR = .true.
end subroutine set_TOPRparams

!==================================================================
!  Subroutine: Init_TOPR_Data
!  Purpose: 
!     1. Allocate internal CUSUM state arrays based on num_allnodes.
!     2. Reset state to initial values.
!==================================================================
subroutine Init_TOPR_Data()
    implicit none
    
    ! 1. Allocate Dynamic State
    if (allocated(TOPR_cusumdouble_statistic)) deallocate(TOPR_cusumdouble_statistic)
    allocate(TOPR_cusumdouble_statistic(num_allnodes, 2))

    if (allocated(TOPR_cusum_statistic)) deallocate(TOPR_cusum_statistic)
    allocate(TOPR_cusum_statistic(num_allnodes))

    ! 2. Initialize State
    call Reset_TOPR_State()

    data_loaded_TOPR = .true.

end subroutine Init_TOPR_Data

!==================================================================
!  Subroutine: Reset_TOPR_State
!  Purpose: Reset CUSUM statistics to zero
!==================================================================
subroutine Reset_TOPR_State()
    implicit none
    if (allocated(TOPR_cusumdouble_statistic)) TOPR_cusumdouble_statistic = 0.0d0
    if (allocated(TOPR_cusum_statistic))       TOPR_cusum_statistic       = 0.0d0
end subroutine Reset_TOPR_State

!==================================================================
!  Subroutine: Clean_TOPR
!  Purpose: Free all internal memory
!==================================================================
subroutine Clean_TOPR()
    implicit none
    if (allocated(TOPR_cusumdouble_statistic)) deallocate(TOPR_cusumdouble_statistic)
    if (allocated(TOPR_cusum_statistic))       deallocate(TOPR_cusum_statistic)
    data_loaded_TOPR = .false.
end subroutine Clean_TOPR

!==================================================================
!  Subroutine: get_TOPRparams
!  Purpose: Print parameter status
!==================================================================
subroutine get_TOPRparams(fid)
    implicit none
    integer, intent(in) :: fid

    write(fid,'(A)') "==============================================="
    write(fid,'(A)') "           TOPR PARAMETER CONFIGURATION"
    write(fid,'(A)') "==============================================="
    write(fid,'(A,I8)') "num_samplingnodes   = ", num_samplingnodes
    write(fid,'(A,F12.6)') "miu_min_TO          = ", miu_min_TO
    write(fid,'(A,F12.6)') "delta_TO            = ", delta_TO
    write(fid,'(A,I8)')   "Topr_TO             = ", Topr_TO
    write(fid,'(A,L1)')   "data_loaded_TOPR    = ", data_loaded_TOPR
    write(fid,'(A)') "==============================================="
end subroutine get_TOPRparams

!======================================================================
!  Subroutine: TOPR
!  Purpose:
!    Select top-r sampling nodes based on current CUSUM statistics.
!    Uses internal state variables for statistics.
!======================================================================
subroutine TOPR(OnlineSample, sampling_index, charting_statistic)
    implicit none

    !===========================================================
    ! Input parameters
    !===========================================================
    real(dp), intent(in) :: OnlineSample(num_samplingnodes)

    !===========================================================
    ! In/out parameters
    ! sampling_index: Input = Current nodes, Output = Next nodes
    !===========================================================
    integer, intent(inout) :: sampling_index(num_samplingnodes)

    !===========================================================
    ! Output parameter
    !===========================================================
    real(dp), intent(out) :: charting_statistic

    !===========================================================
    ! Local variables
    !===========================================================
    real(dp) :: temp_vector(num_allnodes)
    integer  :: sampling_iperm(num_allnodes)
    integer  :: i, j, idx

    if (.not. data_loaded_TOPR) then
        print *, "Error: TOPR data not initialized."
        stop
    end if

    !===========================================================
    ! Step 1. Update CUSUM for sampled nodes (signal contribution)
    !===========================================================
    sampling_iperm = 0
    do j = 1, num_samplingnodes
        idx = sampling_index(j)
        
        TOPR_cusumdouble_statistic(idx,1) = max(TOPR_cusumdouble_statistic(idx,1) + &
                                                (miu_min_TO * OnlineSample(j) - cuparame_TO), 0.0d0)
        TOPR_cusumdouble_statistic(idx,2) = max(TOPR_cusumdouble_statistic(idx,2) + &
                                                (-miu_min_TO * OnlineSample(j) - cuparame_TO), 0.0d0)
                                                
        TOPR_cusum_statistic(idx) = max(TOPR_cusumdouble_statistic(idx,1), TOPR_cusumdouble_statistic(idx,2))
        
        sampling_iperm(idx) = 1
    end do

    !===========================================================
    ! Step 2. Update CUSUM for non-sampled nodes (drift increment)
    !===========================================================
    do j = 1, num_allnodes
        if (sampling_iperm(j) == 0) then
            ! Logic correction: Use 'j' as index, not 'sampling_iperm(j)' (which is 0)
            TOPR_cusumdouble_statistic(j,1) = TOPR_cusumdouble_statistic(j,1) + delta_TO
            TOPR_cusumdouble_statistic(j,2) = TOPR_cusumdouble_statistic(j,2) + delta_TO
            
            TOPR_cusum_statistic(j) = max(TOPR_cusumdouble_statistic(j,1), TOPR_cusumdouble_statistic(j,2))
        end if
    end do

    !===========================================================
    ! Step 3. Compute charting statistic = sum of top-r CUSUM values
    !===========================================================
    temp_vector = TOPR_cusum_statistic
    
    ! Partial sort to find top K values
    call partial_quickselect(values=temp_vector, num_all=num_allnodes, &
                             top_k=Topr_TO, order_index=-1, &
                             values_sub=temp_vector(1:Topr_TO))
                             
    charting_statistic = sum(temp_vector(1:Topr_TO))

    !===========================================================
    ! Step 4. Select Next Sampling Nodes (Sort all by CUSUM)
    !         (Add small random perturbation to break ties)
    !===========================================================
    call RNNOR(temp_vector)
    temp_vector = TOPR_cusum_statistic + eps_perturb * temp_vector
    
    call partial_quickselect(values=temp_vector, num_all=num_allnodes, &
                             top_k=num_samplingnodes, order_index=-1, &
                             iperm_sub=sampling_index)

end subroutine TOPR
                
end module TOPR_mod