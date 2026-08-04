module NAS_mod
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
    
    !---------------------- NAS parameters ---------------------
    real(dp) :: gone_NAS=0.000000000000001
    real(dp) :: gp_NAS=0.0001
    real(dp) :: delta_NAS=0.0001, lambda_NAS=0.002, lambda0_NAS=0.00001
    real(dp) :: k_allowance_NAS=0.01
    
    logical :: initialized_NAS = .false.
    logical :: data_loaded_NAS = .false.

    !---------------------- INTERNAL DYNAMIC STATE ------------------
    ! Encapsulated state variables to avoid passing them as arguments
    real(dp), allocatable :: NAS_Smin_one(:), NAS_Smin_two(:)
    real(dp), allocatable :: NAS_Smax_one(:), NAS_Smax_two(:)
    
    contains
    
!==================================================================
!  Subroutine: set_NASparams
!  Purpose: Initialize or update NAS-specific parameters
!==================================================================
subroutine set_NASparams(gp_NAS_in, gone_NAS_in, delta_NAS_in, &
                         lambda0_NAS_in, lambda_NAS_in, k_allowance_NAS_in)
    implicit none
    real(dp), intent(in), optional :: gp_NAS_in, gone_NAS_in, delta_NAS_in
    real(dp), intent(in), optional :: lambda0_NAS_in, lambda_NAS_in, k_allowance_NAS_in

    if (present(gp_NAS_in))          gp_NAS          = gp_NAS_in
    if (present(gone_NAS_in))        gone_NAS        = gone_NAS_in
    if (present(delta_NAS_in))       delta_NAS       = delta_NAS_in
    if (present(lambda0_NAS_in))     lambda0_NAS     = lambda0_NAS_in
    if (present(lambda_NAS_in))      lambda_NAS      = lambda_NAS_in
    if (present(k_allowance_NAS_in)) k_allowance_NAS = k_allowance_NAS_in

    initialized_NAS = .true.
end subroutine set_NASparams

!==================================================================
!  Subroutine: Init_NAS_Data
!  Purpose: Allocate internal state arrays based on global node count
!==================================================================
subroutine Init_NAS_Data()
    implicit none
    
    ! Check allocation and allocate
    if (allocated(NAS_Smin_one)) deallocate(NAS_Smin_one)
    allocate(NAS_Smin_one(num_allnodes_addone))

    if (allocated(NAS_Smin_two)) deallocate(NAS_Smin_two)
    allocate(NAS_Smin_two(num_allnodes_addone))

    if (allocated(NAS_Smax_one)) deallocate(NAS_Smax_one)
    allocate(NAS_Smax_one(num_allnodes_addone))

    if (allocated(NAS_Smax_two)) deallocate(NAS_Smax_two)
    allocate(NAS_Smax_two(num_allnodes_addone))

    ! Initialize values
    call Reset_NAS_State()

    data_loaded_NAS = .true.
end subroutine Init_NAS_Data

!==================================================================
!  Subroutine: Reset_NAS_State
!  Purpose: Reset statistics to initial zero state
!==================================================================
subroutine Reset_NAS_State()
    implicit none
    if (allocated(NAS_Smin_one)) NAS_Smin_one = 0.0d0
    if (allocated(NAS_Smin_two)) NAS_Smin_two = 0.0d0
    if (allocated(NAS_Smax_one)) NAS_Smax_one = 0.0d0
    if (allocated(NAS_Smax_two)) NAS_Smax_two = 0.0d0
end subroutine Reset_NAS_State

!==================================================================
!  Subroutine: Clean_NAS
!  Purpose: Free memory
!==================================================================
subroutine Clean_NAS()
    implicit none
    if (allocated(NAS_Smin_one)) deallocate(NAS_Smin_one)
    if (allocated(NAS_Smin_two)) deallocate(NAS_Smin_two)
    if (allocated(NAS_Smax_one)) deallocate(NAS_Smax_one)
    if (allocated(NAS_Smax_two)) deallocate(NAS_Smax_two)
    data_loaded_NAS = .false.
end subroutine Clean_NAS

!==================================================================
!  Subroutine: get_NASparams
!==================================================================
subroutine get_NASparams(fid)
    implicit none
    integer, intent(in) :: fid

    write(fid,'(A)') "==============================================="
    write(fid,'(A)') "             NAS PARAMETER CONFIGURATION"
    write(fid,'(A)') "==============================================="
    write(fid,'(A,I8)') "num_samplingnodes   = ", num_samplingnodes
    write(fid,'(A,F12.6)') "gp_NAS              = ", gp_NAS
    write(fid,'(A,F12.6)') "gone_NAS            = ", gone_NAS
    write(fid,'(A,F12.6)') "delta_NAS           = ", delta_NAS
    write(fid,'(A,F12.6)') "lambda_NAS          = ", lambda_NAS
    write(fid,'(A,F12.6)') "lambda0_NAS         = ", lambda0_NAS
    write(fid,'(A,F12.6)') "k_allowance_NAS     = ", k_allowance_NAS
    write(fid,'(A,L1)')   "data_loaded_NAS     = ", data_loaded_NAS
    write(fid,'(A)') "==============================================="
end subroutine get_NASparams

!==================================================================
!  Subroutine: NAS
!  Dynamic state is stored in module-level arrays.
!==================================================================
subroutine NAS(OnlineSample, sampling_index, charting_statistic)
    implicit none

    !==============================================================
    ! Input Parameters
    !==============================================================
    real(dp), intent(in) :: OnlineSample(num_samplingnodes)

    !==============================================================
    ! Input/Output Parameters
    ! sampling_index: Input = Current nodes, Output = Next nodes
    !==============================================================
    integer, intent(inout) :: sampling_index(num_samplingnodes)

    !==============================================================
    ! Output Parameters
    !==============================================================
    real(dp), intent(out) :: charting_statistic

    !==============================================================
    ! Local Variables
    !==============================================================
    ! Combined arrays including the dummy node
    real(dp) :: tOnlineSample(num_samplingnodes_addone)
    integer  :: tsampling_index(num_samplingnodes_addone)
    
    ! Rank update vector (reused)
    real(dp) :: rankindex_s(num_allnodes_addone)
    
    ! Temporary calculation variables
    real(dp) :: S_one(num_allnodes) 
    real(dp) :: temp_vector(num_allnodes)
    real(dp) :: C_index, test_min, test_max, scale_factor
    integer  :: i, min_index, max_index, global_idx

    if (.not. data_loaded_NAS) then
        print *, "Error: NAS data not initialized. Call Init_NAS_Data() first."
        stop
    end if

    !==============================================================
    ! Step 1: Initialize Extended Sample and Index Arrays
    !==============================================================
    tOnlineSample(1:num_samplingnodes) = OnlineSample
    tOnlineSample(num_samplingnodes_addone) = 0.0d0
    
    tsampling_index(1:num_samplingnodes) = sampling_index
    tsampling_index(num_samplingnodes_addone) = num_allnodes_addone

    !==============================================================
    ! Step 2: Update Minimum Statistics (Smin)
    !==============================================================
    ! 2.1 Find location of the minimum value
    min_index = minloc(tOnlineSample, dim=1)
    
    ! 2.2 Construct Rank Index Vector
    rankindex_s = delta_NAS  
    rankindex_s(sampling_index) = 0.0d0        
    rankindex_s(num_allnodes_addone) = 0.0d0   

    global_idx = tsampling_index(min_index)

    if (global_idx == num_allnodes_addone) then
        rankindex_s(num_allnodes_addone) = lambda0_NAS
    else
        rankindex_s(global_idx) = lambda_NAS
    end if

    ! 2.3 Update Smin vectors (Using Internal State)
    NAS_Smin_one = NAS_Smin_one + rankindex_s
    NAS_Smin_two(1:num_allnodes) = NAS_Smin_two(1:num_allnodes) + gone_NAS
    NAS_Smin_two(num_allnodes_addone) = NAS_Smin_two(num_allnodes_addone) + gp_NAS

    ! 2.4 Compute Statistic C_index (Chi-square like)
    C_index = sum(((NAS_Smin_one - NAS_Smin_two)**2.0d0) / NAS_Smin_two)

    ! 2.5 Normalization / Reset
    if (C_index <= k_allowance_NAS) then
        NAS_Smin_one = 0.0d0
        NAS_Smin_two = 0.0d0
        test_min = 0.0d0
    else
        scale_factor = (C_index - k_allowance_NAS) / C_index
        NAS_Smin_one = NAS_Smin_one * scale_factor
        NAS_Smin_two = NAS_Smin_two * scale_factor
        test_min = C_index - k_allowance_NAS
    end if

    !==============================================================
    ! Step 3: Update Maximum Statistics (Smax)
    !==============================================================
    ! 3.1 Find location of the maximum value
    max_index = maxloc(tOnlineSample, dim=1)

    ! 3.2 Construct Rank Index Vector
    rankindex_s = delta_NAS
    rankindex_s(sampling_index) = 0.0d0
    rankindex_s(num_allnodes_addone) = 0.0d0

    global_idx = tsampling_index(max_index)

    if (global_idx == num_allnodes_addone) then
        rankindex_s(num_allnodes_addone) = lambda0_NAS
    else
        rankindex_s(global_idx) = lambda_NAS
    end if

    ! 3.3 Update Smax vectors (Using Internal State)
    NAS_Smax_one = NAS_Smax_one + rankindex_s
    NAS_Smax_two(1:num_allnodes) = NAS_Smax_two(1:num_allnodes) + gone_NAS
    NAS_Smax_two(num_allnodes_addone) = NAS_Smax_two(num_allnodes_addone) + gp_NAS

    ! 3.4 Compute Statistic C_index
    C_index = sum(((NAS_Smax_one - NAS_Smax_two)**2.0d0) / NAS_Smax_two)

    ! 3.5 Normalization / Reset
    if (C_index <= k_allowance_NAS) then
        NAS_Smax_one = 0.0d0
        NAS_Smax_two = 0.0d0
        test_max = 0.0d0
    else
        scale_factor = (C_index - k_allowance_NAS) / C_index
        NAS_Smax_one = NAS_Smax_one * scale_factor
        NAS_Smax_two = NAS_Smax_two * scale_factor
        test_max = C_index - k_allowance_NAS
    end if

    !==============================================================
    ! Step 4: Compute Final Charting Statistic
    !==============================================================
    charting_statistic = max(test_min, test_max)

    !==============================================================
    ! Step 5: Adaptive Node Selection
    !==============================================================
    ! Use the max of Min-stat and Max-stat for selection priority
    S_one(1:num_allnodes) = max(NAS_Smin_one(1:num_allnodes), NAS_Smax_one(1:num_allnodes))
    
    ! Add random perturbation for tie-breaking
    call RNNOR(temp_vector)
    temp_vector = S_one(1:num_allnodes) + eps_perturb * temp_vector
    
    ! Select top-k nodes indices
    call partial_quickselect(values=temp_vector, num_all=num_allnodes, &
                             top_k=num_samplingnodes, order_index=-1, &
                             iperm_sub=sampling_index)

end subroutine NAS

end module NAS_mod
