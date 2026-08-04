module SASAM_mod
    use GlobalSettings_mod
    use RNNOF_INT
    use RNUND_INT
    use RNBIN_INT
    use SVRGP_INT
    use linear_operators
    use utils_mod
    use SparseMatrix_mod
    implicit none
    public 
    save
    
    !---------------------- SASAM Parameters ---------------------
    real(dp) :: Epah_SA = 0.02d0, miu_min_SA=1.0d0
    real(dp) :: theta_1_SA=0.1d0, theta_2_SA=0.7d0
    integer  :: Topr_SA=1    
    integer  :: num_kelnear_SA = 15

    
    real(dp) :: Kernel_threshold_SA=0.0d0
    
    logical  :: initialized_SASAM = .false. 
    logical  :: data_loaded_SASAM = .false.

    !---------------------- INTERNAL STATIC DATA --------------------
    ! Geometry Data (Read-only after init)
    real(dp), allocatable :: SA_Kernel_NeighborDis(:,:) ! (num_allnodes, num_kelnear_SA)
    integer,  allocatable :: SA_Kernel_NeighborId(:,:)  ! (num_allnodes, num_kelnear_SA)
    integer,  allocatable :: SA_NeighborId(:,:)         ! (num_allnodes, num_disnear)

    !---------------------- INTERNAL DYNAMIC STATE ------------------
    ! CUSUM Statistics (Updated iteratively)
    real(dp), allocatable :: SA_cusumdouble_statistic(:,:) ! (num_allnodes, 2)

    contains
    
!==================================================================
!  Subroutine: set_SASAMparams
!  Purpose: Initialize or update SASAM-specific scalar parameters
!==================================================================
subroutine set_SASAMparams( &
    Epah_SA_in, miu_min_SA_in, theta_1_SA_in, theta_2_SA_in, &
    Topr_SA_in, Kernel_threshold_SA_in, num_kelnear_SA_in)

    implicit none
    real(dp), intent(in), optional :: Epah_SA_in, miu_min_SA_in
    real(dp), intent(in), optional :: theta_1_SA_in, theta_2_SA_in, Kernel_threshold_SA_in
    integer, intent(in), optional :: Topr_SA_in, num_kelnear_SA_in

    if (present(Epah_SA_in))            Epah_SA             = Epah_SA_in
    if (present(miu_min_SA_in))         miu_min_SA          = miu_min_SA_in
    if (present(theta_1_SA_in))         theta_1_SA          = theta_1_SA_in
    if (present(theta_2_SA_in))         theta_2_SA          = theta_2_SA_in
    if (present(Topr_SA_in))            Topr_SA             = Topr_SA_in
    if (present(num_kelnear_SA_in))     num_kelnear_SA      = num_kelnear_SA_in
    if (present(Kernel_threshold_SA_in)) Kernel_threshold_SA = Kernel_threshold_SA_in

    initialized_SASAM = .true.

end subroutine set_SASAMparams

!==================================================================
!  Subroutine: Init_SASAM_Data
!  Purpose: 
!     1. Allocate and copy Geometry Data ONE TIME.
!     2. Allocate internal CUSUM state arrays.
!     3. Reset state to initial values.
!==================================================================
subroutine Init_SASAM_Data(Kernel_NeighborDis_in, Kernel_NeighborId_in, NeighborId_in)
    implicit none
    real(dp), intent(in) :: Kernel_NeighborDis_in(:,:)
    integer, intent(in)  :: Kernel_NeighborId_in(:,:)
    integer, intent(in)  :: NeighborId_in(:,:)

    ! Safety Check on Dimensions (Optional but recommended)
    if (size(Kernel_NeighborDis_in, 2) /= num_kelnear_SA) then
        print *, "Warning: Input Kernel Neighbor dimension mismatch with num_kelnear_SA"
    end if

    ! 1. Allocate & Copy Static Data
    if (allocated(SA_Kernel_NeighborDis)) deallocate(SA_Kernel_NeighborDis)
    allocate(SA_Kernel_NeighborDis(num_allnodes, num_kelnear_SA))
    SA_Kernel_NeighborDis = Kernel_NeighborDis_in

    if (allocated(SA_Kernel_NeighborId)) deallocate(SA_Kernel_NeighborId)
    allocate(SA_Kernel_NeighborId(num_allnodes, num_kelnear_SA))
    SA_Kernel_NeighborId = Kernel_NeighborId_in

    if (allocated(SA_NeighborId)) deallocate(SA_NeighborId)
    ! Assuming NeighborId second dimension matches global 'num_disnear' or input size
    allocate(SA_NeighborId(num_allnodes, size(NeighborId_in, 2)))
    SA_NeighborId = NeighborId_in

    ! 2. Allocate Dynamic State
    if (allocated(SA_cusumdouble_statistic)) deallocate(SA_cusumdouble_statistic)
    allocate(SA_cusumdouble_statistic(num_allnodes, 2))

    ! 3. Initialize State
    call Reset_SASAM_State()

    data_loaded_SASAM = .true.

end subroutine Init_SASAM_Data

!==================================================================
!  Subroutine: Reset_SASAM_State
!  Purpose: Reset CUSUM statistics to zero
!==================================================================
subroutine Reset_SASAM_State()
    implicit none
    if (allocated(SA_cusumdouble_statistic)) SA_cusumdouble_statistic = 0.0d0
end subroutine Reset_SASAM_State

!==================================================================
!  Subroutine: Clean_SASAM
!  Purpose: Free all internal memory
!==================================================================
subroutine Clean_SASAM()
    implicit none
    if (allocated(SA_Kernel_NeighborDis))    deallocate(SA_Kernel_NeighborDis)
    if (allocated(SA_Kernel_NeighborId))     deallocate(SA_Kernel_NeighborId)
    if (allocated(SA_NeighborId))            deallocate(SA_NeighborId)
    if (allocated(SA_cusumdouble_statistic)) deallocate(SA_cusumdouble_statistic)
    data_loaded_SASAM = .false.
end subroutine Clean_SASAM

!==================================================================
!  Subroutine: get_SASAMparams
!  Purpose: Print status
!==================================================================
subroutine get_SASAMparams(fid)
    implicit none
    integer, intent(in) :: fid  

    write(fid,'(A)') "==============================================="
    write(fid,'(A)') "          SASAM PARAMETER CONFIGURATION"
    write(fid,'(A)') "==============================================="
    write(fid,'(A,I8)') "num_samplingnodes   = ", num_samplingnodes
    write(fid,'(A,F12.6)') "Epah_SA             = ", Epah_SA
    write(fid,'(A,F12.6)') "miu_min_SA          = ", miu_min_SA
    write(fid,'(A,F12.6)') "theta_1_SA          = ", theta_1_SA
    write(fid,'(A,F12.6)') "theta_2_SA          = ", theta_2_SA
    write(fid,'(A,I8)')   "Topr_SA             = ", Topr_SA
    write(fid,'(A,F12.6)') "Kernel_threshold_SA = ", Kernel_threshold_SA
    write(fid,'(A,I8)')   "num_kelnear_SA      = ", num_kelnear_SA
    write(fid,'(A,L1)')   "data_loaded_SASAM   = ", data_loaded_SASAM
    write(fid,'(A)') "==============================================="
end subroutine get_SASAMparams

!======================================================================
!  Subroutine: SASAM
!  Geometry and CUSUM state are stored at module scope.
!======================================================================
subroutine SASAM(OnlineSample, limit, sampling_index, charting_statistic)
    implicit none

    !---------------- Inputs ----------------
    real(dp), intent(in) :: OnlineSample(num_samplingnodes)
    real(dp), intent(in) :: limit

    !---------------- InOut / Outputs ----------------
    ! sampling_index: Input = Current nodes, Output = Next nodes
    integer , intent(inout) :: sampling_index(num_samplingnodes)
    real(dp), intent(out)   :: charting_statistic

    !---------------- Locals ----------------
    integer :: i, j, k, Nodes
    integer :: one_vec(1), max_index
    integer :: D_nodes, W_nodes
    integer :: sampling_id(num_allnodes)
    real(dp) :: cusum_statistic(num_allnodes)
    real(dp) :: max_value
    real(dp) :: temp_vector(num_allnodes)

    if (.not. data_loaded_SASAM) then
        print *, "Error: SASAM data not initialized."
        stop
    end if

    !======================================================================
    ! Step 1: Update CUSUM statistics within kernel neighborhoods
    !         (Uses SA_Kernel_NeighborDis, SA_Kernel_NeighborId)
    !======================================================================
    do j = 1, num_samplingnodes
        do i = 1, num_kelnear_SA
            if (SA_Kernel_NeighborDis(sampling_index(j), i) <= Kernel_threshold_SA) exit
            
            ! Update double-directional CUSUM statistics
            k = SA_Kernel_NeighborId(sampling_index(j), i)
            
            SA_cusumdouble_statistic(k, 1) = SA_cusumdouble_statistic(k, 1) + &
                 SA_Kernel_NeighborDis(sampling_index(j),i) * (miu_min_SA * OnlineSample(j) - 0.5d0 * miu_min_SA**2.0d0)
            SA_cusumdouble_statistic(k, 2) = SA_cusumdouble_statistic(k, 2) + &
                 SA_Kernel_NeighborDis(sampling_index(j),i) * (-miu_min_SA * OnlineSample(j) - 0.5d0 * miu_min_SA**2.0d0)
        end do
    end do

    !======================================================================
    ! Step 2: Rectify negatives and compute global charting statistic
    !======================================================================
    do i = 1, num_allnodes
        SA_cusumdouble_statistic(i, 1) = max(0.0d0, SA_cusumdouble_statistic(i, 1))
        SA_cusumdouble_statistic(i, 2) = max(0.0d0, SA_cusumdouble_statistic(i, 2))
        cusum_statistic(i) = max(SA_cusumdouble_statistic(i, 1), SA_cusumdouble_statistic(i, 2))
    end do

    ! Compute top-r cumulative charting statistic
    temp_vector = cusum_statistic
    call partial_quickselect(values=temp_vector, num_all=num_allnodes, &
                             top_k=Topr_SA, order_index=-1, &
                             values_sub=temp_vector(1:Topr_SA))
                             
    charting_statistic = sum(temp_vector(1:Topr_SA))

    !======================================================================
    ! Step 3: Initialize node allocation
    !======================================================================
    sampling_index = 0
    sampling_id    = 0

    ! Determine the directed vs. random node counts
    max_value = charting_statistic - limit * theta_1_SA
    if (max_value < 0.0d0) max_value = 0.0d0

    ! Use limit here (ensure limit > 0 to avoid division by zero if applicable)
    if (limit > 1.0d-6) then
        D_nodes = int(max_value * dble(num_samplingnodes) * theta_2_SA / &
                      (limit * (1.0d0 - theta_1_SA)))
    else
        D_nodes = 0
    end if
    
    D_nodes = max(0, min(D_nodes, num_samplingnodes))
    W_nodes = num_samplingnodes - D_nodes

    !======================================================================
    ! Step 4: Directed allocation (focus near the maximum CUSUM region)
    !         (Uses SA_NeighborId)
    !======================================================================
    if (D_nodes > 0) then
        max_index = maxloc(cusum_statistic, 1)
        do i = 1, D_nodes
            ! Use internal NeighborId
            ! Warning: Ensure D_nodes <= num_disnear (size of SA_NeighborId's 2nd dim)
            ! Usually num_disnear > num_samplingnodes, so this is safe.
            sampling_id(SA_NeighborId(max_index, i)) = 1
            sampling_index(i) = SA_NeighborId(max_index, i)
        end do
    end if

    !======================================================================
    ! Step 5: Random fill for remaining W_nodes
    !======================================================================
    Nodes = D_nodes
    do while (Nodes < num_samplingnodes)
        call RNUND(num_allnodes, one_vec)
        if (sampling_id(one_vec(1)) == 1) cycle
        Nodes = Nodes + 1
        sampling_index(Nodes) = one_vec(1)
        sampling_id(one_vec(1)) = 1
    end do

end subroutine SASAM
                            
end module SASAM_mod
