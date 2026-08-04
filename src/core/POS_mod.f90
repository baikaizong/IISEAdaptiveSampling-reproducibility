module POS_mod
    use GlobalSettings_mod
    use RNNOF_INT
    use RNUND_INT
    use RNBIN_INT
    use SVRGP_INT
    use utils_mod
    implicit none
    public 
    save 
    !---------------------- POS parameters --------------------- 
    integer :: omega_POS=5,statisticwindowsize_POS=100,dim_POS=2
    real(dp) :: c1_POS=1.5d0, c2_POS=0.1d0
    real(dp) :: he_POS=0.2255, ht_POS=0.0438, zetap_POS=5.085E-8
    real(dp) :: Kernel_threshold_POS = 0.0d0
    integer ::  num_tkelnear_POS = 500,num_ekelnear_POS = 500
    logical :: initialized_POS = .false.
    
    contains
    
!==================================================================
!  Subroutine: set_POSparams
!  Purpose:
!     Initialize or update POS-specific parameters
!==================================================================
subroutine set_POSparams( &
    omega_POS_in, statisticwindowsize_POS_in, dim_POS_in, &
    c1_POS_in, c2_POS_in, he_POS_in, ht_POS_in, zetap_POS_in,Kernel_threshold_POS_in,num_tkelnear_POS_in,num_ekelnear_POS_in)
    implicit none

    !------------------ Optional Inputs ------------------
    integer, intent(in), optional :: omega_POS_in, statisticwindowsize_POS_in, dim_POS_in,num_tkelnear_POS_in,num_ekelnear_POS_in
    real(dp), intent(in), optional :: c1_POS_in, c2_POS_in, he_POS_in, ht_POS_in, zetap_POS_in,Kernel_threshold_POS_in

    !------------------ Parameter Assignments ------------------
    if (present(omega_POS_in))              omega_POS = omega_POS_in
    if (present(statisticwindowsize_POS_in)) statisticwindowsize_POS = statisticwindowsize_POS_in
    if (present(dim_POS_in))                dim_POS = dim_POS_in

    if (present(c1_POS_in)) c1_POS = c1_POS_in
    if (present(c2_POS_in)) c2_POS = c2_POS_in
    if (present(he_POS_in)) he_POS = he_POS_in
    if (present(ht_POS_in)) ht_POS = ht_POS_in
    if (present(zetap_POS_in)) zetap_POS = zetap_POS_in
    if (present(Kernel_threshold_POS_in)) Kernel_threshold_POS = Kernel_threshold_POS_in
    if (present(num_tkelnear_POS_in)) num_tkelnear_POS = num_tkelnear_POS_in
    if (present(num_ekelnear_POS_in)) num_ekelnear_POS = num_ekelnear_POS_in
    he_POS = c1_POS * sqrt(1.0d0/12.0d0) * (num_samplingnodes**(-1.0d0/(4.0d0 + dim_POS)))
    ht_POS = he_POS**(2.0d0 + c2_POS)
    zetap_POS = 0.01d0 / dble(num_allnodes) * (he_POS**2.0d0)

    !------------------ Status ------------------
    initialized_POS = .true.
end subroutine set_POSparams

!==================================================================
!  Subroutine: get_POSparams
!  Purpose:
!     Print all POS parameters to a given file handle
!==================================================================
subroutine get_POSparams(fid)
    implicit none
    integer, intent(in) :: fid  ! File handle for output

    write(fid,'(A)') "==============================================="
    write(fid,'(A)') "           POS PARAMETER CONFIGURATION"
    write(fid,'(A)') "==============================================="
    write(fid,'(A,I8)') "num_samplingnodes   = ", num_samplingnodes 
    write(fid,'(A,I8)')   "omega_POS           = ", omega_POS
    write(fid,'(A,I8)')   "statisticwindowsize_POS = ", statisticwindowsize_POS
    write(fid,'(A,I8)')   "dim_POS             = ", dim_POS

    write(fid,'(A,F12.6)') "c1_POS              = ", c1_POS
    write(fid,'(A,F12.6)') "c2_POS              = ", c2_POS
    write(fid,'(A,F12.6)') "he_POS              = ", he_POS
    write(fid,'(A,F12.6)') "ht_POS              = ", ht_POS
    write(fid,'(A,F12.6)') "zetap_POS           = ", zetap_POS

    write(fid,'(A,L1)')   "initialized_POS     = ", initialized_POS

    write(fid,'(A)') "==============================================="
    write(fid,'(A)') "End of POS PARAMETER listing"
    write(fid,'(A)') "==============================================="
end subroutine get_POSparams

!======================================================================
!  Subroutine: POS
!  Purpose:
!     Perform Probability-Oriented Sampling (POS) procedure:
!       (1) Compute spatial-temporal correlations based on kernel neighbors
!       (2) Update temporal windowed mean (miu_omega)
!       (3) Recalculate sampling density
!       (4) Perform weighted resampling based on updated density
!
!  Note:
!     Module state persists between successive online observations.
!======================================================================
subroutine POS(OnlineSample, &
                 Kernelhe_NeighborDis, Kernelhe_NeighborId, &
                 Kernelht_NeighborDis, Kernelht_NeighborId, &
                 miu_omega, statistic_window, &
                 mean_charting, Sigma_charting, &
                 sampling_density,sampling_index,charting_statistic, test_statistic)
    implicit none

    !-----------------------------
    ! Input variables
    !-----------------------------
    ! Note: OnlineSample contains data in 1:num_samplingnodes. 
    ! Declared size num_allnodes to match interface, but treated as compact.
    real(dp), intent(in) :: OnlineSample(num_allnodes)
    
    ! Kernel Neighbor Matrices (Read-only, large)
    real(dp), intent(in) :: Kernelhe_NeighborDis(num_allnodes, num_ekelnear_POS)
    real(dp), intent(in) :: Kernelht_NeighborDis(num_allnodes, num_tkelnear_POS)
    integer, intent(in) :: Kernelhe_NeighborId(num_allnodes, num_ekelnear_POS)
    integer, intent(in) :: Kernelht_NeighborId(num_allnodes, num_tkelnear_POS)

    !-----------------------------
    ! In/Out variables
    !-----------------------------
    real(dp), intent(inout) :: miu_omega(num_allnodes, omega_POS)
     real(dp), intent(inout) :: sampling_density(num_allnodes)
    integer, intent(inout) :: sampling_index(num_samplingnodes)

    !-----------------------------
    ! Optional / Output variables
    !-----------------------------
    real(dp), intent(inout), optional :: statistic_window(statisticwindowsize_POS)
    real(dp), intent(in),    optional :: mean_charting, Sigma_charting
    real(dp), intent(out),   optional :: charting_statistic
    real(dp), intent(out),   optional :: test_statistic

    !-----------------------------
    ! Local variables (Automatic Arrays)
    !-----------------------------
    integer :: i, j, k, near_id, curr_node
    real(dp) :: temp_value, temp_stat_accum
    real(dp) :: norm_factor
    
    ! Automatic arrays for better performance on stack/heap
    real(dp) :: inv_sqrt_dens(num_allnodes)
    real(dp) :: kernel_sum(num_allnodes)
    real(dp) :: sampling_observation(num_allnodes)
    real(dp) :: row_sum_miu(num_allnodes)

    
    ! Precompute inverse sqrt to replace division/sqrt in loops
    inv_sqrt_dens = 1.0d0 / sqrt(sampling_density)

    !==============================================================
    ! Step 1: Map compact observations to sparse global vector
    !==============================================================
    sampling_observation = 0.0d0
    do j = 1, num_samplingnodes
        sampling_observation(sampling_index(j)) = OnlineSample(j)
    end do

    !==============================================================
    ! Step 2: Compute temporal correlation statistic
    ! Reuse the precomputed inverse square-root densities and scalars.
    !==============================================================
    temp_stat_accum = 0.0d0
    
    do i = 1, num_samplingnodes
        curr_node = sampling_index(i)
        temp_value = 0.0d0

        ! Inner loop: Check temporal neighbors
        do j = 2, num_tkelnear_POS
            ! Check threshold first to avoid unnecessary memory access
            if (Kernelht_NeighborDis(curr_node, j) <= Kernel_threshold_POS) exit
            
            near_id = Kernelht_NeighborId(curr_node, j)
            
            ! Skip if neighbor was not sampled (observation is 0)
            if (sampling_observation(near_id) == 0.0d0) cycle

            temp_value = temp_value + Kernelht_NeighborDis(curr_node, j) * &
                         sampling_observation(near_id) * inv_sqrt_dens(near_id)
        end do
        
        temp_stat_accum = temp_stat_accum + &
                          temp_value * sampling_observation(curr_node) * inv_sqrt_dens(curr_node)
    end do

    temp_stat_accum = temp_stat_accum / dble(num_samplingnodes - 1)
    if (present(test_statistic)) test_statistic = temp_stat_accum

    !==============================================================
    ! Step 3: Compute charting statistic (Optional CUSUM-like)
    !==============================================================
    if (present(charting_statistic) .and. present(mean_charting) .and. &
        present(Sigma_charting) .and. present(statistic_window)) then

        ! Shift window efficiently
        statistic_window(2:statisticwindowsize_POS) = statistic_window(1:statisticwindowsize_POS - 1)
        statistic_window(1) = (temp_stat_accum - mean_charting) / Sigma_charting

        ! Find max of accumulated window statistics (Prefix Sum logic)
        charting_statistic = 0.0d0
        temp_value = 0.0d0
        
        do i = 1, statisticwindowsize_POS
            temp_value = temp_value + statistic_window(i)
            ! Optimization: compare squared values to avoid repeated sqrt inside loop
            ! assuming statistic calculation intends positive magnitude
            if (temp_value > 0.0d0) then 
                 temp_stat_accum = (temp_value**2.0d0) / dble(i)
                 if (charting_statistic < temp_stat_accum) charting_statistic = temp_stat_accum
            endif
        end do
        charting_statistic = sqrt(charting_statistic)
    end if

    !==============================================================
    ! Step 4: Update the spatial-temporal mean field (miu_omega)
    !==============================================================
    ! Shift history columns to the right
    miu_omega(:, 2:omega_POS) = miu_omega(:, 1:omega_POS - 1)
    
    ! Reset current column and accumulator
    miu_omega(:, 1) = 0.0d0
    kernel_sum = 0.0d0

    ! Scatter-Add: Distribute observed values to spatial neighbors
    do i = 1, num_samplingnodes
        curr_node = sampling_index(i)
        temp_value = OnlineSample(i) ! scalar cache

        do j = 1, num_ekelnear_POS
            if (Kernelhe_NeighborDis(curr_node, j) <= Kernel_threshold_POS) exit            
            
            near_id = Kernelhe_NeighborId(curr_node, j)
            
            miu_omega(near_id, 1) = miu_omega(near_id, 1) + &
                                    Kernelhe_NeighborDis(curr_node, j) * temp_value
            
            kernel_sum(near_id) = kernel_sum(near_id) + Kernelhe_NeighborDis(curr_node, j)
        end do
    end do

    ! Normalize the new column by kernel weights
    where (kernel_sum > 1.0d-12) ! Use a small epsilon instead of 0.0
        miu_omega(:, 1) = miu_omega(:, 1) / kernel_sum
    elsewhere
        miu_omega(:, 1) = 0.0d0
    end where

    !==============================================================
    ! Step 5: Update Sampling Density for NEXT Step
    !==============================================================
    ! Re-calculate based on updated miu_omega
    row_sum_miu = sum(miu_omega, dim=2)
    sampling_density = (row_sum_miu / dble(omega_POS))**2.0d0
    sampling_density = max(sampling_density, zetap_POS)
    
    ! Normalize density
    norm_factor = 1.0d0 / sum(sampling_density)
    sampling_density = sampling_density * norm_factor

    !==============================================================
    ! Step 6: Select new nodes (Weighted Sampling)
    !==============================================================
    call weighted_sampling(sampling_density, num_allnodes, num_samplingnodes, sampling_index)

end subroutine POS

                 
                 
end module POS_mod
