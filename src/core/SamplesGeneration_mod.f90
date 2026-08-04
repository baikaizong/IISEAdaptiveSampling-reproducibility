!==============================================================
!  Online sample generation for all simulation methods.
!
!  Module structure:
!    1. Load static matrices and nodes with Init_SamplesGeneration.
!    2. Configure the next sample with Set_SampleParams.
!    3. Generate the sample with GenerateOnlineSample.
!==============================================================
module SamplesGeneration_mod
    use GlobalSettings_mod
    use SparseMatrix_mod
    use RNMVN_INT
    use CHFAC_INT
    use RNNOF_INT
    use RNCHI_INT
    use RNEXP_INT
    use RNNOR_INT
    implicit none
    private 

    ! Public Interfaces
    public :: GenerateOnlineSample
    public :: Init_SamplesGeneration
    public :: Set_SampleParams
    public :: Clean_SamplesGeneration
    public :: AnomalyParams_type
    public :: generate_unique_random_points



    !==============================================================
    !  INTERNAL MODULE STATE VARIABLES
    !==============================================================
    
    ! 1. Static Heavy Data (Loaded once)
    real(dp), allocatable, save :: SG_Nodesset(:,:)
    type(Bbasis_sparse_type), save :: SG_B0_all
    type(Bbasis_sparse_type), save :: SG_B1_all
    type(Bbasis_sparse_type), save :: SG_CovMat
    
    ! 2. Dynamic control parameters
    integer, save :: SG_ICorOC
    integer, save :: SG_noise_type
    integer, save :: SG_back_index      ! Now acts as Background Type ID
    type(AnomalyParams_type), save :: SG_AnoParams

    ! 3. State Flags
    logical, save :: data_loaded_SG = .false.
    logical, save :: params_set_SG  = .false.

contains

    !==============================================================
    !  Init_SamplesGeneration
    !  Purpose: Load static geometry and matrices ONCE.
    !==============================================================
    subroutine Init_SamplesGeneration(Nodesset_in, B0_in, B1_in, Cov_in)
        implicit none
        real(dp), intent(in), optional :: Nodesset_in(:,:)
        type(Bbasis_sparse_type), intent(in), optional :: B0_in
        type(Bbasis_sparse_type), intent(in), optional :: B1_in
        type(Bbasis_sparse_type), intent(in), optional :: Cov_in

        ! Copy Nodesset
        if (present(Nodesset_in)) then
            if (allocated(SG_Nodesset)) deallocate(SG_Nodesset)
            allocate(SG_Nodesset(size(Nodesset_in,1), size(Nodesset_in,2)))
            SG_Nodesset = Nodesset_in
        end if

        ! Copy Sparse Matrices 
        if (present(B0_in))  call copy_Bbasis_sparse(B0_in, SG_B0_all)
        if (present(B1_in))  call copy_Bbasis_sparse(B1_in, SG_B1_all)
        if (present(Cov_in)) call copy_Bbasis_sparse(Cov_in, SG_CovMat)

        data_loaded_SG = .true.

    end subroutine Init_SamplesGeneration

    !==============================================================
    !  Set_SampleParams
    !  Purpose: Update control parameters for the next generation step.
    !==============================================================
    subroutine Set_SampleParams(ICorOC_in, noise_type_in, back_index_in, &
                                AnoParams_in, back_bandwidth_in)
        implicit none
        integer, intent(in), optional :: ICorOC_in, noise_type_in, back_index_in
        type(AnomalyParams_type), intent(in), optional :: AnoParams_in
        real(dp), intent(in), optional :: back_bandwidth_in

        if (present(ICorOC_in))     SG_ICorOC     = ICorOC_in
        if (present(noise_type_in)) SG_noise_type = noise_type_in
        if (present(back_index_in)) SG_back_index = back_index_in
        if (present(AnoParams_in))  SG_AnoParams  = AnoParams_in
        if (present(back_bandwidth_in)) kernel_bandwidth = back_bandwidth_in

        params_set_SG = .true.
    end subroutine Set_SampleParams

    !==============================================================
    !  Clean_SamplesGeneration
    !  Purpose: Free memory.
    !==============================================================
    subroutine Clean_SamplesGeneration()
        implicit none
        if (allocated(SG_Nodesset)) deallocate(SG_Nodesset)
        call free_Bbasis_sparse(SG_B0_all)
        call free_Bbasis_sparse(SG_B1_all)
        call free_Bbasis_sparse(SG_CovMat)
        data_loaded_SG = .false.
    end subroutine Clean_SamplesGeneration

    !==============================================================
    !  GenerateOnlineSample
    !  Main Entry Point: Simplified Interface.
    !==============================================================
    subroutine GenerateOnlineSample(sampling_index, OnlineSample, count_out)
        implicit none
        
        ! --- Input Arguments ---
        integer, intent(in) :: sampling_index(:)       
        
        ! --- Output Arguments ---
        real(dp), intent(out) :: OnlineSample(:) 
        ! Optional count of anomalous nodes in the generated sample
         integer, intent(out), optional :: count_out
        
        ! --- Local Variables ---
        integer :: n_samp
        real(dp), allocatable :: buffer(:)
        
        ! Safety Checks
        if (.not. data_loaded_SG) stop "Error: SamplesGeneration data not initialized."
        if (.not. params_set_SG)  stop "Error: SamplesGeneration params not set."

        n_samp = size(sampling_index)
        
        ! 1. Compute Noise Component
        call ComputeNoise(sampling_index, OnlineSample)

        ! 2. Add Background Component (if required)
        if (SG_back_index /= BACK_NONE) then
            allocate(buffer(n_samp))
            call ComputeBackground(sampling_index, buffer)
            OnlineSample = OnlineSample + buffer
            deallocate(buffer)
        end if

        ! 3. Add Anomaly Component (Only for OC state)
        if (SG_ICorOC == OC) then
            allocate(buffer(n_samp))
            call ComputeAnomaly(sampling_index, buffer,count_out)
            OnlineSample = OnlineSample + buffer
            deallocate(buffer)
        end if

    end subroutine GenerateOnlineSample

    !==============================================================
    !  ComputeNoise
    !  Generates noise using module-level CovMat and settings.
    !==============================================================
    subroutine ComputeNoise(idx, vec)
        implicit none
        integer, intent(in) :: idx(:)
        real(dp), intent(inout) :: vec(:)
        
        integer :: i, n
        real(dp) :: rvec(1)
        type(Bbasis_sparse_type) :: Cov_sub
        real(dp), allocatable :: SubCov(:,:), RSIG(:,:), tmp_noise(:)
        integer :: rank, nz, col, local_col
        integer, allocatable :: map_g2l(:)

        n = size(idx)

        select case (SG_noise_type)
        case (NOISE_GAUSSIAN)
            do i = 1, n
                vec(i) = RNNOF()
                vec(i) =vec(i) * sigma_noise
            end do

        case (NOISE_EXPONENTIAL)
            do i = 1, n
                call RNEXP(rvec)
                vec(i) = rvec(1) - 1.0d0
            end do

        case (NOISE_CHISQUARE)
            do i = 1, n
                call RNCHI(Df, rvec)
                vec(i) = (rvec(1) - Df) / sqrt(2.0d0 * Df)
            end do

        case (NOISE_CORR_GAUSSIAN)
            allocate(tmp_noise(n), SubCov(n,n), RSIG(n,n))
            allocate(map_g2l(num_allnodes))
            
            SubCov = 0.0d0
            map_g2l = 0
            
            do i = 1, n
                map_g2l(idx(i)) = i
            end do
            
            call build_sub_sparse(SG_CovMat, idx, Cov_sub)
            
            do i = 1, n
                nz = Cov_sub%nzcount(i)
                do col = 1, nz
                    local_col = map_g2l(Cov_sub%colind(i, col))
                    if (local_col > 0) SubCov(i, local_col) = Cov_sub%val(i, col)
                end do
            end do
            
            call CHFAC(SubCov, rank, RSIG)
            call DRNMVN(RSIG, tmp_noise)
            vec = tmp_noise
            
            call free_Bbasis_sparse(Cov_sub)
            deallocate(tmp_noise, SubCov, RSIG, map_g2l)
            
        case default
            stop "Unknown Noise Type"
        end select
    end subroutine ComputeNoise

!==============================================================
    !  ComputeBackground
    !  Generates background using Expanded Padding Strategy.
    !  Ensures strict variance consistency without periodic wrapping.
    !==============================================================
    subroutine ComputeBackground(idx, vec)
        implicit none
        integer, intent(in) :: idx(:)
        real(dp), intent(out) :: vec(:)
        
        ! Variables for BSpline
        type(Bbasis_sparse_type) :: B0_sub
        real(dp), allocatable :: theta0(:)
        real(dp), allocatable :: ExtendedRand(:) 
        real(dp), allocatable :: KMask(:,:)
        
        integer :: win_size, wx, wy, sx, sy
        integer :: pad_size, ext_nx, ext_ny, ext_idx
        integer :: center_x_in_ext, center_y_in_ext, cur_x, cur_y
        integer :: i, n, nid
        real(dp) :: accum, scale_factor, sq_sum, dist2, weight

        n = size(idx)

        select case (SG_back_index)
        
        case (BACK_BSPLINE)
            allocate(theta0(SG_B0_all%ncol))
            call RNNOR(theta0)
            theta0 = theta0 * sigma_ground
            
            call build_sub_sparse(SG_B0_all, idx, B0_sub)
            call spmv_Bbasis_sparse(B0_sub, theta0, vec)
            
            call free_Bbasis_sparse(B0_sub)
            deallocate(theta0)

        case (BACK_KERNEL_RANDOM)
            ! 1. Determine Kernel Size & Padding
            win_size = ceiling(3.0d0 * kernel_bandwidth)
            pad_size = win_size  ! Pad enough so kernel never goes out of bounds
            
            ! 2. Define Extended Grid Dimensions
            ext_nx = num_x + 2 * pad_size
            ext_ny = num_y + 2 * pad_size
            
            ! 3. Generate Noise on Expanded Grid
            allocate(ExtendedRand(ext_nx * ext_ny))
            call RNNOR(ExtendedRand) 
            ! ExtendedRand ~ N(0,1). We scale via kernel later.

            ! 4. Build Variance-Preserving Kernel (L2 Normalized)
            allocate(KMask(-win_size:win_size, -win_size:win_size))
            sq_sum = 0.0d0
            do wx = -win_size, win_size
                do wy = -win_size, win_size
                    dist2 = dble(wx*wx + wy*wy)
                    weight = exp(-dist2 / (2.0d0 * (kernel_bandwidth**2.0d0)))
                    KMask(wx, wy) = weight
                    sq_sum = sq_sum + weight**2.0d0
                end do
            end do
            
            ! Normalize: Sum(K^2) = 1
            if (sq_sum > 0.0d0) then
                KMask = KMask / sqrt(sq_sum)
            endif
            
            ! Target Scale Factor
            scale_factor = sigma_kernel

            ! 5. Convolution 
            do i = 1, n
                nid = idx(i)
                sx = (nid - 1)/num_x + 1
                sy = mod(nid - 1, num_x) + 1
                center_x_in_ext = sx + pad_size
                center_y_in_ext = sy + pad_size
                
                accum = 0.0d0
                
                do wx = -win_size, win_size
                    do wy = -win_size, win_size
                        ! Look up neighbor in Extended Grid
                        cur_x = center_x_in_ext + wx
                        cur_y = center_y_in_ext + wy
                        ext_idx = (cur_x - 1) * num_x + cur_y 
                        ext_idx = (cur_x - 1) * ext_ny + cur_y
                        
                        accum = accum + ExtendedRand(ext_idx) * KMask(wx, wy)
                    end do
                end do
                
                ! Apply final variance scaling
                vec(i) = accum * scale_factor
            end do
            
            deallocate(ExtendedRand)
            deallocate(KMask)

        case default
            vec = 0.0d0

        end select
    end subroutine ComputeBackground

!==============================================================
    !  ComputeAnomaly
    !  Generates anomaly using module-level Params.
    !==============================================================
    subroutine ComputeAnomaly(idx, vec, count_out)
    implicit none
    integer, intent(in) :: idx(:)
    real(dp), intent(out) :: vec(:)
    ! Optional count of anomalous nodes in the generated sample
    integer, intent(out), optional :: count_out
    
    integer :: i, n, nid, k, nz, j
    real(dp) :: cx, cy, nx_coord, ny_coord, dx, dy, dist2, rad2
    real(dp) :: val_base, cur_time
    real(dp) :: semi_a, semi_b     
    real(dp) :: inv_a2, inv_b2 
    real(dp) :: R_outer, d_shift, area_intersect, area_unit, shift_ratio, R_sq, term_d_2, dx_in, cx_shift
    type(Bbasis_sparse_type) :: B1_sub
    
    ! Variables for Convolution
    real(dp), allocatable :: KMask(:,:)
    integer :: win_size, wx, wy, sx, sy
    integer :: wrap_x, wrap_y, neighbor_k
    real(dp) :: accum, sq_sum, weight, scale_factor
    
    real(dp) :: l1_sum, min_dist
    real(dp), allocatable :: base_signal(:)
    integer :: best_node, xi, yi
    
    type(AnomalyParams_type) :: p 
    integer :: local_count

    n = size(idx)
    p = SG_AnoParams
    vec = 0.0d0 
    local_count = 0

    select case (p%type_id)
    
    case (TYPE_CIRCLE)
    ! Apply the selected anomaly model.
        if (.not. allocated(SG_Nodesset)) stop "Nodes required for TYPE_CIRCLE"
        cx = p%center_idx(1)
        cy = p%center_idx(2)
        rad2 = p%radius**2
        do i = 1, n
            nid = idx(i)
            dx = dist_period(SG_Nodesset(nid,1), cx)
            dy = dist_period(SG_Nodesset(nid,2), cy)
            if (dx*dx + dy*dy <= rad2) then
                vec(i) = p%value
                local_count = local_count + 1
            end if
        end do

    case (TYPE_ELLIPSE)
        if (.not. allocated(SG_Nodesset)) stop "Nodes required for TYPE_ELLIPSE"
        cx = p%center_idx(1)
        cy = p%center_idx(2)

        semi_b = sqrt(p%area / (pi_value * p%ellipticity))
        semi_a = semi_b * p%ellipticity

        inv_a2 = 1.0d0 / (semi_a**2)
        inv_b2 = 1.0d0 / (semi_b**2)

        do i = 1, n
            nid = idx(i)

            dx = dist_period(SG_Nodesset(nid,1), cx)
            dy = dist_period(SG_Nodesset(nid,2), cy)
    
            if ( (dx*dx * inv_a2) + (dy*dy * inv_b2) <= 1.0d0 ) then
                vec(i) = vec(i) + p%value 
                local_count = local_count + 1
            end if
        end do

    case (TYPE_CRESCENT)
        if (.not. allocated(SG_Nodesset)) stop "Nodes required for TYPE_CRESCENT"

        cx = p%center_idx(1)
        cy = p%center_idx(2)
        shift_ratio = p%ellipticity

        if (shift_ratio > 2.0d0) shift_ratio = 2.0d0 
        if (shift_ratio < 0.05d0) shift_ratio = 0.05d0
        
        if (shift_ratio > 1.999d0) then
            ! Almost a full circle (Ring shape logic if intended, but here seemingly simpler)
            R_outer = sqrt(p%area / pi_value)
            inv_a2  = 1.0d0 / (R_outer**2) 
    
            do i = 1, n
                nid = idx(i)
                dx = dist_period(SG_Nodesset(nid,1), cx)
                dy = dist_period(SG_Nodesset(nid,2), cy)

                if ( (dx*dx + dy*dy) * inv_a2 <= 1.0d0 ) then
                    vec(i) = vec(i) + p%value 
                    local_count = local_count + 1
                end if
            end do
    
        else
            term_d_2 = shift_ratio * 0.5d0
            if (term_d_2 > 1.0d0) term_d_2 = 1.0d0

            area_intersect = 2.0d0 * acos(term_d_2) - term_d_2 * sqrt(4.0d0 - shift_ratio*shift_ratio)
            area_unit = pi_value - area_intersect

            R_outer = sqrt(p%area / area_unit)
            R_sq = R_outer**2

            d_shift = shift_ratio * R_outer
            cx_shift = cx + d_shift
    
            do i = 1, n
                nid = idx(i)
                dx = dist_period(SG_Nodesset(nid,1), cx)
                dy = dist_period(SG_Nodesset(nid,2), cy)

                if ( (dx*dx + dy*dy) <= R_sq ) then
                    dx_in = dist_period(SG_Nodesset(nid,1), cx_shift)

                    if ( (dx_in*dx_in + dy*dy) > R_sq ) then
                        vec(i) = vec(i) + p%value 
                        local_count = local_count + 1
                    end if
                end if
            end do
        end if

    case (TYPE_SPATIAL)
        if (.not. allocated(SG_Nodesset)) stop "Nodes required for TYPE_SPATIAL"
        do i = 1, n
            nid = idx(i)
            nx_coord = SG_Nodesset(nid,1)
            ny_coord = SG_Nodesset(nid,2)
            rad2 = nx_coord*nx_coord + ny_coord*ny_coord
            if (rad2 >= 0.40d0 .and. rad2 <= 0.50d0) then
                vec(i) = p%theta * ((nx_coord - 2.0d0*ny_coord + 2.0d0)**2)
                local_count = local_count + 1
            end if
        end do

    case (TYPE_ST)
        if (.not. allocated(SG_Nodesset)) stop "Nodes required for TYPE_ST"
        cx = p%center_idx(1)
        cy = p%center_idx(2)
        rad2 = (p%radius + p%delta_r * dble(p%time_idx))**2
        cur_time = p%delta_t * dble(p%time_idx)
        
        do i = 1, n
            nid = idx(i)
            dx = dist_period(SG_Nodesset(nid,1), cx)
            dy = dist_period(SG_Nodesset(nid,2), cy)
            dist2 = dx*dx + dy*dy
            if (dist2 <= rad2) then
                vec(i) = (1.0d0 - dist2/rad2) * cur_time
                local_count = local_count + 1
            end if
        end do

    case (TYPE_BSPLINE)
        call build_sub_sparse(SG_B1_all, idx, B1_sub)
        do i = 1, n
            nz = B1_sub%nzcount(i)
            do k = 1, nz
                if (B1_sub%colind(i, k) == p%bspline_idx) then
                    vec(i) = p%value * B1_sub%val(i, k)
                    local_count = local_count + 1
                    exit
                end if
            end do
        end do
        call free_Bbasis_sparse(B1_sub)

case (TYPE_RANDOM_POINTS_CONV) 
    if (.not. allocated(SG_Nodesset)) stop "Nodes required for convolution"
    if (.not. allocated(p%points_coords)) stop "points_coords not allocated"


    win_size = ceiling(3.0d0 * kernel_bandwidth)
    allocate(KMask(-win_size:win_size, -win_size:win_size))
    
    l1_sum = 0.0d0
    do wx = -win_size, win_size
        do wy = -win_size, win_size
            dist2 = dble(wx*wx + wy*wy)
            weight = exp(-dist2 / (2.0d0 * (kernel_bandwidth**2.0d0)))
            KMask(wx, wy) = weight
            l1_sum = l1_sum + weight
        end do
    end do

    if (l1_sum > 0.0d0) KMask = KMask / l1_sum


    allocate(base_signal(num_x * num_y))
    base_signal = 0.0d0

    do k = 1, p%num_points
        min_dist = huge(1.0_dp)
        best_node = 1
        
        do j = 1, num_x * num_y
            dx = dist_period(SG_Nodesset(j, 1), p%points_coords(1, k))
            dy = dist_period(SG_Nodesset(j, 2), p%points_coords(2, k))
            dist2 = dx*dx + dy*dy
            
            if (dist2 < min_dist) then
                min_dist = dist2
                best_node = j
            end if
        end do
        base_signal(best_node) = base_signal(best_node) + p%value
    end do



    do i = 1, n
        nid = idx(i)
        xi = mod(nid - 1, num_x) + 1
        yi = (nid - 1) / num_x + 1
        
        accum = 0.0d0
        
        do wx = -win_size, win_size
            do wy = -win_size, win_size
                wrap_x = mod(xi + wx - 1, num_x)
                if (wrap_x < 0) wrap_x = wrap_x + num_x
                wrap_x = wrap_x + 1
                wrap_y = mod(yi + wy - 1, num_y)
                if (wrap_y < 0) wrap_y = wrap_y + num_y
                wrap_y = wrap_y + 1
                neighbor_k = (wrap_y - 1) * num_x + wrap_x
                accum = accum + base_signal(neighbor_k) * KMask(wx, wy)
            end do
        end do
        
        vec(i) = accum 
        if (abs(accum) > 1.0d-8) then
            local_count = local_count + 1
        end if
    end do

    deallocate(KMask)
    deallocate(base_signal)    

    case default
        print *, "Warning: Unknown Anomaly Type", p%type_id
    end select
    
    if (present(count_out)) then
        count_out = local_count
    end if
    
end subroutine ComputeAnomaly

    !==============================================================
    !  Auxiliary Routines
    !==============================================================
    
    subroutine BuildGaussianKernel(bw, Mask, w_size)
        real(dp), intent(in) :: bw
        real(dp), allocatable, intent(out) :: Mask(:,:)
        integer, intent(out) :: w_size
        integer :: wx, wy
        real(dp) :: dist2, w_sum
        
        w_size = ceiling(3.0d0 * bw)
        allocate(Mask(-w_size:w_size, -w_size:w_size))
        
        w_sum = 0.0d0
        do wx = -w_size, w_size
            do wy = -w_size, w_size
                dist2 = dble(wx*wx + wy*wy)
                Mask(wx, wy) = exp(-dist2 / (2.0d0 * bw**2))
                w_sum = w_sum + Mask(wx, wy)
            end do
        end do
        Mask = Mask / w_sum
    end subroutine BuildGaussianKernel

    elemental function dist_period(v1, v2) result(d)
        real(dp), intent(in) :: v1, v2
        real(dp) :: d
        d = abs(v1 - v2)
    end function dist_period
    
    subroutine generate_unique_random_points(num_points, coords)
    implicit none

    ! --- Arguments ---
    integer, intent(in) :: num_points
    real(dp), allocatable, intent(out) :: coords(:,:)

    ! --- Local Variables ---
    integer :: pts_generated, linear_idx, rand_x, rand_y
    real(dp) :: rand_val
    logical, allocatable :: used_mask(:)

    if (num_points > num_x * num_y) then
        stop "Error: Requested random points exceed grid size!"
    end if

    if (allocated(coords)) deallocate(coords)
    allocate(coords(2, num_points))

    allocate(used_mask(num_x * num_y))
    used_mask = .false.
    
    pts_generated = 0
    do while (pts_generated < num_points)
        call random_number(rand_val)

        linear_idx = int(rand_val * dble(num_x * num_y)) + 1
        
        if (linear_idx < 1) linear_idx = 1
        if (linear_idx > num_x * num_y) linear_idx = num_x * num_y

        if (.not. used_mask(linear_idx)) then
            used_mask(linear_idx) = .true.
            pts_generated = pts_generated + 1

            rand_x = mod(linear_idx - 1, num_x) + 1
            rand_y = (linear_idx - 1) / num_x + 1

            coords(1, pts_generated) = dble(rand_x)
            coords(2, pts_generated) = dble(rand_y)
        end if
    end do

    deallocate(used_mask)

end subroutine generate_unique_random_points

end module SamplesGeneration_mod
