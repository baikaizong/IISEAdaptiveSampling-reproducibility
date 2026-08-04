!==============================================================
!  Numerical, ranking, and sampling utility routines.
!  - Cubic B-spline: set degx=3, degy=3
!  - Open uniform knots
!  - Optional first derivatives w.r.t x and y
!  - Sampling: "midpoints" (default) or "endpoints"
!==============================================================
module utils_mod
    use GlobalSettings_mod
    use linear_operators
    use SparseMatrix_mod
    use SVRGP_INT
    implicit none
    public 
    contains

    subroutine ComputeCovMatrix(Nodesset, CovMat)
        implicit none
        real(dp), intent(in)  :: Nodesset(num_allnodes, 2)
        real(dp), intent(out) :: CovMat(num_allnodes, num_allnodes)
        integer :: i, j
        real(dp) :: dx, dy, r

        do i = 1, num_allnodes
            CovMat(i,i) = 1.0d0
            do j = i + 1, num_allnodes
                dx = Nodesset(i,1) - Nodesset(j,1)
                dy = Nodesset(i,2) - Nodesset(j,2)
                r  = abs(dx) + abs(dy)
                CovMat(i,j) = exp(-r / 0.02d0)
                CovMat(j,i) = CovMat(i,j)
            end do
        end do
    end subroutine ComputeCovMatrix
    
subroutine partial_quickselect(values, num_all, top_k, order_index, values_sub, iperm_sub)
    implicit none
    integer, intent(in) :: num_all, top_k, order_index
    real(dp), intent(in) :: values(num_all)
    real(dp), intent(out), optional :: values_sub(top_k)
    integer, intent(out), optional :: iperm_sub(top_k)

    real(dp):: values_subtem(top_k)
    integer :: iperm_subtem(top_k)
    real(dp) :: val_cmp_topk(num_all)
    real(dp) :: val_cmp_i, sgn
    integer :: i, j

    ! order_index = +1 selects the smallest top_k values.
    ! order_index = -1 selects the largest top_k values.
    sgn = real(order_index, dp)

    ! Initialize: store raw values but sort by val_cmp = sgn*value
    iperm_subtem  = [(i, i=1,top_k)]
    val_cmp_topk = sgn * values
    values_subtem = val_cmp_topk(1:top_k)
    ! Sort by val_cmp (ascending)
    call SVRGP(values_subtem, values_subtem, iperm_subtem)
    do i = top_k + 1, num_all
        if (val_cmp_topk(i) > values_subtem(top_k)) cycle
        ! From right to left find insertion point
        do j = top_k-1,1,-1
            if (val_cmp_topk(i) > values_subtem(j)) then
                ! Insert at position j+1
                values_subtem(j+2:top_k) = values_subtem(j+1:top_k-1)
                iperm_subtem(j+2:top_k)    = iperm_subtem(j+1:top_k-1)
                values_subtem(j+1) = val_cmp_topk(i)
                iperm_subtem(j+1)    = i
                exit
            end if
        end do

        ! Insert val_cmp_topk(i) at position 1 when it is the minimum.
        if (val_cmp_topk(i) <= values_subtem(1)) then
            values_subtem(2:top_k) = values_subtem(1:top_k-1)
            iperm_subtem(2:top_k)    = iperm_subtem(1:top_k-1)
            values_subtem(1) = val_cmp_topk(i)
            iperm_subtem(1)    = i
        end if
    end do
    if (present(values_sub)) then
        values_sub = sgn *values_subtem
    end if

    if (present(iperm_sub)) then
        iperm_sub=iperm_subtem
    end if    
end subroutine partial_quickselect


subroutine NearestNeighbors(Nodesset, k_neigh, board_index, AllNeighborDis, AllNeighborId)
    use omp_lib
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none
    
    ! --- Arguments ---
    integer, intent(in) :: k_neigh, board_index
    real(dp), intent(in)  :: Nodesset(num_allnodes,2)
    real(dp), intent(out) :: AllNeighborDis(num_allnodes, k_neigh)
    integer, intent(out)  :: AllNeighborId(num_allnodes, k_neigh)

    ! --- Local Variables ---
    integer :: i, j
    real(dp) :: dx, dy
    real(dp), allocatable :: DistMat(:) 

    ! --- Progress Bar Variables ---
    integer :: prog_cnt, my_cnt
    integer :: step_size, bar_len, filled_len, pct
    character(len=1) :: cr
    character(len=20) :: bar_str

    !-------------------- Basic check --------------------
    if (k_neigh <= 0 .or. k_neigh > num_allnodes) then
        print *, "Error: invalid k_neigh"
        return
    end if

    ! --- Init Progress Params ---
    prog_cnt = 0
    bar_len = 20
    cr = char(13)
    step_size = max(10, num_allnodes / 100) 

    print *, ">> Calculating Nearest Neighbors (OpenMP)..."

    !-------------------- Parallel Computation --------------------
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, dx, dy, DistMat, my_cnt, pct, filled_len, bar_str)
    
        allocate(DistMat(num_allnodes))

        !$OMP DO SCHEDULE(DYNAMIC)
        do i = 1, num_allnodes
            
            DistMat = 0.0d0
            do j = 1, num_allnodes    
                dx = abs(Nodesset(i,1) - Nodesset(j,1))
                dy = abs(Nodesset(i,2) - Nodesset(j,2))
                if (board_index == 1) then
                    dx = min(dx, (1.0d0 - dx))
                    dy = min(dy, (1.0d0 - dy))
                end if
                
                DistMat(j) = dx*dx + dy*dy
            end do
            
            call partial_quickselect(DistMat, num_allnodes, k_neigh, 1, &
                                     AllNeighborDis(i,:), AllNeighborId(i,:))
            

            !$OMP ATOMIC CAPTURE
            prog_cnt = prog_cnt + 1
            my_cnt = prog_cnt
            !$OMP END ATOMIC

            if (mod(my_cnt, step_size) == 0 .or. my_cnt == num_allnodes) then
                

                !$OMP CRITICAL(print_prog)
                    pct = int(real(my_cnt) / real(num_allnodes) * 100.0)
                    
                    filled_len = int(real(pct)/100.0 * real(bar_len))
                    if (filled_len < 0) filled_len = 0
                    if (filled_len > bar_len) filled_len = bar_len
                    bar_str = repeat('#', filled_len) // repeat('.', bar_len - filled_len)
                    
                    write(output_unit, '(A, A, I3, A, A, A)', advance='no') &
                        cr, " Progress: ", pct, "% [", trim(bar_str), "]"
                    flush(output_unit)
                !$OMP END CRITICAL(print_prog)
                
            end if

        end do
        !$OMP END DO
        
        deallocate(DistMat)
        
    !$OMP END PARALLEL

    write(output_unit, *) 

end subroutine NearestNeighbors
    !====================================================================
    !  Subroutine: EpaKernelTopK
    !  Purpose   : Compute pairwise Epanechnikov kernel values once,
    !              then select top-k for each node (descending order).
    !====================================================================
    subroutine EpaKernelTopK(Nodesset, Epah, k_top, AllKernelDis, AllKernelId)
        implicit none
        real(dp), intent(in)  :: Nodesset(num_allnodes,2), Epah
        integer,  intent(in)  :: k_top
        real(dp), intent(out) :: AllKernelDis(num_allnodes,k_top)
        integer,  intent(out) :: AllKernelId(num_allnodes,k_top)

        real(dp), allocatable :: Kmat(:)
        real(dp) :: dx, dy, r2, h2
        integer  :: i, j

        if (k_top <= 0 .or. k_top > num_allnodes) then
            print *, "Error: invalid k_top.";  return
        end if

        h2 = Epah**2.0d0
        allocate(Kmat(num_allnodes))

        !-------------------- precompute symmetric kernel matrix --------------------
        do i = 1, num_allnodes
            Kmat=0.0_dp
            do j = 1, num_allnodes
                dx = abs(Nodesset(i,1)-Nodesset(j,1))
                dy = abs(Nodesset(i,2)-Nodesset(j,2))
                r2 = (dx*dx + dy*dy)/h2
                Kmat(j) =1.0d0 - r2
                Kmat(j)=MAX(Kmat(j),0.0d0)
            end do
            call partial_quickselect( Kmat, num_allnodes,k_top, -1, AllKernelDis(i,:),AllKernelId(i,:)) 
        end do
        deallocate(Kmat)
    end subroutine EpaKernelTopK
                     
    !==============================================================
    ! Efraimidis-Spirakis weighted sampling without replacement
    ! Efficient O(n log m) implementation using random priority keys
    !
    ! Input:
    !   weights(n)     - non-negative sampling weights (probability density)
    !   n              - total number of candidates
    !   m              - number of samples to draw (m <= n)
    ! Output:
    !   sample_index(m)- selected indices (1-based)
    !
    ! Notes:
    !   - Works without replacement.
    !   - Handles zero weights safely (ignored automatically).
    !==============================================================                           
    subroutine weighted_sampling(weights, n, m, sample_index)

        implicit none
        integer, intent(in) :: n, m
        real(dp), intent(in) :: weights(n)
        integer, intent(out) :: sample_index(m)

        real(dp), allocatable :: key(:)
        real(dp) :: u, logu
        integer :: i, j, k
        integer, allocatable :: idx(:)
        integer :: temp_idx
        real(dp) :: temp_key

        !--------------------------------------------------------------
        ! Step 1. Allocate arrays
        !--------------------------------------------------------------
        allocate(key(n))
        allocate(idx(n))

        !--------------------------------------------------------------
        ! Step 2. Generate priority keys based on weights
        !   key_i = log(U) / weight_i  (log form for numerical stability)
        !--------------------------------------------------------------
        do i = 1, n
            if (weights(i) <= 0.0d0) then
                key(i) = -1.0d300   ! effectively ignore zero-weight elements
            else
                call random_number(u)
                if (u <= 1.0d-300) u = 1.0d-300
                logu = log(u)
                key(i) = logu / weights(i)
            end if
            idx(i) = i
        end do

        !--------------------------------------------------------------
        ! Step 3. Partial sort (descending order of key)
        !   We only need the top m largest keys.
        !   Simple quickselect or partial bubble since n not huge.
        !--------------------------------------------------------------
        do i = 1, m
            do j = i + 1, n
                if (key(j) > key(i)) then
                    temp_key = key(i); key(i) = key(j); key(j) = temp_key
                    temp_idx = idx(i); idx(i) = idx(j); idx(j) = temp_idx
                end if
            end do
        end do

        !--------------------------------------------------------------
        ! Step 4. Output the selected indices
        !--------------------------------------------------------------
        sample_index(:) = idx(1:m)

        !--------------------------------------------------------------
        ! Step 5. Cleanup
        !--------------------------------------------------------------
        deallocate(key)
        deallocate(idx)
    end subroutine weighted_sampling

    !======================================================================
    !  Subroutine: compute_sample_cov
    !  Purpose:
    !     Compute the sample covariance matrix of a data matrix X (n x p).
    !
    !  Input:
    !     X      - Data matrix (n x p), each row is one sample
    !     n      - Number of samples
    !     p      - Number of variables
    !
    !  Output:
    !     CovMat - Sample covariance matrix (p x p)
    !     MeanVec (optional) - Sample mean vector (p x 1)
    !
    !
    !======================================================================
    subroutine Covestimate(X, n, p, CovMat, MeanVec)
        implicit none
        integer, intent(in) :: n, p
        real(dp), intent(in) :: X(n, p)
        real(dp), intent(out) :: CovMat(p, p)
        real(dp), intent(out), optional :: MeanVec(p)

        !---------------- Local variables ----------------
        real(dp) :: Mean_local(p)
        real(dp) :: X_centered(n, p)
        integer :: i

        !==================================================
        ! Step 1. Compute sample mean for each variable
        !==================================================
        Mean_local = sum(X, dim=1) / dble(n)
        if (present(MeanVec)) MeanVec = Mean_local

        !==================================================
        ! Step 2. Center the data
        !==================================================
        do i = 1, p
            X_centered(:, i) = X(:, i) - Mean_local(i)
        end do

        !==================================================
        ! Step 3. Compute covariance matrix
        !==================================================
        CovMat =1.0 / dble(n- 1)*(X_centered .tx. X_centered)
    end subroutine Covestimate

    !====================================================================
    ! Subroutine: process_sparse_neighbors
    ! Purpose: Extracts top-K nearest neighbors from a Sparse Covariance Matrix.
    ! A top-k buffer avoids sorting the full input array.
    !            This is much faster when rows are dense.
    !====================================================================
    subroutine process_sparse_neighbors(S, OutId, OutVal)
        implicit none
        
        ! Inputs
        type(Bbasis_sparse_type), intent(in) :: S
        
        ! Outputs (Pre-allocated: num_allnodes x num_varnear_CDS)
        integer, intent(out) :: OutId(:,:)   
        real(dp), intent(out) :: OutVal(:,:) 
        
        ! Local variables
        integer :: i, j, k, col, row_nz
        integer :: max_k
        real(dp) :: val_curr
        
        ! Small temporary buffers for top-k maintenance
        ! These replace the huge temp arrays from before.
        integer :: top_ids(num_varnear_CDS)
        real(dp) :: top_vals(num_varnear_CDS)
        
        max_k = num_varnear_CDS
        
        ! Loop over all nodes
        do i = 1, num_allnodes
            row_nz = S%nzcount(i)
            
            ! Reset Top-K buffer
            top_ids = 0
            top_vals = 0.0d0
            
            ! Iterate through the sparse row
            do j = 1, row_nz
                col = S%colind(i, j)
                val_curr = S%val(i, j)
                
                ! Logic: Exclude self (diagonal)
                if (col /= i) then
                    ! Attempt to insert into the Top-K buffer
                    call update_top_k(max_k, top_vals, top_ids, val_curr, col)
                end if
            end do
            
            ! Copy result to output
            OutId(i, :)  = top_ids(:)
            OutVal(i, :) = top_vals(:)
        end do

    end subroutine process_sparse_neighbors

    !====================================================================
    ! Subroutine: update_top_k
    ! Purpose: Inserts a value into a sorted array ONLY IF it belongs 
    !          in the top K (based on absolute magnitude).
    !          The array is kept sorted descending by ABS(value).
    ! Efficiency: Avoids sorting the entire row. Complexity O(K).
    !====================================================================
    subroutine update_top_k(k, vals, ids, new_val, new_id)
        implicit none
        integer, intent(in) :: k
        real(dp), intent(inout) :: vals(k)
        integer, intent(inout) :: ids(k)
        real(dp), intent(in) :: new_val
        integer, intent(in) :: new_id
        
        integer :: i, j
        real(dp) :: abs_new, abs_curr
        
        abs_new = abs(new_val)
        
        ! 1. Quick Check: Is the new value smaller than the smallest value in the buffer?
        !    The buffer is sorted descending, so the last element (k) is the smallest.
        !    If new value is smaller, we don't need to do anything.
        if (abs_new <= abs(vals(k))) return
        
        ! 2. Find insertion position
        !    We look from the beginning (largest) to find where this fits.
        do i = 1, k
            abs_curr = abs(vals(i))
            
            if (abs_new > abs_curr) then
                ! Found the spot at index 'i'.
                ! Shift elements from i to k-1 downwards to make room.
                ! (The last element vals(k) falls off the cliff)
                do j = k, i + 1, -1
                    vals(j) = vals(j-1)
                    ids(j)  = ids(j-1)
                end do
                
                ! Insert new element
                vals(i) = new_val
                ids(i)  = new_id
                return
            end if
        end do
        
    end subroutine update_top_k
    

subroutine Compute_Kernel_Covariance(Cov_dense)
    use GlobalSettings_mod
    implicit none

    !-------------------- Outputs --------------------
    real(dp), allocatable, intent(out) :: Cov_dense(:,:)
    
    !-------------------- Local Variables --------------------
    integer :: num_nodes
    integer :: i, j, xi, yi, xj, yj, dx, dy, wx, wy
    integer :: win_size, cov_win_size
    real(dp) :: bw, l1_sum, weight, dist2, cov_val
    
    real(dp), allocatable :: KernelMask(:,:)
    real(dp), allocatable :: AutoCovKernel(:,:)
    

    num_nodes = num_x * num_y
    bw = kernel_bandwidth
    win_size = ceiling(3.0d0 * bw)  
    
    if (allocated(Cov_dense)) deallocate(Cov_dense)
    allocate(Cov_dense(num_nodes, num_nodes))
    Cov_dense = 0.0d0
    

    allocate(KernelMask(-win_size:win_size, -win_size:win_size))
    l1_sum = 0.0d0
    
    do wx = -win_size, win_size
        do wy = -win_size, win_size
            dist2 = dble(wx*wx + wy*wy)
            weight = exp(-dist2 / (2.0d0 * (bw**2.0d0)))
            KernelMask(wx, wy) = weight
            l1_sum = l1_sum + weight
        end do
    end do
    

    if (l1_sum > 0.0d0) then
        KernelMask = (KernelMask / l1_sum) * sigma_kernel
    else
        KernelMask = 0.0d0
    endif


    cov_win_size = 2 * win_size
    allocate(AutoCovKernel(-cov_win_size:cov_win_size, -cov_win_size:cov_win_size))
    AutoCovKernel = 0.0d0
    
    do dx = -cov_win_size, cov_win_size
        do dy = -cov_win_size, cov_win_size
            cov_val = 0.0d0
            do wx = max(-win_size, dx - win_size), min(win_size, dx + win_size)
                do wy = max(-win_size, dy - win_size), min(win_size, dy + win_size)
                    cov_val = cov_val + KernelMask(wx, wy) * KernelMask(wx - dx, wy - dy)
                end do
            end do
            AutoCovKernel(dx, dy) = cov_val
        end do
    end do
    
    do i = 1, num_nodes
        xi = mod(i - 1, num_x) + 1
        yi = (i - 1) / num_x + 1
        do dy = -cov_win_size, cov_win_size
            yj = yi + dy
            if (yj < 1 .or. yj > num_y) cycle 
            
            do dx = -cov_win_size, cov_win_size
                xj = xi + dx
                if (xj < 1 .or. xj > num_x) cycle 
                j = (yj - 1) * num_x + xj
                Cov_dense(i, j) = AutoCovKernel(dx, dy)
            end do
        end do
    end do

    deallocate(KernelMask)
    deallocate(AutoCovKernel)

end subroutine Compute_Kernel_Covariance

subroutine compute_conv_stddev( output_std)
    implicit none

    ! --- Arguments ---
    real(dp), intent(out) :: output_std        ! Standard deviation after convolution

    ! --- Local Variables ---
    integer  :: win_size, wx, wy
    real(dp) :: dist2, weight, l1_sum
    real(dp), allocatable :: KMask(:,:)

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

    output_std = sqrt((sigma_kernel**2.0) * sqrt(sum(KMask**2))+sigma_noise**2.0)
    deallocate(KMask)

end subroutine compute_conv_stddev

end module utils_mod
