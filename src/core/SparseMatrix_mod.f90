module SparseMatrix_mod
    use GlobalSettings_mod
    implicit none    
    ! Define the sparse matrix structure (ELLPACK format)
    type :: Bbasis_sparse_type
        integer :: nrow      
        integer :: ncol      
        integer :: max_nz   
        integer, allocatable :: nzcount(:)   ! Dimension: (nrow)
        integer, allocatable :: colind(:,:)  ! Dimension: (nrow, max_nz)
        real(dp), allocatable :: val(:,:)    ! Dimension: (nrow, max_nz)
    end type Bbasis_sparse_type

    contains
    
!==================================================================
!  Subroutine: copy_Bbasis_sparse
!  Purpose: Deep copy of Bbasis_sparse_type structure.
!==================================================================
subroutine copy_Bbasis_sparse(src, dest)
    implicit none
    type(Bbasis_sparse_type), intent(in)  :: src
    type(Bbasis_sparse_type), intent(out) :: dest
    dest%nrow   = src%nrow
    dest%ncol   = src%ncol
    dest%max_nz = src%max_nz
    if (allocated(src%nzcount)) then
        ! Allocate with the size of source (usually nrow)
        allocate(dest%nzcount(size(src%nzcount)))
        dest%nzcount = src%nzcount
    end if
    if (allocated(src%colind)) then
        ! Allocate with exact dimensions of source
        allocate(dest%colind(size(src%colind, 1), size(src%colind, 2)))
        dest%colind = src%colind
    end if
    if (allocated(src%val)) then
        allocate(dest%val(size(src%val, 1), size(src%val, 2)))
        dest%val = src%val
    end if

end subroutine copy_Bbasis_sparse

    ! ==================================================================
    ! Subroutine: build_Bbasis_sparse
    ! Description: Converts a dense matrix into the sparse structure.
    ! ==================================================================
    subroutine build_Bbasis_sparse(Bbasis_all, Bbasiss)
        implicit none
        real(dp), intent(in)  :: Bbasis_all(:,:)   
        type(Bbasis_sparse_type), intent(out) :: Bbasiss

        integer :: nrow, ncol
        integer :: i, j, k
        integer, allocatable :: tmp_count(:)

        nrow = size(Bbasis_all, 1)
        ncol = size(Bbasis_all, 2)

        Bbasiss%nrow = nrow
        Bbasiss%ncol = ncol

        ! 1. First pass: Count non-zeros per row to find max_nz
        allocate(tmp_count(nrow))
        tmp_count = 0
        
        do i = 1, nrow
            do j = 1, ncol
                if (abs(Bbasis_all(i,j)) > eps_tol) then
                    tmp_count(i) = tmp_count(i) + 1
                end if
            end do
        end do
        
        Bbasiss%max_nz = maxval(tmp_count)
        
        ! specific check to prevent allocation error if matrix is all zeros
        if (Bbasiss%max_nz == 0) Bbasiss%max_nz = 1 

        ! 2. Allocate memory
        allocate(Bbasiss%nzcount(nrow))
        allocate(Bbasiss%colind(nrow, Bbasiss%max_nz))
        allocate(Bbasiss%val   (nrow, Bbasiss%max_nz))

        Bbasiss%nzcount = tmp_count
        Bbasiss%colind  = 0
        Bbasiss%val     = 0.0_dp

        ! 3. Second pass: Fill the structure
        do i = 1, nrow
            k = 0
            do j = 1, ncol
                if (abs(Bbasis_all(i,j)) > eps_tol) then
                    k = k + 1
                    Bbasiss%colind(i,k) = j
                    Bbasiss%val(i,k)    = Bbasis_all(i,j)
                end if
            end do
        end do

        deallocate(tmp_count)

    end subroutine build_Bbasis_sparse
    
    ! ==================================================================
! Subroutine: recover_Bbasis_dense
! Purpose: Reconstruct a dense matrix from the sparse representation.
!          (Inverse operation of build_Bbasis_sparse)
! ==================================================================
subroutine recover_Bbasis_dense(Bbasiss, Bbasis_all)
    implicit none
    
    ! Input: The sparse structure
    type(Bbasis_sparse_type), intent(in) :: Bbasiss
    
    ! Output: The dense matrix (Allocatable)
    real(dp), allocatable, intent(out) :: Bbasis_all(:,:)

    integer :: nrow, ncol
    integer :: i, k, col_idx
    real(dp) :: val

    ! 1. Retrieve dimensions from the sparse structure
    nrow = Bbasiss%nrow
    ncol = Bbasiss%ncol

    ! 2. Handle allocation
    ! If Bbasis_all is already allocated, deallocate it to ensure correct size
    if (allocated(Bbasis_all)) deallocate(Bbasis_all)
    allocate(Bbasis_all(nrow, ncol))

    ! 3. Initialize to Zero
    ! Crucial step: The sparse structure only stores non-zeros. 
    ! Everything else must be 0.0.
    Bbasis_all = 0.0_dp

    ! 4. Populate the dense matrix
    ! Iterate over each row
    do i = 1, nrow
        ! Iterate over the non-zero elements recorded for this row
        do k = 1, Bbasiss%nzcount(i)
            
            col_idx = Bbasiss%colind(i, k)  ! Source column index
            val     = Bbasiss%val(i, k)     ! Get the value
            
            ! Assign to the dense matrix
            Bbasis_all(i, col_idx) = val
            
        end do
    end do

end subroutine recover_Bbasis_dense
! ==================================================================
    ! Subroutine: build_sub_sparse
    ! Description: Extracts a sub-matrix based on row indices (Row Slicing).
    !              The output keeps the input column count.
    ! ==================================================================
subroutine build_sub_sparse(B_all_sparse, sampling_index, B_sub_sparse)
        implicit none
        type(Bbasis_sparse_type), intent(in)  :: B_all_sparse
        integer,                  intent(in)  :: sampling_index(:)
        type(Bbasis_sparse_type), intent(out) :: B_sub_sparse

        integer :: i, r, nz, max_nz_local, num_samplingnodes

        num_samplingnodes = size(sampling_index)

        ! 1. Set basic dimensions
        B_sub_sparse%nrow = num_samplingnodes
        B_sub_sparse%ncol = B_all_sparse%ncol

        ! 2. Determine local max_nz for the subset of rows
        !    This avoids excessive allocation when the selected
        !    rows are sparse compared to the global matrix.
        max_nz_local = 0
        do i = 1, num_samplingnodes
            r = sampling_index(i)
            
            ! Bounds check
            if (r > B_all_sparse%nrow .or. r < 1) then
                print *, "Error: Sampling index out of bounds at i=", i, " value=", r
                stop "build_sub_sparse: index error"
            end if
            
            nz = B_all_sparse%nzcount(r)
            if (nz > max_nz_local) max_nz_local = nz
        end do
        
        ! Ensure max_nz is at least 1 to avoid allocation errors
        if (max_nz_local == 0) max_nz_local = 1
        B_sub_sparse%max_nz = max_nz_local

        ! 3. Allocate arrays for the sub-matrix
        !    Note: Intent(out) usually deallocates automatically, but explicit
        !    allocation is standard practice here.
        if (allocated(B_sub_sparse%nzcount)) deallocate(B_sub_sparse%nzcount)
        if (allocated(B_sub_sparse%colind))  deallocate(B_sub_sparse%colind)
        if (allocated(B_sub_sparse%val))     deallocate(B_sub_sparse%val)

        allocate(B_sub_sparse%nzcount(num_samplingnodes))
        allocate(B_sub_sparse%colind(num_samplingnodes, max_nz_local))
        allocate(B_sub_sparse%val   (num_samplingnodes, max_nz_local))

        ! Initialize to zero
        B_sub_sparse%nzcount = 0
        B_sub_sparse%colind  = 0
        B_sub_sparse%val     = 0.0_dp

        ! 4. Copy data from global matrix to sub-matrix
        do i = 1, num_samplingnodes
            r = sampling_index(i)
            nz = B_all_sparse%nzcount(r)

            B_sub_sparse%nzcount(i) = nz
            
            if (nz > 0) then
                ! Fortran array slicing handles the copy efficiently
                B_sub_sparse%colind(i, 1:nz) = B_all_sparse%colind(r, 1:nz)
                B_sub_sparse%val(i, 1:nz)    = B_all_sparse%val(r, 1:nz)
            end if
        end do

    end subroutine build_sub_sparse

    ! ==================================================================
    ! Subroutine: spmv_Bbasis_sparse (Right Multiplication)
    ! Operation: y = A * x
    ! ==================================================================
    subroutine spmv_Bbasis_sparse(A, x, y)
        implicit none
        type(Bbasis_sparse_type), intent(in) :: A
        real(dp), intent(in)  :: x(:)     ! Input vector (length A%ncol)
        real(dp), intent(out) :: y(:)     ! Output vector (length A%nrow)

        integer :: i, k, nz, col

        ! Dimension Check
        if (size(x) /= A%ncol) stop "spmv_Bbasis_sparse: x dimension mismatch"
        if (size(y) /= A%nrow) stop "spmv_Bbasis_sparse: y dimension mismatch"

        y = 0.0_dp

        ! OpenMP parallelization for speedup on multi-core CPUs
        !$omp parallel do private(i, k, nz, col) shared(A, x, y)
        do i = 1, A%nrow
            nz = A%nzcount(i)
            do k = 1, nz
                col = A%colind(i,k)
                y(i) = y(i) + A%val(i,k) * x(col)
            end do
        end do
        !$omp end parallel do

end subroutine spmv_Bbasis_sparse

    ! ==================================================================
    ! Subroutine: spmv_Bbasis_sparse_transpose (Left Multiplication)
    ! Operation: y = x^T * A  (equivalent to y = A^T * x)
    ! ==================================================================
    subroutine spmv_Bbasis_sparse_transpose(A, x, y)
        implicit none
        type(Bbasis_sparse_type), intent(in) :: A
        real(dp), intent(in)  :: x(:)     ! Input vector (length A%nrow)
        real(dp), intent(out) :: y(:)     ! Output vector (length A%ncol)

        integer :: i, k, nz, col
        real(dp) :: val_x

        ! Dimension Check
        if (size(x) /= A%nrow) stop "spmv_sparse_transpose: x dimension mismatch"
        if (size(y) /= A%ncol) stop "spmv_sparse_transpose: y dimension mismatch"

        y = 0.0_dp

        ! Loop over rows of A and scatter the values to y
        ! Note: This creates dependencies on y, so simple OpenMP 'parallel for' 
        ! needs atomic operations or reduction, which might be slower.
        do i = 1, A%nrow
            val_x = x(i)
            if (abs(val_x) > eps_tol) then
                nz = A%nzcount(i)
                do k = 1, nz
                    col = A%colind(i,k)
                    y(col) = y(col) + val_x * A%val(i,k)
                end do
            end if
        end do

    end subroutine spmv_Bbasis_sparse_transpose

! ==================================================================
    ! Subroutine: make_csr_transpose
    ! Purpose:    Builds a CSR Transpose structure for fast column access.
    ! Inputs:     A (Sparse Matrix)
    ! Outputs:    at_rowptr, at_colind, at_val (CSR arrays representing A^T)
    ! ==================================================================
    subroutine make_csr_transpose(A, at_rowptr, at_colind, at_val)
        type(Bbasis_sparse_type), intent(in) :: A
        integer, allocatable, intent(out) :: at_rowptr(:), at_colind(:)
        real(dp), allocatable, intent(out) :: at_val(:)

        integer :: n, m, i, j, c, idx
        integer, allocatable :: count(:)

        n = A%ncol
        m = A%nrow

        ! 1. Count non-zeros per column (Row length of A^T)
        allocate(count(n))
        count = 0
        do i = 1, m
            do j = 1, A%nzcount(i)
                c = A%colind(i, j)
                if (c > 0) count(c) = count(c) + 1
            end do
        end do

        ! 2. Construct row pointers
        allocate(at_rowptr(n + 1))
        at_rowptr(1) = 1
        do i = 1, n
            at_rowptr(i+1) = at_rowptr(i) + count(i)
        end do

        ! 3. Fill column indices and values
        allocate(at_colind(at_rowptr(n+1)-1))
        allocate(at_val(at_rowptr(n+1)-1))
        
        count = 0 
        do i = 1, m
            do j = 1, A%nzcount(i)
                c = A%colind(i, j)
                if (c > 0) then
                    idx = at_rowptr(c) + count(c)
                    at_colind(idx) = i          
                    at_val(idx) = A%val(i, j)
                    count(c) = count(c) + 1
                end if
            end do
        end do
    end subroutine make_csr_transpose

    ! ==================================================================
    ! Subroutine: Sparse_MatMul_AB
    ! Purpose:    Computes C = A * B.
    ! Inputs:     A, B (Sparse Matrices)
    ! Outputs:    C_dense (Optional, 2D array) OR C_sparse (Optional, Type)
    ! ==================================================================
    subroutine Sparse_MatMul_AB(A, B, C_dense, C_sparse)
        type(Bbasis_sparse_type), intent(in) :: A, B
        real(dp), allocatable, intent(out), optional :: C_dense(:,:)
        type(Bbasis_sparse_type), intent(out), optional :: C_sparse

        integer :: m, n, k_dim
        integer :: i, j, k, col_B, k_idx
        real(dp) :: val_A
        
        ! SPA (Sparse Accumulator) variables
        real(dp), allocatable :: spa(:)
        integer, allocatable :: spa_marker(:)
        integer, allocatable :: idx_list(:) 
        integer :: list_cnt
        integer, allocatable :: row_nz_counts(:)
        integer :: max_nz_C

        ! Check dimensions
        if (A%ncol /= B%nrow) then
            print *, "Error: Sparse_MatMul_AB dimension mismatch"
            return
        end if

        m = A%nrow
        n = B%ncol
        k_dim = A%ncol

        ! Initialize SPA
        allocate(spa(n)); spa = zero
        allocate(spa_marker(n)); spa_marker = 0
        allocate(idx_list(n))

        ! --------------------------------------------------
        ! Path A: Output is Dense Matrix
        ! --------------------------------------------------
        if (present(C_dense)) then
            if (allocated(C_dense)) deallocate(C_dense)
            allocate(C_dense(m, n))
            C_dense = zero

            ! Compute C(i, :) = sum( A(i,k) * Row_k(B) )
            do i = 1, m
                list_cnt = 0
                do k_idx = 1, A%nzcount(i)
                    k = A%colind(i, k_idx)      
                    val_A = A%val(i, k_idx)

                    ! Accumulate row k of B into SPA
                    do j = 1, B%nzcount(k)
                        col_B = B%colind(k, j)
                        
                        if (spa_marker(col_B) /= i) then
                            spa_marker(col_B) = i
                            spa(col_B) = zero
                            list_cnt = list_cnt + 1
                            idx_list(list_cnt) = col_B
                        end if
                        spa(col_B) = spa(col_B) + val_A * B%val(k, j)
                    end do
                end do
                
                ! Scatter SPA to Dense Output
                do j = 1, list_cnt
                    col_B = idx_list(j)
                    C_dense(i, col_B) = spa(col_B)
                end do
            end do
            return 
        end if

        ! --------------------------------------------------
        ! Path B: Output is Sparse Matrix (Two-Pass Method)
        ! --------------------------------------------------
        if (present(C_sparse)) then
            allocate(row_nz_counts(m))
            
            ! Pass 1: Symbolic Analysis (Count NNZ per row)
            spa_marker = 0
            do i = 1, m
                list_cnt = 0
                do k_idx = 1, A%nzcount(i)
                    k = A%colind(i, k_idx)
                    do j = 1, B%nzcount(k)
                        col_B = B%colind(k, j)
                        if (spa_marker(col_B) /= i) then
                            spa_marker(col_B) = i
                            list_cnt = list_cnt + 1
                        end if
                    end do
                end do
                row_nz_counts(i) = list_cnt
            end do

            ! Allocate Sparse Structure
            max_nz_C = maxval(row_nz_counts)
            if (max_nz_C == 0) max_nz_C = 1

            C_sparse%nrow = m
            C_sparse%ncol = n
            C_sparse%max_nz = max_nz_C
            allocate(C_sparse%nzcount(m))
            allocate(C_sparse%colind(m, max_nz_C))
            allocate(C_sparse%val(m, max_nz_C))
            C_sparse%nzcount = row_nz_counts
            C_sparse%colind = 0
            C_sparse%val = zero

            ! Pass 2: Numeric Computation
            spa_marker = 0 
            do i = 1, m
                list_cnt = 0
                
                do k_idx = 1, A%nzcount(i)
                    k = A%colind(i, k_idx)
                    val_A = A%val(i, k_idx)

                    do j = 1, B%nzcount(k)
                        col_B = B%colind(k, j)
                        if (spa_marker(col_B) /= i) then
                            spa_marker(col_B) = i
                            spa(col_B) = zero
                            list_cnt = list_cnt + 1
                            idx_list(list_cnt) = col_B
                        end if
                        spa(col_B) = spa(col_B) + val_A * B%val(k, j)
                    end do
                end do

                ! Sort indices for deterministic sparse output
                call sort_indices(list_cnt, idx_list, spa)

                ! Fill Sparse Structure
                do j = 1, list_cnt
                    col_B = idx_list(j)
                    if (abs(spa(col_B)) > eps_tol) then 
                        C_sparse%colind(i, j) = col_B
                        C_sparse%val(i, j) = spa(col_B)
                    else
                         C_sparse%colind(i, j) = col_B 
                         C_sparse%val(i, j) = zero
                    endif
                end do
            end do
            deallocate(row_nz_counts)
        end if
        
        deallocate(spa, spa_marker, idx_list)

    end subroutine Sparse_MatMul_AB

    ! ==================================================================
    ! Subroutine: Sparse_MatMul_AtA
    ! Purpose:    Computes C = A^T * A. Exploits symmetry.
    ! Inputs:     A (Sparse Matrix)
    ! Outputs:    C_dense (Optional) OR C_sparse (Optional)
    ! ==================================================================
    subroutine Sparse_MatMul_AtA(A, C_dense, C_sparse)
        type(Bbasis_sparse_type), intent(in) :: A
        real(dp), allocatable, intent(out), optional :: C_dense(:,:)
        type(Bbasis_sparse_type), intent(out), optional :: C_sparse

        integer :: n, i, j, k, r, c, idx
        
        ! A^T Structure
        integer, allocatable :: at_rowptr(:), at_colind(:)
        real(dp), allocatable :: at_val(:)

        real(dp), allocatable :: spa(:)
        integer, allocatable :: spa_marker(:)
        integer, allocatable :: idx_list(:)
        integer :: list_cnt
        integer, allocatable :: row_nz_counts(:)
        integer :: max_nz_C
        real(dp) :: v_At 

        n = A%ncol 

        ! 1. Build A^T CSR structure for fast column access
        call make_csr_transpose(A, at_rowptr, at_colind, at_val)

        allocate(spa(n)); spa = zero
        allocate(spa_marker(n)); spa_marker = 0
        allocate(idx_list(n))

        ! --------------------------------------------------
        ! Path A: Output is Dense Matrix
        ! --------------------------------------------------
        if (present(C_dense)) then
            allocate(C_dense(n, n))
            C_dense = zero
            
            do i = 1, n
                list_cnt = 0
                ! Iterate over A^T row i (which is A col i)
                do k = at_rowptr(i), at_rowptr(i+1)-1
                    r = at_colind(k)     
                    v_At = at_val(k) 

                    ! Accumulate row r of A
                    do j = 1, A%nzcount(r)
                        c = A%colind(r, j)
                        ! Symmetry Optimization: Compute Upper Triangle Only
                        if (c >= i) then
                            if (spa_marker(c) /= i) then
                                spa_marker(c) = i
                                spa(c) = zero
                                list_cnt = list_cnt + 1
                                idx_list(list_cnt) = c
                            end if
                            spa(c) = spa(c) + v_At * A%val(r, j)
                        end if
                    end do
                end do
                
                ! Fill Upper and Mirror to Lower Triangle
                do j = 1, list_cnt
                    c = idx_list(j)
                    C_dense(i, c) = spa(c)
                    C_dense(c, i) = spa(c) 
                end do
            end do
        end if

        ! --------------------------------------------------
        ! Path B: Output is Sparse Matrix (Two-Pass Method)
        ! --------------------------------------------------
        if (present(C_sparse)) then
            allocate(row_nz_counts(n))
            
            ! Pass 1: Symbolic Analysis (Full matrix for simplicity in sparse structure)
            spa_marker = 0
            do i = 1, n
                list_cnt = 0
                do k = at_rowptr(i), at_rowptr(i+1)-1
                    r = at_colind(k)
                    do j = 1, A%nzcount(r)
                        c = A%colind(r, j)
                        if (spa_marker(c) /= i) then
                            spa_marker(c) = i
                            list_cnt = list_cnt + 1
                        end if
                    end do
                end do
                row_nz_counts(i) = list_cnt
            end do

            ! Allocation
            max_nz_C = maxval(row_nz_counts)
            if (max_nz_C==0) max_nz_C=1
            C_sparse%nrow = n; C_sparse%ncol = n; C_sparse%max_nz = max_nz_C
            allocate(C_sparse%nzcount(n))
            allocate(C_sparse%colind(n, max_nz_C))
            allocate(C_sparse%val(n, max_nz_C))
            C_sparse%nzcount = row_nz_counts
            C_sparse%val = zero

            ! Pass 2: Numeric Computation
            spa_marker = 0
            do i = 1, n
                list_cnt = 0
                do k = at_rowptr(i), at_rowptr(i+1)-1
                    r = at_colind(k)
                    v_At = at_val(k)
                    do j = 1, A%nzcount(r)
                        c = A%colind(r, j)
                        if (spa_marker(c) /= i) then
                            spa_marker(c) = i
                            spa(c) = zero
                            list_cnt = list_cnt + 1
                            idx_list(list_cnt) = c
                        end if
                        spa(c) = spa(c) + v_At * A%val(r, j)
                    end do
                end do
                
                call sort_indices(list_cnt, idx_list, spa)
                
                do j = 1, list_cnt
                    c = idx_list(j)
                    C_sparse%colind(i, j) = c
                    C_sparse%val(i, j) = spa(c)
                end do
            end do
            deallocate(row_nz_counts)
        end if

        deallocate(spa, spa_marker, idx_list)
        deallocate(at_rowptr, at_colind, at_val)

    end subroutine Sparse_MatMul_AtA

    ! ==================================================================
    ! Subroutine: Sparse_MatMul_AAt
    ! Purpose:    Computes C = A * A^T. Exploits symmetry.
    ! Inputs:     A (Sparse Matrix)
    ! Outputs:    C_dense (Optional) OR C_sparse (Optional)
    ! ==================================================================
    subroutine Sparse_MatMul_AAt(A, C_dense, C_sparse)
        type(Bbasis_sparse_type), intent(in) :: A
        real(dp), allocatable, intent(out), optional :: C_dense(:,:)
        type(Bbasis_sparse_type), intent(out), optional :: C_sparse

        integer :: m, i, j, k, k_idx, r, r_idx
        integer :: col_in_A
        real(dp) :: val_A_ik

        ! A^T Structure variables
        integer, allocatable :: at_rowptr(:), at_colind(:)
        real(dp), allocatable :: at_val(:)

        real(dp), allocatable :: spa(:)
        integer, allocatable :: spa_marker(:)
        integer, allocatable :: idx_list(:)
        integer :: list_cnt
        integer, allocatable :: row_nz_counts(:)
        integer :: max_nz_C

        m = A%nrow

        ! 1. Build A^T (Map: Column -> List of Rows)
        call make_csr_transpose(A, at_rowptr, at_colind, at_val)

        allocate(spa(m)); spa = zero
        allocate(spa_marker(m)); spa_marker = 0
        allocate(idx_list(m))

        ! --------------------------------------------------
        ! Path A: Output is Dense Matrix
        ! --------------------------------------------------
        if (present(C_dense)) then
            allocate(C_dense(m, m))
            C_dense = zero
            
            do i = 1, m
                list_cnt = 0
                ! C(i, :) = sum_k A(i,k) * Row_k(A^T)
                
                do k_idx = 1, A%nzcount(i)
                    col_in_A = A%colind(i, k_idx) 
                    val_A_ik = A%val(i, k_idx)    
                    
                    ! Look up which rows contain col_in_A using A^T structure
                    do r_idx = at_rowptr(col_in_A), at_rowptr(col_in_A+1)-1
                        r = at_colind(r_idx)     
                        
                        ! Symmetry Optimization: Upper Triangle Only
                        if (r >= i) then
                            if (spa_marker(r) /= i) then
                                spa_marker(r) = i
                                spa(r) = zero
                                list_cnt = list_cnt + 1
                                idx_list(list_cnt) = r
                            end if
                            ! C(i, r) += A(i,k) * A(r,k)
                            spa(r) = spa(r) + val_A_ik * at_val(r_idx)
                        end if
                    end do
                end do
                
                ! Fill Upper and Mirror to Lower
                do j = 1, list_cnt
                    r = idx_list(j)
                    C_dense(i, r) = spa(r)
                    C_dense(r, i) = spa(r)
                end do
            end do
        end if

        ! --------------------------------------------------
        ! Path B: Output is Sparse Matrix (Two-Pass Method)
        ! --------------------------------------------------
        if (present(C_sparse)) then
            allocate(row_nz_counts(m))
            
            ! Pass 1: Symbolic Analysis
            spa_marker = 0
            do i = 1, m
                list_cnt = 0
                do k_idx = 1, A%nzcount(i)
                    col_in_A = A%colind(i, k_idx)
                    do r_idx = at_rowptr(col_in_A), at_rowptr(col_in_A+1)-1
                        r = at_colind(r_idx)
                        if (spa_marker(r) /= i) then
                            spa_marker(r) = i
                            list_cnt = list_cnt + 1
                        end if
                    end do
                end do
                row_nz_counts(i) = list_cnt
            end do

            ! Allocation
            max_nz_C = maxval(row_nz_counts)
            if (max_nz_C==0) max_nz_C=1
            C_sparse%nrow = m; C_sparse%ncol = m; C_sparse%max_nz = max_nz_C
            allocate(C_sparse%nzcount(m))
            allocate(C_sparse%colind(m, max_nz_C))
            allocate(C_sparse%val(m, max_nz_C))
            C_sparse%nzcount = row_nz_counts
            C_sparse%val = zero

            ! Pass 2: Numeric Computation
            spa_marker = 0
            do i = 1, m
                list_cnt = 0
                do k_idx = 1, A%nzcount(i)
                    col_in_A = A%colind(i, k_idx)
                    val_A_ik = A%val(i, k_idx)
                    
                    do r_idx = at_rowptr(col_in_A), at_rowptr(col_in_A+1)-1
                        r = at_colind(r_idx)
                        if (spa_marker(r) /= i) then
                            spa_marker(r) = i
                            spa(r) = zero
                            list_cnt = list_cnt + 1
                            idx_list(list_cnt) = r
                        end if
                        spa(r) = spa(r) + val_A_ik * at_val(r_idx)
                    end do
                end do
                
                call sort_indices(list_cnt, idx_list, spa)
                
                do j = 1, list_cnt
                    r = idx_list(j)
                    C_sparse%colind(i, j) = r
                    C_sparse%val(i, j) = spa(r)
                end do
            end do
            deallocate(row_nz_counts)
        end if

        deallocate(spa, spa_marker, idx_list)
        deallocate(at_rowptr, at_colind, at_val)

    end subroutine Sparse_MatMul_AAt

    ! ==================================================================
    ! Subroutine: sort_indices
    ! Purpose:    Insertion sort to order column indices.
    ! Inputs:     inds (indices), vals (values - read only)
    ! ==================================================================
    subroutine sort_indices(n, inds, vals)
        integer, intent(in) :: n
        integer, intent(inout) :: inds(:)
        real(dp), intent(in) :: vals(:) 
        integer :: i, j, temp
        
        do i = 2, n
            temp = inds(i)
            j = i - 1
            do while (j >= 1)
                if (inds(j) <= temp) exit
                inds(j+1) = inds(j)
                j = j - 1
            end do
            inds(j+1) = temp
        end do
    end subroutine sort_indices
    
    ! ==================================================================
    ! Subroutine: free_Bbasis_sparse
    ! Description: Deallocates memory within the sparse structure.
    ! ==================================================================
    subroutine free_Bbasis_sparse(S)
        implicit none
        type(Bbasis_sparse_type), intent(inout) :: S
        
        if (allocated(S%nzcount)) deallocate(S%nzcount)
        if (allocated(S%colind))  deallocate(S%colind)
        if (allocated(S%val))     deallocate(S%val)
        
        S%nrow = 0
        S%ncol = 0
        S%max_nz = 0
    end subroutine free_Bbasis_sparse

! ==================================================================
    ! Subroutine: mpi_bcast_sparse
    ! Description: Efficiently broadcasts a sparse matrix structure
    !              without converting to dense format.
    ! ==================================================================
    subroutine mpi_bcast_sparse(S, root, comm)
        use mpi
        implicit none
        type(Bbasis_sparse_type), intent(inout) :: S
        integer, intent(in) :: root
        integer, intent(in) :: comm
        
        integer :: ierr, my_rank
        
        call MPI_Comm_rank(comm, my_rank, ierr)

        ! 1. Broadcast Metadata (Scalars)
        !    MPI requires separate calls or packing. Doing separate for clarity.
        call MPI_Bcast(S%nrow,   1, MPI_INTEGER, root, comm, ierr)
        call MPI_Bcast(S%ncol,   1, MPI_INTEGER, root, comm, ierr)
        call MPI_Bcast(S%max_nz, 1, MPI_INTEGER, root, comm, ierr)

        ! 2. Allocate memory on non-root processes
        if (my_rank /= root) then
            ! Clean up if already allocated
            if (allocated(S%nzcount)) deallocate(S%nzcount)
            if (allocated(S%colind))  deallocate(S%colind)
            if (allocated(S%val))     deallocate(S%val)

            ! Allocate based on received dimensions
            allocate(S%nzcount(S%nrow))
            ! Handle case where max_nz might be 0 (empty matrix)
            if (S%max_nz > 0) then
                allocate(S%colind(S%nrow, S%max_nz))
                allocate(S%val(S%nrow, S%max_nz))
            else
                ! Allocate size 1 to avoid issues, though it won't be used
                allocate(S%colind(S%nrow, 1))
                allocate(S%val(S%nrow, 1))
            end if
        end if

        ! 3. Broadcast the Arrays
        !    nzcount (Size: nrow)
        call MPI_Bcast(S%nzcount, S%nrow, MPI_INTEGER, root, comm, ierr)

        if (S%max_nz > 0) then
            ! colind (Size: nrow * max_nz)
            call MPI_Bcast(S%colind, S%nrow * S%max_nz, MPI_INTEGER, root, comm, ierr)
            
            ! val (Size: nrow * max_nz)
            call MPI_Bcast(S%val,    S%nrow * S%max_nz, MPI_DOUBLE_PRECISION, root, comm, ierr)
        end if

    end subroutine mpi_bcast_sparse    
end module SparseMatrix_mod
