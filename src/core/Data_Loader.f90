module Data_Loader
    use mpi
    use GlobalSettings_mod
    use mMSTD_mod
    use Bsplines_mod
    use SamplesGeneration_mod
    use utils_mod          
    use SparseMatrix_mod   
    implicit none

    public 

    contains
    
subroutine Prepare_Spatial_Neighbors(ranks, output_dir, NeighborDis, NeighborId,Nodesset)
    implicit none
    
    ! --- Inputs ---
    integer, intent(in) :: ranks
    character(len=*), intent(in) :: output_dir
    
    ! --- Outputs ---
    real(dp), allocatable, intent(out) :: NeighborDis(:,:)
    integer, allocatable, intent(out) :: NeighborId(:,:)
    
    ! --- Optional Input (used for generation if needed) ---
    real(dp), allocatable, intent(in), optional :: Nodesset(:,:)
    
    ! --- Local Variables ---
    integer :: ierr, stat
    logical :: exist_dis, exist_id
    character(len=256) :: file_path_dis, file_path_id
    real(dp), allocatable :: Dense_Tmp(:,:)
    
    ! 1. Allocate Memory
    allocate(NeighborDis(num_allnodes, num_disnearspatial), stat=stat)
    allocate(NeighborId(num_allnodes, num_disnearspatial), stat=stat)
    if (stat /= 0) call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
    
    file_path_dis = trim(output_dir)//"/NeighborDis_m0.dat"
    file_path_id  = trim(output_dir)//"/NeighborId_m0.dat"

    ! 2. Rank 0 handles Logic
    if (ranks == 0) then
        ! Check existence of m0 files
        inquire(file=file_path_dis, exist=exist_dis)
        inquire(file=file_path_id,  exist=exist_id)
        
        ! 3. Compute if missing
        if (.not. (exist_dis .and. exist_id)) then
            write(*,'(a)') "Computing Mode 0 (Euclidean) Neighbors..."
            
            ! Ensure we have Nodes
            if (present(Nodesset)) then
                call NearestNeighbors(Nodesset, num_disnearspatial, 0, NeighborDis, NeighborId)
            else
                allocate(Dense_Tmp(num_allnodes,2))
                call load_dense_array(trim(output_dir)//"/Nodesset.dat", Dense_Tmp)
                call NearestNeighbors(Dense_Tmp, num_disnearspatial, 0, NeighborDis, NeighborId)
                deallocate(Dense_Tmp)
            end if
            
            ! Save newly computed m0 files
            call save_dense_array(file_path_dis, NeighborDis)
            call save_dense_array_int(file_path_id, NeighborId)
        else
            ! 4. Load existing files
            write(*,'(a)') "Loading Mode 0 Neighbors..."
            call load_dense_array(file_path_dis, NeighborDis)
            call load_dense_array_int(file_path_id, NeighborId)
        end if
    end if

    ! 5. Broadcast to all ranks
    call MPI_Bcast(NeighborDis, size(NeighborDis), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(NeighborId, size(NeighborId), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

end subroutine Prepare_Spatial_Neighbors  

    !====================================================================
    ! Subroutine: Prepare_Simulation_Data
    ! Purpose: Allocates, Loads/Computes, and Broadcasts simulation data.
    !          Supports separate file I/O and fully optional arguments.
    !====================================================================
    subroutine Prepare_Simulation_Data(ranks, output_dir, &
    ! --- Basic Outputs (Optional) ---
    Nodesset, NeighborDis, NeighborId, &
    Kernel_NeighborDis, Kernel_NeighborId, &
    ! --- B-Spline Basis (Optional) ---
    B0_out, B1_out, &
    ! --- Covariance Matrices & Neighbors (Optional) ---
    CovMat_out, &
    ! --- Kernel Method Matrices & Neighbors (Optional) ---
    Kernel_Cov_NeighborId, Kernel_Cov_NeighborVal, Kernel_Cov_diag, Cov_dense_global)

    implicit none
    
    ! --- Inputs ---
    integer, intent(in) :: ranks
    character(len=*), intent(in) :: output_dir 
    
    ! --- Outputs ---
    real(dp), allocatable, intent(out), optional :: Nodesset(:,:)
    real(dp), allocatable, intent(out), optional :: NeighborDis(:,:), Kernel_NeighborDis(:,:)
    integer, allocatable, intent(out), optional :: NeighborId(:,:), Kernel_NeighborId(:,:)

    ! --- Optional Sparse Structures ---
    type(Bbasis_sparse_type), intent(out), optional :: B0_out, B1_out
    type(Bbasis_sparse_type), intent(out), optional :: CovMat_out

    ! --- Kernel Method group ---
    integer, allocatable, intent(out), optional  :: Kernel_Cov_NeighborId(:,:)
    real(dp), allocatable, intent(out), optional :: Kernel_Cov_NeighborVal(:,:)
    real(dp), allocatable, intent(out), optional :: Kernel_Cov_diag(:)
    real(dp), allocatable, intent(out), optional :: Cov_dense_global(:,:) 
    real(dp), allocatable :: Cov_dense_local(:,:)
    real(dp), allocatable :: temp_diag(:,:)



    integer :: i, j, k, max_idx
    real(dp) :: max_val

    ! --- Local Variables --- 
    type(Bbasis_sparse_type) :: Cov0_Local, KernelCov_Local
    logical :: need_kcov_mem, file_exists, exist_id, exist_val,exist_dense,exist_diag 
    integer :: ierr, stat
    character(len=256) :: file_path_dis, file_path_dense, file_path_id, file_path_base, file_path_diag, path_kcov
    real(dp), allocatable :: Dense_Tmp(:,:)
    real(dp), allocatable :: Nodes_Tmp(:,:)

    !================================================================
    ! 0. Nodesset (Optional)
    !================================================================
    if (present(Nodesset)) then
        allocate(Nodesset(num_allnodes,2), stat=stat)
        if (stat/=0) call MPI_Abort(MPI_COMM_WORLD, 1, ierr)

        file_path_base = trim(output_dir) // "/Nodesset.dat"
        if (ranks==0) then
            inquire(file=file_path_base, exist=file_exists)
            if (file_exists) then
                write(*,'(a)') "Loading Nodesset..."
                call load_dense_array(file_path_base, Nodesset)
            else
                write(*,'(a)') "Computing Nodesset..."
                call GeneratePoints(Nodesset)
                call save_dense_array(file_path_base, Nodesset)
            end if
        end if
        call MPI_Bcast(Nodesset, size(Nodesset), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    end if

    !================================================================
    ! 1. Neighbors (Spatial) - Default Mode 1 Logic
    !================================================================
    if (present(NeighborDis) .and. present(NeighborId)) then
        allocate(NeighborDis(num_allnodes,num_disnearspatial), stat=stat)
        allocate(NeighborId(num_allnodes,num_disnearspatial), stat=stat)
        if (stat/=0) call MPI_Abort(MPI_COMM_WORLD, 2, ierr)

        file_path_dis = trim(output_dir)//"/NeighborDis_m1.dat"
        file_path_id  = trim(output_dir)//"/NeighborId_m1.dat"

        if (ranks == 0) then
            inquire(file=file_path_dis, exist=exist_val)
            inquire(file=file_path_id,  exist=exist_id)
            
            ! If files missing, compute Mode 1
            if (.not. (exist_val .and. exist_id)) then
                ! Need Nodes
                if (present(Nodesset)) then
                    call NearestNeighbors(Nodesset, num_disnearspatial, 1, NeighborDis, NeighborId)
                else
                    allocate(Dense_Tmp(num_allnodes,2))
                    call load_dense_array(trim(output_dir)//"/Nodesset.dat", Dense_Tmp)
                    call NearestNeighbors(Dense_Tmp, num_disnearspatial, 1, NeighborDis, NeighborId)
                    deallocate(Dense_Tmp)
                end if
                
                call save_dense_array(file_path_dis, NeighborDis)
                call save_dense_array_int(file_path_id, NeighborId)
            else
                ! Load Mode 1
                call load_dense_array(file_path_dis, NeighborDis)
                call load_dense_array_int(file_path_id, NeighborId)
            end if
        end if
        call MPI_Bcast(NeighborDis, size(NeighborDis), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(NeighborId, size(NeighborId), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    end if

    !================================================================
    ! 2. Kernels (Optional)
    !================================================================
    if (present(Kernel_NeighborDis) .and. present(Kernel_NeighborId)) then
        allocate(Kernel_NeighborDis(num_allnodes,num_kelnear_mM), stat=stat)
        allocate(Kernel_NeighborId(num_allnodes,num_kelnear_mM), stat=stat)
        
        file_path_dis = trim(output_dir) // "/Kernel_NeighborDis.dat"
        file_path_id  = trim(output_dir) // "/Kernel_NeighborId.dat"
        
        if (ranks==0) then
            inquire(file=file_path_dis, exist=file_exists)
            if (file_exists) then
                call load_dense_array(file_path_dis, Kernel_NeighborDis)
                call load_dense_array_int(file_path_id, Kernel_NeighborId)
            else
                if (present(Nodesset)) then
                    call EpaKernelTopK(Nodesset, Epah_mM, num_kelnear_mM, Kernel_NeighborDis, Kernel_NeighborId)
                else
                    allocate(Dense_Tmp(num_allnodes,2))
                    call load_dense_array(trim(output_dir)//"/Nodesset.dat", Dense_Tmp)
                    call EpaKernelTopK(Dense_Tmp, Epah_mM, num_kelnear_mM, Kernel_NeighborDis, Kernel_NeighborId)
                    deallocate(Dense_Tmp)
                end if
                call save_dense_array(file_path_dis, Kernel_NeighborDis)
                call save_dense_array_int(file_path_id, Kernel_NeighborId)
            end if
        end if
        call MPI_Bcast(Kernel_NeighborDis, size(Kernel_NeighborDis), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(Kernel_NeighborId, size(Kernel_NeighborId), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    end if

    !================================================================
    ! 3. B0 Basis (Optional)
    !================================================================
    if (present(B0_out)) then
        file_path_base = trim(output_dir) // "/B0_sparse.dat"
        if (ranks == 0) then
            inquire(file=file_path_base, exist=file_exists)
            if (file_exists) then
                call load_sparse_type(file_path_base, B0_out)
            else
                allocate(Dense_Tmp(num_allnodes, knot0_square_TS))
                call gen_bspline2D(num_x, num_y, kx0_TS, ky0_TS, degx0_TS, degy0_TS, Dense_Tmp)
                call build_Bbasis_sparse(Dense_Tmp, B0_out)
                deallocate(Dense_Tmp)
                call save_sparse_type(file_path_base, B0_out)
            end if
        end if
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        call mpi_bcast_sparse(B0_out, 0, MPI_COMM_WORLD)
    end if

    !================================================================
    ! 4. B1 Basis (Optional)
    !================================================================
    if (present(B1_out)) then
        file_path_base = trim(output_dir) // "/B1_sparse.dat"
        if (ranks == 0) then
            inquire(file=file_path_base, exist=file_exists)
            if (file_exists) then
                call load_sparse_type(file_path_base, B1_out)
            else
                allocate(Dense_Tmp(num_allnodes, knot1_square_TS))
                call gen_bspline2D(num_x, num_y, kx1_TS, ky1_TS, degx1_TS, degy1_TS, Dense_Tmp)
                call build_Bbasis_sparse(Dense_Tmp, B1_out)
                deallocate(Dense_Tmp)
                call save_sparse_type(file_path_base, B1_out)
            end if
        end if
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        call mpi_bcast_sparse(B1_out, 0, MPI_COMM_WORLD)
    end if

    !================================================================
    ! 5. Global CovMat (Optional)
    !================================================================
    if (present(CovMat_out)) then
        file_path_base = trim(output_dir) // "/CovMat_sparse.dat"
        if (ranks == 0) then
            inquire(file=file_path_base, exist=file_exists)
            if (file_exists) then
                call load_sparse_type(file_path_base, CovMat_out)
            else
                allocate(Dense_Tmp(num_allnodes, num_allnodes))
                if(present(Nodesset)) then
                    call ComputeCovMatrix(Nodesset, Dense_Tmp)
                else
                    allocate(Nodes_Tmp(num_allnodes,2))
                    call load_dense_array(trim(output_dir)//"/Nodesset.dat", Nodes_Tmp)
                    call ComputeCovMatrix(Nodes_Tmp, Dense_Tmp)
                    deallocate(Nodes_Tmp)
                end if
                call build_Bbasis_sparse(Dense_Tmp, CovMat_out)
                deallocate(Dense_Tmp)
                call save_sparse_type(file_path_base, CovMat_out)
            end if
        end if
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        call mpi_bcast_sparse(CovMat_out, 0, MPI_COMM_WORLD)
    end if

!================================================================
    ! 6. Kernel Method (Optimized Dense Matrix & Top-K Extraction)
    !================================================================
    if (present(Kernel_Cov_NeighborId) .or. present(Kernel_Cov_NeighborVal) .or. &
        present(Kernel_Cov_diag)) then
        
        ! Allocate output arrays if not present
        if (present(Kernel_Cov_NeighborId) .and. .not. allocated(Kernel_Cov_NeighborId)) &
            allocate(Kernel_Cov_NeighborId(num_allnodes, num_varnear_CDS))
        if (present(Kernel_Cov_NeighborVal) .and. .not. allocated(Kernel_Cov_NeighborVal)) &
            allocate(Kernel_Cov_NeighborVal(num_allnodes, num_varnear_CDS))
        if (present(Kernel_Cov_diag) .and. .not. allocated(Kernel_Cov_diag)) &
            allocate(Kernel_Cov_diag(num_allnodes))
            
        if (.not. allocated(Cov_dense_global)) &
            allocate(Cov_dense_global(num_allnodes, num_allnodes))

        ! Define cache file paths
        file_path_id    = trim(output_dir) // "/KernelCov_NeighborId.dat"
        file_path_dis   = trim(output_dir) // "/KernelCov_NeighborVal.dat"
        file_path_diag  = trim(output_dir) // "/KernelCov_diag.dat"
        file_path_dense = trim(output_dir) // "/KernelCov_dense.dat" 
        
        if (ranks == 0) then
            need_kcov_mem = .false.
            
            ! Check if cache files exist
            if (present(Kernel_Cov_NeighborId)) then
                inquire(file=file_path_id, exist=exist_id)
                inquire(file=file_path_dis, exist=exist_val)
                if (.not. (exist_id .and. exist_val)) need_kcov_mem = .true.
            end if
            if (present(Kernel_Cov_diag)) then
                inquire(file=file_path_diag, exist=exist_diag)
                if (.not. exist_diag) need_kcov_mem = .true.
            end if
            

            inquire(file=file_path_dense, exist=exist_dense)
            if (.not. exist_dense) need_kcov_mem = .true.

            if (need_kcov_mem) then    
    
                call Compute_Kernel_Covariance(Cov_dense_global)
                

                allocate(temp_diag(num_allnodes, 1))
                do i = 1, num_allnodes
                    Cov_dense_global(i, i) = Cov_dense_global(i, i) + sigma_noise**2.0d0
                    
                    if (present(Kernel_Cov_diag)) then
                        Kernel_Cov_diag(i) = Cov_dense_global(i, i)
                        temp_diag(i, 1) = Kernel_Cov_diag(i)
                    end if
                end do
                
    
                call save_dense_array(file_path_dense, Cov_dense_global)
                if (present(Kernel_Cov_diag)) then
                    call save_dense_array(file_path_diag, temp_diag)
                end if
                deallocate(temp_diag)

                if (present(Kernel_Cov_NeighborId) .and. present(Kernel_Cov_NeighborVal)) then
                    allocate(Cov_dense_local(num_allnodes, num_allnodes))
                    Cov_dense_local = Cov_dense_global

                    do i = 1, num_allnodes
                        ! Mask self-correlation to exclude it from neighbors
                        Cov_dense_local(i, i) = -huge(1.0_dp)
                        
                        do k = 1, num_varnear_CDS
                            max_val = -huge(1.0_dp)
                            max_idx = 0
                            
                            ! Search column 'i' for max covariance
                            do j = 1, num_allnodes
                                if (Cov_dense_local(j, i) > max_val) then
                                    max_val = Cov_dense_local(j, i)
                                    max_idx = j
                                end if
                            end do
                            
                            ! Record the top neighbor
                            Kernel_Cov_NeighborId(i, k)  = max_idx
                            Kernel_Cov_NeighborVal(i, k) = max_val
                            
                            ! Mask the found neighbor to prevent re-selection
                            Cov_dense_local(max_idx, i) = -huge(1.0_dp)
                        end do
                    end do
                    
                    call save_dense_array_int(file_path_id, Kernel_Cov_NeighborId)
                    call save_dense_array(file_path_dis, Kernel_Cov_NeighborVal)
                    deallocate(Cov_dense_local)
                end if
            else
                ! Load existing data from cache
                if (present(Kernel_Cov_NeighborId)) then
                    call load_dense_array_int(file_path_id, Kernel_Cov_NeighborId)
                    call load_dense_array(file_path_dis, Kernel_Cov_NeighborVal)
                end if
                if (present(Kernel_Cov_diag)) then
                    ! Load 2D temp array and map back to 1D
                    allocate(temp_diag(num_allnodes, 1))
                    call load_dense_array(file_path_diag, temp_diag)
                    Kernel_Cov_diag(:) = temp_diag(:, 1)
                    deallocate(temp_diag)
                end if

                call load_dense_array(file_path_dense, Cov_dense_global)
            end if
        end if


        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        
        if (present(Kernel_Cov_NeighborId)) then
            call MPI_Bcast(Kernel_Cov_NeighborId, size(Kernel_Cov_NeighborId), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(Kernel_Cov_NeighborVal, size(Kernel_Cov_NeighborVal), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        end if
        
        if (present(Kernel_Cov_diag)) then
            call MPI_Bcast(Kernel_Cov_diag, size(Kernel_Cov_diag), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        end if
        

        call MPI_Bcast(Cov_dense_global, size(Cov_dense_global), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    end if



end subroutine Prepare_Simulation_Data

 
!====================================================================
    ! Subroutine: compute_and_save_diag
    ! Purpose: Extracts diagonal from Sparse S, applies scaling, adds noise_var.
    ! Formula: OutDiag(i) = noise_var + scale_factor * S(i,i)
    !====================================================================
    subroutine compute_and_save_diag(S, filename, noise_var,  OutDiag)
        implicit none
        type(Bbasis_sparse_type), intent(in) :: S
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: noise_var 
        real(dp), intent(out) :: OutDiag(:)
        
        integer :: i, k
        logical :: file_exists
        
        inquire(file=filename, exist=file_exists)
        
        if (file_exists) then
            write(*,'(a)') "Loading Diagonal from " // trim(filename) // "..."
            call load_dense_vec(filename, OutDiag)
        else
            write(*,'(a)') "Computing Diagonal (Scaled) + Noise..."

            OutDiag = noise_var 

            do i = 1, S%nrow
                do k = 1, S%nzcount(i)
                    if (S%colind(i, k) == i) then
                        OutDiag(i) = OutDiag(i) + S%val(i, k) 
                        exit
                    end if
                end do
            end do
            
            call save_dense_vec(filename, OutDiag)
        end if
    end subroutine compute_and_save_diag

    !====================================================================
    ! Helper: Vector I/O
    !====================================================================
    subroutine save_dense_vec(filename, Vec)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: Vec(:)
        integer :: unit=35
        open(unit=unit, file=filename, status='replace', action='write')
        write(unit, *) Vec
        close(unit)
    end subroutine save_dense_vec

    subroutine load_dense_vec(filename, Vec)
        character(len=*), intent(in) :: filename
        real(dp), intent(out) :: Vec(:)
        integer :: unit=35
        open(unit=unit, file=filename, status='old', action='read')
        read(unit, *) Vec
        close(unit)
    end subroutine load_dense_vec


    subroutine save_dense_array(filename, Arr)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: Arr(:,:)
        integer :: unit=30
        open(unit=unit, file=filename, status='replace', action='write')
        write(unit, *) Arr
        close(unit)
    end subroutine save_dense_array

    subroutine load_dense_array(filename, Arr)
        character(len=*), intent(in) :: filename
        real(dp), intent(out) :: Arr(:,:)
        integer :: unit=30
        open(unit=unit, file=filename, status='old', action='read')
        read(unit, *) Arr
        close(unit)
    end subroutine load_dense_array

    subroutine save_dense_array_int(filename, Arr)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: Arr(:,:)
        integer :: unit=30
        open(unit=unit, file=filename, status='replace', action='write')
        write(unit, *) Arr
        close(unit)
    end subroutine save_dense_array_int

    subroutine load_dense_array_int(filename, Arr)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: Arr(:,:)
        integer :: unit=30
        open(unit=unit, file=filename, status='old', action='read')
        read(unit, *) Arr
        close(unit)
    end subroutine load_dense_array_int

    subroutine save_sparse_type(filename, S)
        character(len=*), intent(in) :: filename
        type(Bbasis_sparse_type), intent(in) :: S
        integer :: unit_num, i
        if (.not. allocated(S%nzcount)) return 
        unit_num = 20
        open(unit=unit_num, file=filename, status='replace', action='write', form='formatted')
        write(unit_num, *) S%nrow, S%ncol, S%max_nz
        write(unit_num, *) S%nzcount
        do i = 1, S%nrow
            write(unit_num, *) S%colind(i, 1:S%nzcount(i))
        end do
        do i = 1, S%nrow
            write(unit_num, *) S%val(i, 1:S%nzcount(i))
        end do
        close(unit_num)
    end subroutine save_sparse_type

    subroutine load_sparse_type(filename, S)
        character(len=*), intent(in) :: filename
        type(Bbasis_sparse_type), intent(out) :: S
        integer :: unit_num, i, nr, nc, mnz
        unit_num = 21
        open(unit=unit_num, file=filename, status='old', action='read', form='formatted')
        read(unit_num, *) nr, nc, mnz
        S%nrow = nr; S%ncol = nc; S%max_nz = mnz
        if (allocated(S%nzcount)) deallocate(S%nzcount)
        if (allocated(S%colind))  deallocate(S%colind)
        if (allocated(S%val))     deallocate(S%val)
        allocate(S%nzcount(nr))
        allocate(S%colind(nr, mnz))
        allocate(S%val(nr, mnz))
        S%colind = 0; S%val = 0.0d0
        read(unit_num, *) S%nzcount
        do i = 1, nr
            if (S%nzcount(i) > 0) read(unit_num, *) S%colind(i, 1:S%nzcount(i))
        end do
        do i = 1, nr
            if (S%nzcount(i) > 0) read(unit_num, *) S%val(i, 1:S%nzcount(i))
        end do
        close(unit_num)
    end subroutine load_sparse_type

end module Data_Loader
