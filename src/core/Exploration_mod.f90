module Exploration_mod
    use GlobalSettings_mod
    use mMSTD_mod
    implicit none
    
    ! Set default visibility to private
    private 
    public :: Init_Exploration_Data, Reset_Exploration_State, &
              ExplorationEvaluation, Clean_Exploration_Data
    
    ! --- Internal Static Data ---
    ! Geometry and topology data
    real(dp), allocatable :: NeighborDis_m0(:,:)
    real(dp), allocatable :: Nodesset_EX(:,:)
    integer,  allocatable :: NeighborId_m0(:,:)
    
    ! Calculation state (Current min distance for each node)
    real(dp), allocatable :: Min_Disterm(:)

contains   

    ! =========================================================
    ! Initialize data: Allocate memory and copy geometry
    ! =========================================================
   subroutine Init_Exploration_Data(NeighborDis_m0_in, NeighborId_m0_in, Nodesset_in)
        implicit none
        real(dp), intent(in) :: NeighborDis_m0_in(:,:)
        integer, intent(in)  :: NeighborId_m0_in(:,:)
        real(dp), intent(in) :: Nodesset_in(:,:)
        integer :: err

        ! 1. Free existing memory if allocated
        if (allocated(NeighborDis_m0)) deallocate(NeighborDis_m0)
        if (allocated(NeighborId_m0))  deallocate(NeighborId_m0)
        if (allocated(Min_Disterm))    deallocate(Min_Disterm)
        if (allocated(Nodesset_EX))  deallocate(Nodesset_EX)
        ! 2. Allocate memory
        allocate(NeighborDis_m0(num_allnodes, num_disnearspatial), stat=err)
        allocate(NeighborId_m0(num_allnodes, num_disnearspatial), stat=err)
        allocate(Nodesset_EX(num_allnodes, 2), stat=err)
        allocate(Min_Disterm(num_allnodes), stat=err)

        if (err /= 0) then
            print *, "Error: Memory allocation failed in Init_Exploration_Data"
            stop
        end if

        ! 3. Copy data and initialize state
        NeighborDis_m0 = NeighborDis_m0_in
        NeighborId_m0  = NeighborId_m0_in
        Min_Disterm    = huge(1.0_dp)
        Nodesset_EX=Nodesset_in

    end subroutine Init_Exploration_Data

    ! =========================================================
    ! Reset state: Set distances to infinity without reallocation
    ! =========================================================
    subroutine Reset_Exploration_State()
        implicit none
        
        if (.not. allocated(Min_Disterm)) then
            allocate(Min_Disterm(num_allnodes))
        end if
        
        Min_Disterm = huge(1.0_dp) 
    end subroutine Reset_Exploration_State

    ! =========================================================
    ! Clean up: Free memory
    ! =========================================================
    subroutine Clean_Exploration_Data()
        implicit none
        if (allocated(NeighborDis_m0)) deallocate(NeighborDis_m0)
        if (allocated(NeighborId_m0))  deallocate(NeighborId_m0)
        if (allocated(Min_Disterm))    deallocate(Min_Disterm)
    end subroutine Clean_Exploration_Data

    ! =========================================================
    ! Core Evaluation: Update distances and return MaxMin value
    ! =========================================================
    subroutine ExplorationEvaluation(sampling_index, maxmin_dis)
    implicit none
    
    integer,  intent(in)  :: sampling_index(:) 
    real(dp), intent(out) :: maxmin_dis

    ! Local variables
    integer  :: i, j, k, s
    integer  :: pick_idx, neighbor_idx, near_index
    integer  :: n_samples
    real(dp) :: temp_val, dx, dy, dist_sq, cur_min_dist

    logical, allocatable :: is_visited(:)

    n_samples = size(sampling_index)
    

    allocate(is_visited(num_allnodes))
    is_visited = .false.


    do i = 1, n_samples
        pick_idx = sampling_index(i)
        near_index = 0
        
        do j = 1, num_disnearspatial
            neighbor_idx = NeighborId_m0(pick_idx, j)

            is_visited(neighbor_idx) = .true.
            
            temp_val = NeighborDis_m0(pick_idx, j)

            if (temp_val < Min_Disterm(neighbor_idx)) then
                Min_Disterm(neighbor_idx) = temp_val
                near_index = 0 
            else
                near_index = near_index + 1
            end if

            if (near_index > near_Threshold_mM) exit
        end do
    end do


    do k = 1, num_allnodes
        if (.not. is_visited(k)) then
            cur_min_dist = huge(1.0_dp)
            
            do s = 1, n_samples
                pick_idx = sampling_index(s)

                dx = abs(Nodesset_EX(k,1) - Nodesset_EX(pick_idx,1))
                dy = abs(Nodesset_EX(k,2) - Nodesset_EX(pick_idx,2))
                dist_sq = dx*dx + dy*dy                
                if (dist_sq < cur_min_dist) cur_min_dist = dist_sq
            end do
            if (cur_min_dist < Min_Disterm(k)) then
                Min_Disterm(k) = cur_min_dist
            end if
            
        end if
    end do

    maxmin_dis = maxval(Min_Disterm)
    
    deallocate(is_visited)

end subroutine ExplorationEvaluation

end module Exploration_mod