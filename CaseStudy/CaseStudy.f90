module constants
    use iso_fortran_env, only: int64
    implicit none
    ! Case-study grid and statistic constants.
    integer, parameter :: num_x=232,num_y=292, Time_window = 299
    integer, parameter :: num_allnode=num_x*num_y
    real(8), parameter :: x_step=0.01d0, y_step=0.01d0, miu_min=1.0d0
    real(8), parameter :: cuparame=(miu_min**2.0d0)/2.0d0
    ! Grid and allocation constants.
    real(8), parameter :: Kernel_threshold=0.0d0
    integer, parameter ::center_point=100, near_Threshold = 10,Top_r = 1
!=================================================================================================================================
    real(8), parameter :: gone=1.0d0/dble(num_allnode)

contains

    subroutine load_normalized_data(path, norm_data)
        character(len=*), intent(in) :: path
        real(8), intent(out) :: norm_data(num_allnode, Time_window)
        integer :: io_status, unit_in
        integer(int64) :: actual_size, expected_size

        expected_size = int(num_allnode, int64) * int(Time_window, int64) * 8_int64
        inquire(file=trim(path), size=actual_size, iostat=io_status)
        if (io_status /= 0 .or. actual_size /= expected_size) then
            write(*, '(A)') 'Invalid or missing normalized-data cache: ' // trim(path)
            write(*, '(A,I0)') 'Expected bytes: ', expected_size
            if (io_status == 0) write(*, '(A,I0)') 'Actual bytes:   ', actual_size
            error stop 2
        end if

        open(newunit=unit_in, file=trim(path), status='old', action='read', &
             access='stream', form='unformatted', iostat=io_status)
        if (io_status /= 0) then
            write(*, '(A)') 'Unable to open normalized-data cache: ' // trim(path)
            error stop 2
        end if
        read(unit_in, iostat=io_status) norm_data
        close(unit_in)
        if (io_status /= 0) then
            write(*, '(A)') 'Unable to read normalized-data cache: ' // trim(path)
            error stop 2
        end if
        write(*, '(A)') 'Loaded normalized solar data from ' // trim(path)
    end subroutine load_normalized_data

    subroutine write_statistics(output_dir, method, statistics, wide)
        character(len=*), intent(in) :: output_dir, method
        real(8), intent(in) :: statistics(Time_window)
        logical, intent(in), optional :: wide
        character(len=2048) :: path
        integer :: i, io_status, unit_out
        logical :: use_wide

        path = trim(output_dir) // '/' // trim(method) // '_ChartingStatistics-5.txt'
        open(newunit=unit_out, file=trim(path), status='replace', action='write', &
             iostat=io_status)
        if (io_status /= 0) then
            write(*, '(A)') 'Unable to create result file: ' // trim(path)
            error stop 3
        end if
        use_wide = .false.
        if (present(wide)) use_wide = wide
        do i = 1, Time_window
            if (use_wide) then
                write(unit_out, '(F20.4)') statistics(i)
            else
                write(unit_out, '(F10.4)') statistics(i)
            end if
        end do
        close(unit_out)
        write(*, '(A)') 'Wrote ' // trim(path)
    end subroutine write_statistics

end module constants

program main
    use constants
    use RNSET_INT
    implicit none
    real(8), allocatable :: norm_data(:,:)
    real(8), allocatable :: NeighborDis(:), Kernel_NeighborDis(:)
    real(8), allocatable :: Nodesset(:,:), Statistics(:)
    integer, allocatable :: NeighborId(:), Kernel_NeighborId(:)
    integer :: num_samplingnodes, seed, arg_status, allocation_status
    real(8) :: Limit
    real(8) :: pho, sigmaprio, lambda, Epah, alpha
    real(8) :: theta_1, theta_2, delta, k_allowance
    character(len=1024) :: input_path, output_dir, seed_arg
    num_samplingnodes=1000
    allocate(norm_data(num_allnode, Time_window), &
             NeighborDis(num_allnode), Kernel_NeighborDis(num_allnode), &
             Nodesset(num_allnode, 2), Statistics(Time_window), &
             NeighborId(num_allnode), Kernel_NeighborId(num_allnode), &
             stat=allocation_status)
    if (allocation_status /= 0) then
        write(*, '(A)') 'Unable to allocate case-study arrays.'
        error stop 2
    end if
    input_path = 'cache/case-study/norm_data.bin'
    output_dir = 'results/case-study'
    seed_arg = ''
    call get_command_argument(1, input_path, status=arg_status)
    if (arg_status /= 0 .or. len_trim(input_path) == 0) &
        input_path = 'cache/case-study/norm_data.bin'
    call get_command_argument(2, output_dir, status=arg_status)
    if (arg_status /= 0 .or. len_trim(output_dir) == 0) &
        output_dir = 'results/case-study'
    call get_command_argument(3, seed_arg, status=arg_status)
    if (arg_status == 0 .and. len_trim(seed_arg) > 0) then
        read(seed_arg, *, iostat=arg_status) seed
        if (arg_status /= 0 .or. seed < 0) then
            write(*, '(A)') 'Seed must be a nonnegative integer.'
            error stop 2
        end if
        call RNSET(seed)
        write(*, '(A,I0)') 'IMSL random seed: ', seed
    end if
    call load_normalized_data(input_path, norm_data)

    pho=0.00001
    sigmaprio=0.1
    lambda=0.8
    Epah=0.02
    alpha=1.0
    CALL GeneratePoints(Nodesset)
    CALL EDistance(Nodesset, NeighborDis, NeighborId)
    CALL EpaKernelDistance(Nodesset, Epah, Kernel_NeighborDis, Kernel_NeighborId)
    call AS_Charting(pho, sigmaprio, alpha, lambda, num_samplingnodes, &
                     NeighborDis, NeighborId, Kernel_NeighborDis, &
                     Kernel_NeighborId, norm_data, Statistics)
    call write_statistics(output_dir, 'mM-KST', Statistics)

    theta_1=0.1
    theta_2=0.7
    Epah=0.02
    Limit = 35.6875

    call SASAM_Charting(Kernel_NeighborDis, Kernel_NeighborId, NeighborId, &
                        num_samplingnodes, Limit, theta_1, theta_2, &
                        norm_data, Statistics)
    call write_statistics(output_dir, 'SASAM', Statistics)

    delta=0.00001
    call Topr_Charting(num_samplingnodes, delta, norm_data, Statistics)
    call write_statistics(output_dir, 'TOPR', Statistics)

    call TSSPR_Charting(num_samplingnodes, norm_data, Statistics)
    call write_statistics(output_dir, 'TSSPR', Statistics, wide=.true.)
    delta=1.0/(dble(num_allnode)-dble(num_samplingnodes)/2)
    lambda=(1.0-dble(num_allnode-num_samplingnodes)*delta)/(1.0-(0.5**dble(num_samplingnodes)))
    k_allowance=dble(num_allnode-num_samplingnodes)*((delta**2.0)/gone-2.0*gone)+1.0
    call NAS_Charting(lambda,k_allowance,delta, num_samplingnodes, norm_data, Statistics)
    call write_statistics(output_dir, 'NAS', Statistics)

end program main
subroutine GeneratePoints(Nodesset)
    use constants
    implicit none
    real(8), intent(out) :: Nodesset(num_allnode, 2)  ! Output: 2D array storing (x, y) coordinates
    integer :: i, j
    real(8) :: x_start, y_start



    ! Compute the starting positions (center of the first grid cell)
    x_start = x_step / 2.0
    y_start = y_step / 2.0

    ! Generate a uniform grid of points within the unit square (0,0) to (1,1)
    do i = 1, num_y
        do j = 1, num_x
            Nodesset((i - 1) * num_x + j, 1) = x_start + (i - 1) * x_step  ! X-coordinate
            Nodesset((i - 1) * num_x + j, 2) = y_start + (j - 1) * y_step  ! Y-coordinate
        end do
    end do
end subroutine GeneratePoints

subroutine sortDis(Dis, ID, order)
    use constants
    implicit none

    real(8), intent(inout) :: Dis(num_allnode)
    integer, intent(out) :: ID(num_allnode)
    integer, intent(in) :: order
    real(8), allocatable :: work_dis(:)
    integer, allocatable :: work_id(:)
    integer :: i, j, left, middle, right, out, width
    logical :: take_left

    if (order /= 1 .and. order /= -1) then
        error stop 'sortDis: order must be 1 or -1'
    end if

    allocate(work_dis(num_allnode), work_id(num_allnode))
    ID = [(i, i = 1, num_allnode)]
    width = 1

    do while (width < num_allnode)
        left = 1
        do while (left <= num_allnode)
            middle = min(left + width - 1, num_allnode)
            right = min(left + 2 * width - 1, num_allnode)
            i = left
            j = middle + 1
            out = left

            do while (i <= middle .and. j <= right)
                if (order == 1) then
                    take_left = Dis(i) <= Dis(j)
                else
                    take_left = Dis(i) >= Dis(j)
                end if
                if (take_left) then
                    work_dis(out) = Dis(i)
                    work_id(out) = ID(i)
                    i = i + 1
                else
                    work_dis(out) = Dis(j)
                    work_id(out) = ID(j)
                    j = j + 1
                end if
                out = out + 1
            end do

            do while (i <= middle)
                work_dis(out) = Dis(i)
                work_id(out) = ID(i)
                i = i + 1
                out = out + 1
            end do
            do while (j <= right)
                work_dis(out) = Dis(j)
                work_id(out) = ID(j)
                j = j + 1
                out = out + 1
            end do
            left = right + 1
        end do
        Dis = work_dis
        ID = work_id
        width = 2 * width
    end do

    deallocate(work_dis, work_id)
end subroutine sortDis

! Swap subroutine for real numbers to avoid redundant code


subroutine EDistance(Nodesset, NeighborDis, NeighborId)
    use constants
    implicit none
    real(8), intent(in) :: Nodesset(num_allnode,2)  ! Input: Node coordinates
    real(8), intent(out) :: NeighborDis(num_allnode)  ! Output: Squared Euclidean distances
    integer, intent(out) :: NeighborId(num_allnode)  ! Output: Sorted indices by distance
    integer :: i, center_index
    real(8) :: Temone, Temtwo,x_len,y_len

    ! Compute the index of the center point in the grid
    center_index = (center_point - 1) * num_x + center_point
    x_len=x_step*num_x
    y_len=y_step*num_y
    ! Compute the squared Euclidean distance with periodic boundary conditions
    do i = 1, num_allnode
        Temone = abs(Nodesset(i,1) - Nodesset(center_index,1))
        Temtwo = abs(Nodesset(i,2) - Nodesset(center_index,2))

        ! Apply periodic boundary conditions
        Temone = min(Temone,  x_len - Temone)
        Temtwo = min(Temtwo,  y_len - Temtwo)

        ! Store squared distance
        NeighborDis(i) = (Temone**2) + (Temtwo**2)
        NeighborId(i) = i
    end do

    ! Sort distances in ascending order and update indices
    call sortDis(NeighborDis, NeighborId, 1)  ! 1 for ascending order
end subroutine EDistance


subroutine EpaKernelDistance(Nodesset, Epah, Kernel_NeighborDis, Kernel_NeighborId)
    use constants
    implicit none
    real(8), intent(in) :: Nodesset(num_allnode,2), Epah  ! Input: Node coordinates and Epanechnikov kernel bandwidth
    real(8), intent(out) :: Kernel_NeighborDis(num_allnode)  ! Output: Kernel distance values
    integer, intent(out) :: Kernel_NeighborId(num_allnode)  ! Output: Sorted indices
    integer :: i, center_index
    real(8) :: kernel_b, Temone, Temtwo

    ! Precompute kernel denominator
    kernel_b = Epah**2.0

    ! Compute the index of the center point in the grid
    center_index = (center_point - 1) * num_x + center_point

    ! Compute kernel-based distances
    do i = 1, num_allnode
        Temone = Nodesset(i,1) - Nodesset(center_index,1)
        Temtwo = Nodesset(i,2) - Nodesset(center_index,2)
        ! Compute normalized Epanechnikov kernel distance
        Kernel_NeighborDis(i) = ((Temone**2) + (Temtwo**2)) / kernel_b
        Kernel_NeighborDis(i) = min(Kernel_NeighborDis(i), 1.0)
        Kernel_NeighborId(i) = i
    end do
    ! Transform distances as per Epanechnikov kernel property
    Kernel_NeighborDis = 1.0 - Kernel_NeighborDis
    ! Sort distances in descending order and update indices
    call sortDis(Kernel_NeighborDis, Kernel_NeighborId, -1)  ! -1 for descending order
end subroutine EpaKernelDistance

subroutine AS_Allocation(pho, sigmaprio, alpha, lambda, num_samplingnodes, &
                         NeighborDis, NeighborId, Kernel_NeighborDis, Kernel_NeighborId, &
                         sampleReal, ProbabilityTerm, KernelVarianceTerm, DistanceTerm, &
                         stdTerm, charting_statistic)
    use constants
    use RNUND_INT
    USE, INTRINSIC :: IEEE_ARITHMETIC
    implicit none

    ! Input parameters
    real(8), intent(in) :: pho, sigmaprio, alpha, lambda
    real(8), intent(in) :: sampleReal(num_allnode)
    real(8), intent(in) :: NeighborDis(num_allnode), Kernel_NeighborDis(num_allnode)
    integer, intent(in) :: NeighborId(num_allnode), Kernel_NeighborId(num_allnode), num_samplingnodes

    ! Inout & Output parameters
    real(8), intent(inout) :: ProbabilityTerm(num_allnode), KernelVarianceTerm(num_allnode)
    real(8), intent(inout) :: DistanceTerm(num_allnode), stdTerm(num_allnode)
    real(8), intent(out) :: charting_statistic

    ! Local variables
    integer :: k, i, j, Nodes, max_index, one_vec(1), max_vector(num_allnode), near_index
    real(8) :: samplevalue, max_value, tem_max, tem_maxindex, tem_min, tem_minindex
    real(8) :: statisticmax, statisticmin
    real(8) :: TemVector(num_allnode), MeanTerm(num_allnode)
    integer :: x_index, y_index, x_delta, y_delta, near_id
    real(8) :: lambda_sq


    ! Precompute constants
    lambda_sq = lambda ** 2.0

    ! Update terms
    ProbabilityTerm = ProbabilityTerm * lambda
    KernelVarianceTerm = KernelVarianceTerm * lambda_sq
    DistanceTerm = DistanceTerm + pho

    ! Compute standard term and mean term
    stdTerm = ProbabilityTerm / ((KernelVarianceTerm+sigmaprio)**0.5)
    stdTerm = stdTerm ** 2.0
    stdTerm=MIN(stdTerm,500.0)
    MeanTerm = exp(stdTerm * alpha)

    ! Compute initial weighted distance
    TemVector = DistanceTerm * MeanTerm
    Nodes = 0

    ! Adaptive sampling loop
    do while (Nodes < num_samplingnodes)
        ! Find the maximum value in TemVector
        max_value = maxval(TemVector)
        k = 0
        do i = 1, num_allnode
            if (abs(TemVector(i) - max_value)<1.0E-10) then
                k = k + 1
                max_vector(k) = i
            end if
        end do

        ! Select one of the max candidates randomly if multiple exist
        if (k == 1) then
            max_index = max_vector(k)
        else
            call RNUND(k, one_vec)
            max_index = max_vector(one_vec(1))
        end if
        !print*,max_index
        ! Update sampled nodes count
        Nodes = Nodes + 1

        samplevalue =sampleReal(max_index)

        ! Compute the cordinate difference to the center_point
        y_delta = (max_index - 1) / num_x + 1 - center_point
        x_delta = mod(max_index - 1, num_x) + 1 - center_point

        ! Update DistanceTerm for neighbors
        near_index = 0
        do i = 1, num_allnode
            y_index = (NeighborId(i) - 1) / num_x + 1 + y_delta
            x_index = mod(NeighborId(i) - 1, num_x) + 1 + x_delta

            ! Apply periodic boundary conditions
            if (x_index > num_x) x_index = x_index - num_x
            if (x_index < 1) x_index = num_x + x_index
            if (y_index > num_y) y_index = y_index - num_y
            if (y_index < 1) y_index = num_y + y_index

            near_id = (y_index - 1) * num_x + x_index

            ! Update DistanceTerm and TemVector if necessary
            if (NeighborDis(i) <= DistanceTerm(near_id)) then
                DistanceTerm(near_id) = NeighborDis(i)
                TemVector(near_id) = MeanTerm(near_id) * DistanceTerm(near_id)
                near_index = 0
            else
                near_index = near_index + 1
            end if

            ! Early exit if threshold exceeded
            if (near_index > near_Threshold) exit
        end do

        ! Update KernelVarianceTerm and ProbabilityTerm for neighbors
        do i = 1, num_allnode
            y_index = (Kernel_NeighborId(i) - 1) / num_x + 1 + y_delta
            x_index = mod(Kernel_NeighborId(i) - 1, num_x) + 1 + x_delta
            if (x_index <= 0 .or. x_index > num_x .or. y_index <= 0 .or. y_index > num_y) cycle
            near_id = (y_index - 1) * num_x + x_index
            if (Kernel_NeighborDis(i) > Kernel_threshold) then
                KernelVarianceTerm(near_id) = KernelVarianceTerm(near_id) + (Kernel_NeighborDis(i) ** 2.0)
                ProbabilityTerm(near_id) = ProbabilityTerm(near_id) + Kernel_NeighborDis(i) * samplevalue
            else
                exit
            end if
        end do
    end do
    !print*,"++++++++++++++++++++++++++"
    ! Compute standard term
    stdTerm = ProbabilityTerm / (sqrt(KernelVarianceTerm) + 1.0E-10)

    ! Compute top-r statistical values
    TemVector = stdTerm
    statisticmax = 0.0
    statisticmin = 0.0
    do i = 1, Top_r
        tem_max = -1.0E10
        tem_min = 1.0E10
        do j = 1, num_allnode
            if (TemVector(j) > tem_max) then
                tem_max = TemVector(j)
                tem_maxindex = j
            end if
            if (TemVector(j) < tem_min) then
                tem_min = TemVector(j)
                tem_minindex = j
            end if
        end do
        statisticmax = statisticmax + tem_max
        statisticmin = statisticmin + tem_min
        TemVector(tem_maxindex) = 0.0
        TemVector(tem_minindex) = 0.0
    end do

    ! Final charting statistic
    charting_statistic = max(statisticmax, abs(statisticmin))
end subroutine AS_Allocation

subroutine AS_Charting(pho, sigmaprio, alpha, lambda, num_samplingnodes, &
                       NeighborDis, NeighborId, Kernel_NeighborDis, &
                       Kernel_NeighborId, norm_data, Statistics)
    use constants
    implicit none

    ! Input parameters
    real(8), intent(in) :: pho             ! Parameter for allocation
    real(8), intent(in) :: sigmaprio       ! Prior standard deviation
    real(8), intent(in) :: alpha           ! Learning rate parameter
    real(8), intent(in) :: lambda          ! Regularization parameter
    real(8), intent(in) :: NeighborDis(num_allnode)         ! Distance to neighbors
    real(8), intent(in) :: Kernel_NeighborDis(num_allnode)  ! Distance to kernel neighbors
    integer, intent(in) :: NeighborId(num_allnode)          ! Neighbor node IDs
    integer, intent(in) :: Kernel_NeighborId(num_allnode)   ! Kernel neighbor node IDs
    integer, intent(in) :: num_samplingnodes                ! Number of sampled nodes
    real(8), intent(in) :: norm_data(num_allnode,Time_window)

    ! Output parameter
    real(8), intent(out) :: Statistics(Time_window)

    ! Local variables
    integer :: j                   ! Loop index
    real(8) :: charting_statistic         ! Monitored statistic
    real(8) :: ProbabilityTerm(num_allnode)  ! Probability term for allocation
    real(8) :: KernelVarianceTerm(num_allnode)  ! Kernel variance term
    real(8) :: DistanceTerm(num_allnode)  ! Distance-related term
    real(8) :: stdTerm(num_allnode)       ! Standard deviation term

    ProbabilityTerm = 0.0
    DistanceTerm = 100.0
    KernelVarianceTerm = 0.0
    stdTerm = 0.0
    do j = 1, Time_window
        call AS_Allocation(pho, sigmaprio, alpha, lambda, &
                           num_samplingnodes, NeighborDis, NeighborId, &
                           Kernel_NeighborDis, Kernel_NeighborId, &
                           norm_data(:,j), ProbabilityTerm, KernelVarianceTerm, &
                           DistanceTerm, stdTerm, charting_statistic)
        Statistics(j)=charting_statistic
        print*,j
        print*,charting_statistic
        print*,"----------------------"
    end do

    return
end subroutine AS_Charting



subroutine SASAM_Allocation(Kernel_NeighborDis, Kernel_NeighborId, NeighborId, &
                            num_samplingnodes, limit, theta_1, theta_2, &
                            sampleReal, cusumdouble_statistic, cusum_statistic, charting_statistic)
    USE constants
    use RNUND_INT
    implicit none

    ! Input parameters
    real(8), intent(in) :: sampleReal(num_allnode), Kernel_NeighborDis(num_allnode), limit, theta_1, theta_2
    integer, intent(in) :: Kernel_NeighborId(num_allnode), NeighborId(num_allnode), num_samplingnodes

    ! In-out parameters
    real(8), intent(inout) :: cusumdouble_statistic(num_allnode,2), cusum_statistic(num_allnode)
    real(8), intent(inout) :: charting_statistic

    ! Local variables
    integer :: i, j, k, Nodes, one_vec(1), max_index, D_nodes, W_Nodes
    real(8) :: samplevalue, max_value
    real(8) :: Tem_vector(num_allnode)
    integer :: sampling_index(num_allnode)
    integer :: x_index, y_index, x_delta, y_delta, near_id
    integer :: xx_index, yy_index, xx_delta, yy_delta, nnear_id

    ! Initialize arrays
    sampling_index = 0
    max_value = charting_statistic - limit * theta_1
    max_value = max(max_value, 0.0)

    ! Compute D_nodes and W_Nodes
    D_nodes = max_value * num_samplingnodes * theta_2 / (limit * (1 - theta_1))
    D_nodes = min(D_nodes, num_samplingnodes)
    W_Nodes = num_samplingnodes - D_nodes

    ! Select top D_nodes based on cusum_statistic
    if (D_nodes > 0) then
        max_value = -1
        do i = 1, num_allnode
            if (max_value < cusum_statistic(i)) then
                max_value = cusum_statistic(i)
                max_index = i
            end if
        end do

        y_delta = (max_index - 1) / num_x + 1 - center_point
        x_delta = MOD(max_index - 1, num_x) + 1 - center_point

        do i = 1, D_nodes
            y_index = (NeighborId(i) - 1) / num_x + 1 + y_delta
            x_index = MOD(NeighborId(i) - 1, num_x) + 1 + x_delta

            if (x_index > num_x .OR. x_index <= 0 .OR. y_index > num_y .OR. y_index <= 0) cycle
            near_id = (y_index - 1) * num_x + x_index

            xx_delta = x_index - center_point
            yy_delta = y_index - center_point
            samplevalue =sampleReal(near_id)
            sampling_index(near_id) = 1
            do j = 1, num_allnode
                if (Kernel_NeighborDis(j) <= Kernel_threshold) exit

                xx_index = (Kernel_NeighborId(j) - 1) / num_x + 1 + xx_delta
                yy_index = MOD(Kernel_NeighborId(j) - 1, num_x) + 1 + yy_delta

                if (xx_index > num_x .OR. xx_index <= 0 .OR. yy_index > num_y .OR. yy_index <= 0) cycle
                nnear_id = (yy_index - 1) * num_x + xx_index
                cusumdouble_statistic(nnear_id,1) = cusumdouble_statistic(nnear_id,1) + &
                    Kernel_NeighborDis(j) * (miu_min * samplevalue - (miu_min**2.0) / 2.0)
                cusumdouble_statistic(nnear_id,2) = cusumdouble_statistic(nnear_id,2) + &
                    Kernel_NeighborDis(j) * (-miu_min * samplevalue - (miu_min**2.0) / 2.0)
            end do
        end do
    endif

    ! Select remaining W_Nodes randomly
    Nodes = D_nodes
    do while (Nodes < num_samplingnodes)
        CALL RNUND(num_allnode, one_vec)
        k = one_vec(1)
        if (sampling_index(k) == 1) cycle

        Nodes = Nodes + 1
        samplevalue =sampleReal(k)

        sampling_index(k) = 1

        xx_delta = (k - 1) / num_x + 1 - center_point
        yy_delta = MOD(k - 1, num_x) + 1 - center_point

        do j = 1, num_allnode
            if (Kernel_NeighborDis(j) <= Kernel_threshold) exit

            xx_index = (Kernel_NeighborId(j) - 1) / num_x + 1 + xx_delta
            yy_index = MOD(Kernel_NeighborId(j) - 1, num_x) + 1 + yy_delta

            if (xx_index > num_x .OR. xx_index <= 0 .OR. yy_index > num_y .OR. yy_index <= 0) cycle
            nnear_id = (xx_index - 1) * num_x + yy_index
            cusumdouble_statistic(nnear_id,1) = cusumdouble_statistic(nnear_id,1) + &
                Kernel_NeighborDis(j) * (miu_min * samplevalue - (miu_min**2.0) / 2.0)
            cusumdouble_statistic(nnear_id,2) = cusumdouble_statistic(nnear_id,2) + &
                Kernel_NeighborDis(j) * (-miu_min * samplevalue - (miu_min**2.0) / 2.0)
        end do
    end do

    ! Update CUSUM statistics
    do i = 1, num_allnode
        cusumdouble_statistic(i,1) = max(cusumdouble_statistic(i,1), 0.0)
        cusumdouble_statistic(i,2) = max(cusumdouble_statistic(i,2), 0.0)
        cusum_statistic(i) = max(cusumdouble_statistic(i,1), cusumdouble_statistic(i,2))
    end do

    ! Compute charting statistic based on Top_r highest values
    Tem_vector = cusum_statistic
    charting_statistic = 0.0
    do i = 1, Top_r
        max_value = 0
        do j = 1, num_allnode
            if (max_value <= Tem_vector(j)) then
                max_value = Tem_vector(j)
                max_index = j
            end if
        end do
        charting_statistic = charting_statistic + max_value
        Tem_vector(max_index) = -1
    end do

end subroutine SASAM_Allocation

subroutine SASAM_Charting(Kernel_NeighborDis, Kernel_NeighborId, NeighborId, &
                          num_samplingnodes, Limit, theta_1, theta_2, norm_data, Statistics)
    use constants
    implicit none

    ! Input parameters
    real(8), intent(in) :: theta_1        ! Parameter for SASAM algorithm
    real(8), intent(in) :: theta_2        ! Parameter for SASAM algorithm
    real(8), intent(in) :: Limit          ! Control limit threshold
    real(8), intent(in) :: Kernel_NeighborDis(num_allnode)  ! Distance to kernel neighbors
    integer, intent(in) :: Kernel_NeighborId(num_allnode)   ! Kernel neighbor node IDs
    integer, intent(in) :: NeighborId(num_allnode)          ! Neighbor node IDs
    integer, intent(in) :: num_samplingnodes                ! Number of sampled nodes
    real(8), intent(in) :: norm_data(num_allnode,Time_window)

    ! Output parameter
    real(8), intent(out) :: Statistics(Time_window)
    ! Local variables
    integer :: j
    real(8) :: charting_statistic          ! Monitoring statistic
    real(8) :: cusum_statistic(num_allnode) ! CUSUM statistic for monitoring
    real(8) :: cusumdouble_statistic(num_allnode, 2) ! Double CUSUM statistic

    ! Initialize variables
    cusum_statistic = 0.0
    cusumdouble_statistic = 0.0

    do j = 1, Time_window
        call SASAM_Allocation(Kernel_NeighborDis, Kernel_NeighborId, NeighborId, &
                              num_samplingnodes, Limit, theta_1, theta_2, &
                              norm_data(:,j), cusumdouble_statistic, &
                              cusum_statistic, charting_statistic)
        Statistics(j)=charting_statistic
        print*,j
        print*,charting_statistic
        print*,"----------------------"
    end do
    return
end subroutine SASAM_Charting

subroutine TopR_Allocation(sampleReal, num_samplingnodes, delta, cusumdouble_statistic, cusum_statistic, charting_statistic)
    USE constants
    use RNUND_INT
    implicit none

    ! Input parameters
    real(8), intent(in) :: sampleReal(num_allnode)
    real(8), intent(in) :: delta
    integer, intent(in) :: num_samplingnodes

    ! Inout parameters
    real(8), intent(inout) :: cusumdouble_statistic(num_allnode,2)
    real(8), intent(inout) :: cusum_statistic(num_allnode)

    ! Output parameter
    real(8), intent(out) :: charting_statistic

    ! Local variables
    real(8) :: max_value, samplevalue
    integer :: i, j, k, nodes, max_index, remaining_nodes
    real(8) :: temp_vector(num_allnode)
    integer :: max_vector(num_allnode)
    integer :: one_vec(1)

    ! Copy the CUSUM statistic into a temporary array
    temp_vector = cusum_statistic
    nodes = 0

    ! Select the top 'num_samplingnodes' nodes
    do while (nodes < num_samplingnodes)
        ! Find the maximum value in the temporary vector
        max_value = maxval(temp_vector)
        k = 0

        ! Identify all nodes with the maximum value
        do i = 1, num_allnode
            if (abs(temp_vector(i) - max_value) < 1.0d-12) then
                k = k + 1
                max_vector(k) = i
            end if
        end do

        remaining_nodes = num_samplingnodes - nodes

        if (k <= remaining_nodes) then
            ! If all identified nodes can be selected
            do j = 1, k
                max_index = max_vector(j)
                temp_vector(max_index) = -1.0
                nodes = nodes + 1

            samplevalue =sampleReal(max_index)
            cusumdouble_statistic(max_index,1) = max(cusumdouble_statistic(max_index,1) + (miu_min * samplevalue - cuparame), 0.0)
            cusumdouble_statistic(max_index,2) = max(cusumdouble_statistic(max_index,2) + (-miu_min * samplevalue - cuparame), 0.0)
            cusum_statistic(max_index) = max(cusumdouble_statistic(max_index,1), cusumdouble_statistic(max_index,2))
            end do
        else
            ! Randomly select 'remaining_nodes' from 'k' candidates
            j = 1
            do while (j <= remaining_nodes)
                call RNUND(k, one_vec)
                if (max_vector(one_vec(1)) == 0) then
                    cycle
                end if
                max_index = max_vector(one_vec(1))
                max_vector(one_vec(1)) = 0
                temp_vector(max_index) = -1.0
                nodes = nodes + 1
                samplevalue =sampleReal(max_index)
                cusumdouble_statistic(max_index,1) = max( &
                    cusumdouble_statistic(max_index,1) + (miu_min * samplevalue - cuparame), 0.0)
                cusumdouble_statistic(max_index,2) = max( &
                    cusumdouble_statistic(max_index,2) + (-miu_min * samplevalue - cuparame), 0.0)
                cusum_statistic(max_index) = max(cusumdouble_statistic(max_index,1), cusumdouble_statistic(max_index,2))
                j = j + 1
            end do
        end if
    end do

    ! Update CUSUM statistics for all unselected nodes
    do j = 1, num_allnode
        if ((temp_vector(j)+1.0)<=1.0d-12) then
            cycle
        end if
        cusumdouble_statistic(j,1) = cusumdouble_statistic(j,1) + delta
        cusumdouble_statistic(j,2) = cusumdouble_statistic(j,2) + delta
        cusum_statistic(j) = max(cusumdouble_statistic(j,1), cusumdouble_statistic(j,2))
    end do
    ! Compute the Top-R charting statistic
    temp_vector = cusum_statistic
    charting_statistic = 0.0

    do i = 1, Top_r
        max_value = 0.0
        do j = 1, num_allnode
            if (max_value <= temp_vector(j)) then
                max_value = temp_vector(j)
                max_index = j
            end if
        end do
        charting_statistic = charting_statistic + max_value
        temp_vector(max_index) = -1.0
    end do
end subroutine TopR_Allocation

subroutine Topr_Charting(num_samplingnodes, delta, norm_data, Statistics)
    USE constants
    implicit none

    ! Input variables
    real(8), intent(in) :: delta
    integer, intent(in) :: num_samplingnodes
    real(8), intent(in) :: norm_data(num_allnode,Time_window)

    ! Output parameter
    real(8), intent(out) :: Statistics(Time_window)

    integer :: j
    real(8) :: charting_statistic
    real(8) :: cusumdouble_statistic(num_allnode,2)
    real(8) :: cusum_statistic(num_allnode)


    cusum_statistic = 0.0
    cusumdouble_statistic = 0.0
    do j = 1, Time_window
        CALL TopR_Allocation(norm_data(:,j), num_samplingnodes, delta, cusumdouble_statistic, cusum_statistic, charting_statistic)
        Statistics(j)=charting_statistic
        print*,j
        print*,charting_statistic
        print*,"----------------------"
    end do
    return

    end subroutine Topr_Charting

   subroutine TSSPR_Allocation(num_samplingnodes, sampleReal, L_statistic, R_statistic, charting_statistic)
    USE constants
    use RNUN_INT
    implicit none

    ! Input & Output Variables
    integer, intent(in)  :: num_samplingnodes  ! Number of sampling nodes
    real(8), intent(in)   :: sampleReal(num_allnode)  ! Shift pattern for each node
    real(8), intent(inout):: L_statistic(num_allnode,2), R_statistic(num_allnode,2)  ! Left & Right statistics
    real(8), intent(out)  :: charting_statistic  ! Final charting statistic

    ! Local Variables
    real(8) :: Sample_statistic(num_allnode)  ! Random uniform samples
    real(8) :: TS_statistic(num_allnode)      ! Test statistic for each node
    real(8) :: samplevalue, max_value, Temvalue
    real(8) :: Tem_vector(num_allnode)        ! Temporary vector for ranking
    integer :: i, j, Nodes, max_index
    integer :: max_locvec(1)

    ! Generate uniform random numbers
    CALL RNUN(Sample_statistic, num_allnode)

    ! Compute TS_statistic for each node
    do i=1,num_allnode
        TS_statistic(i) = max(R_statistic(i,1) + L_statistic(i,1) * Sample_statistic(i), &
                              R_statistic(i,2) + L_statistic(i,2) * Sample_statistic(i))
    end do

    ! Initialize R_statistic (increases by 1)
    R_statistic = R_statistic + 1.0

    Nodes = 0
    do while (Nodes < num_samplingnodes)
        ! Find the node with the maximum TS_statistic
        max_locvec = MAXLOC(TS_statistic)
        max_index = max_locvec(1)
        max_value = TS_statistic(max_index)

        ! Mark the selected node as used
        TS_statistic(max_index) = -1.0

        samplevalue =sampleReal(max_index)


        ! Compute transformation factor
        Temvalue = miu_min * samplevalue - cuparame
        Temvalue=MIN(Temvalue,100.0)
        Temvalue = Temvalue
        R_statistic(max_index,1) = R_statistic(max_index,1) * Temvalue
        L_statistic(max_index,1) = L_statistic(max_index,1) * Temvalue

        Temvalue = -miu_min * samplevalue - cuparame
        Temvalue=MIN(Temvalue,100.0)
        Temvalue = Temvalue
        R_statistic(max_index,2) = R_statistic(max_index,2) * Temvalue
        L_statistic(max_index,2) = L_statistic(max_index,2) * Temvalue

        Nodes = Nodes + 1
    end do

    ! Compute max R_statistic across both columns
    do i=1,num_allnode
        Tem_vector(i) = max(R_statistic(i,1), R_statistic(i,2))
    end do

    ! Compute the final charting statistic (sum of top `Top_r` values)
    charting_statistic = 0.0
    do i = 1, Top_r
        max_locvec = MAXLOC(Tem_vector)
        max_index = max_locvec(1)
        max_value = Tem_vector(max_index)

        charting_statistic = charting_statistic + max_value
        Tem_vector(max_index) = -1.0  ! Mark as used
    end do

    return
end subroutine TSSPR_Allocation

subroutine TSSPR_Charting(num_samplingnodes, norm_data, Statistics)
    use constants
    implicit none

    ! Input parameters
    integer, intent(in) :: num_samplingnodes                ! Number of sampled nodes
    real(8), intent(in) :: norm_data(num_allnode,Time_window)

    ! Output parameter
    real(8), intent(out) :: Statistics(Time_window)

    ! Local variables
    integer :: j               ! System clock-based seed
    real(8) :: charting_statistic          ! Monitoring statistic
    real(8):: L_statistic(num_allnode,2), R_statistic(num_allnode,2)

    L_statistic = 1.0
    R_statistic = 0.0
    do j = 1, Time_window
        CALL TSSPR_Allocation(num_samplingnodes, norm_data(:,j), L_statistic, R_statistic, charting_statistic)
        Statistics(j)=charting_statistic
        print*,j
        print*,charting_statistic
        print*,"----------------------"
    end do

    return
end subroutine TSSPR_Charting

subroutine NAS_allocation(num_samplingnodes, lambda, k_allowance, delta, &
                          sampleReal, Smin_one, Smin_two, Smax_one, Smax_two, charting_statistic)
    use constants
    USE RNSRI_INT
    implicit none

    ! Input Parameters
    integer, intent(in) :: num_samplingnodes       ! Number of sampling nodes
    real(8), intent(in) :: lambda,k_allowance,delta     ! Parameters for updating statistics
    real(8), intent(in) :: sampleReal(num_allnode)  ! Shift pattern for nodes

    ! Input/Output Parameters
    real(8), intent(inout) :: Smin_one(num_allnode), Smin_two(num_allnode)
    real(8), intent(inout) :: Smax_one(num_allnode), Smax_two(num_allnode)

    ! Output Parameter
    real(8), intent(out) :: charting_statistic

    ! Local Variables
    real(8) :: S_one(num_allnode)  ! Temporary storage for maximum statistic
    real(8) :: sampling_value(num_samplingnodes)  ! Stores sampled values
    real(8) :: rankindex_s(num_allnode)  ! Rank index for sorting
    real(8) :: min_value, max_value, test_min, test_max, test
    real(8) :: C_index, S_max
    integer :: i, j, k, d, min_index, max_index, tem_index
    integer :: loc_index(num_allnode)  ! Location indices
    integer :: tem(num_samplingnodes)  ! Temporary storage
    integer :: sampling_loc(num_samplingnodes)  ! Sampled node locations

    ! Initialize S_one with the maximum values from Smin_one and Smax_one
    do i = 1, num_allnode
        S_one(i) = max(Smin_one(i), Smax_one(i))
    end do

    ! Sampling procedure
    k = 1
    S_max = 0.0
    do while (k <= num_samplingnodes)
        d = 0
        do i = 1, num_allnode
            if (S_one(i) > S_max) then
                S_max = S_one(i)
                d = 1
                loc_index(d) = i
            else if (S_one(i) == S_max) then
                d = d + 1
                loc_index(d) = i
            end if
        end do

        ! Select a subset of d nodes for sampling
        tem_index = min(d, num_samplingnodes - k + 1)
        CALL RNSRI (d, tem(1:tem_index))

        do i = k, k + tem_index - 1
            sampling_loc(i) = loc_index(tem(i - k + 1))
            sampling_value(i) =sampleReal(sampling_loc(i))
            S_one(sampling_loc(i)) = -1.0  ! Mark the sampled node as visited
        end do

        S_max = 0.0
        k = k + tem_index
    end do
    !print*,"------------"
    !print*,sampling_loc
    ! Compute minimum and maximum charting statistics
    rankindex_s = delta

    ! Update Smin statistics
    do i = 1, num_samplingnodes
        rankindex_s(sampling_loc(i)) = 0
    end do
    min_index = minloc(sampling_value, 1)

    Smin_one = Smin_one + rankindex_s
    Smin_two = Smin_two + gone

    Smin_one(sampling_loc(min_index)) = Smin_one(sampling_loc(min_index)) + lambda

    C_index = 0.0
    do i = 1, num_allnode
        C_index = C_index + ((Smin_one(i) - Smin_two(i))**2.0) / Smin_two(i)
    end do

    if (C_index <= k_allowance) then
        Smin_one = 0.0
        Smin_two = 0.0
        test_min = 0.0
    else
        Smin_one = Smin_one * (C_index - k_allowance) / C_index
        Smin_two = Smin_two * (C_index - k_allowance) / C_index
        test_min = C_index - k_allowance
    end if

    ! Update Smax statistics
    max_index = maxloc(sampling_value, 1)

    Smax_one = Smax_one + rankindex_s
    Smax_two = Smax_two + gone


    Smax_one(sampling_loc(max_index)) = Smax_one(sampling_loc(max_index)) + lambda

    C_index = 0.0
    do i = 1, num_allnode
        C_index = C_index + ((Smax_one(i) - Smax_two(i))**2.0) / Smax_two(i)
    end do

    if (C_index <= k_allowance) then
        Smax_one = 0.0
        Smax_two = 0.0
        test_max = 0.0
    else
        Smax_one = Smax_one * (C_index - k_allowance) / C_index
        Smax_two = Smax_two * (C_index - k_allowance) / C_index
        test_max = C_index - k_allowance
    end if

    ! Final charting statistic
    !print*,test_min
    !print*,test_max
    !print*,"-------"
    !pause
    charting_statistic = max(test_min, test_max)
    return
end subroutine NAS_allocation

subroutine NAS_Charting(lambda,k_allowance,delta, num_samplingnodes, norm_data,Statistics)
    use constants
    implicit none

    ! Input parameters
    real(8), intent(in) :: lambda,k_allowance,delta
    integer, intent(in) :: num_samplingnodes                ! Number of sampled nodes
    real(8), intent(in) :: norm_data(num_allnode,Time_window)

    ! Output parameter
    real(8), intent(out) :: Statistics(Time_window)

    ! Local variables
    integer :: j
    real(8) :: charting_statistic         ! Monitored statistic
    real(8):: Smin_one(num_allnode),Smin_two(num_allnode),Smax_one(num_allnode),Smax_two(num_allnode)


    Smin_one=0.0
    Smin_two=0.0
    Smax_one=0.0
    Smax_two=0.0
    do j = 1, Time_window
        call NAS_allocation(num_samplingnodes, lambda, k_allowance, delta, &
                            norm_data(:,j), Smin_one, Smin_two, &
                            Smax_one, Smax_two, charting_statistic)
        Statistics(j)=charting_statistic
        print*,j
        print*,charting_statistic
        print*,"----------------------"
    end do
    return
end subroutine NAS_Charting

