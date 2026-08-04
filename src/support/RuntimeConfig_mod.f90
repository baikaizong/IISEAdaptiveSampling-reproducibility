module RuntimeConfig_mod
    use mpi
    implicit none
    private

    character(len=1024), save, public :: results_root = "results"
    character(len=1024), save, public :: cache_root = "cache"
    integer, save :: fixed_seed = -1

    public :: init_runtime_config
    public :: initialize_evaluation
    public :: finalize_evaluation
    public :: cache_path
    public :: ensure_directory
    public :: get_clock_seed
    public :: seed_intrinsic_random

contains

subroutine init_runtime_config(first_path_arg, create_directories)
    integer, intent(in), optional :: first_path_arg
    logical, intent(in), optional :: create_directories
    integer :: arg0, argc, ios, env_status
    character(len=1024) :: value
    logical :: make_directories

    arg0 = 1
    if (present(first_path_arg)) arg0 = first_path_arg
    make_directories = .true.
    if (present(create_directories)) make_directories = create_directories
    argc = command_argument_count()

    if (argc >= arg0) then
        call get_command_argument(arg0, value)
        if (len_trim(value) > 0) results_root = trim(value)
    end if
    if (argc >= arg0 + 1) then
        call get_command_argument(arg0 + 1, value)
        if (len_trim(value) > 0) cache_root = trim(value)
    end if

    value = ""
    call get_environment_variable("MMKSTD_SEED", value=value, status=env_status)
    if (env_status == 0 .and. len_trim(value) > 0) then
        read(value, *, iostat=ios) fixed_seed
        if (ios /= 0) fixed_seed = -1
    end if
    if (argc >= arg0 + 2) then
        call get_command_argument(arg0 + 2, value)
        if (len_trim(value) > 0) then
            read(value, *, iostat=ios) fixed_seed
            if (ios /= 0) fixed_seed = -1
        end if
    end if

    if (make_directories) then
        call ensure_directory(trim(results_root))
        call ensure_directory(trim(cache_root))
        call ensure_directory(cache_path("Prepared_settings"))
        call ensure_directory(cache_path("Prepared_settings-1"))
        call ensure_directory(cache_path("Prepared_settings-2"))
        call ensure_directory(cache_path("Prepared_settings-3"))
        call ensure_directory(cache_path("Prepared_settings-4"))
    end if
end subroutine init_runtime_config

subroutine initialize_evaluation(default_experiment, experiment, rank)
    character(len=*), intent(in) :: default_experiment
    character(len=*), intent(out) :: experiment
    integer, intent(out) :: rank
    integer :: ierr, process_count

    experiment = default_experiment
    call get_command_argument(1, experiment)
    if (len_trim(experiment) == 0) experiment = default_experiment
    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, process_count, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    if (process_count < 2) then
        if (rank == 0) write(*, '(A)') "Error: at least two MPI processes are required."
        call MPI_Finalize(ierr)
        stop 2
    end if

    call init_runtime_config(2, rank == 0)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
end subroutine initialize_evaluation

subroutine finalize_evaluation()
    integer :: ierr

    call MPI_Finalize(ierr)
end subroutine finalize_evaluation

function cache_path(name) result(path)
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: path

    path = trim(cache_root) // "/" // trim(name)
end function cache_path

subroutine ensure_directory(path)
    character(len=*), intent(in) :: path
    character(len=2300) :: command
    character(len=64) :: operating_system
    integer :: env_status, exit_status

    if (len_trim(path) == 0) return

    operating_system = ""
    call get_environment_variable("OS", value=operating_system, status=env_status)
    if (env_status == 0 .and. index(operating_system, "Windows_NT") > 0) then
        command = 'if not exist "' // trim(path) // '" mkdir "' // trim(path) // '"'
    else
        command = 'mkdir -p "' // trim(path) // '"'
    end if

    call execute_command_line(trim(command), wait=.true., exitstat=exit_status)
    if (exit_status /= 0) then
        write(*, '(A)') "Warning: could not create directory: " // trim(path)
    end if
end subroutine ensure_directory

subroutine get_clock_seed(seed)
    integer, intent(out) :: seed

    if (fixed_seed >= 0) then
        seed = fixed_seed
    else
        call system_clock(count=seed)
    end if
end subroutine get_clock_seed

subroutine seed_intrinsic_random(stream_id)
    integer, intent(in), optional :: stream_id
    integer :: i, n, base_seed, stream
    integer, allocatable :: seed_values(:)

    stream = 0
    if (present(stream_id)) stream = stream_id
    call get_clock_seed(base_seed)
    call random_seed(size=n)
    allocate(seed_values(n))
    do i = 1, n
        seed_values(i) = modulo(base_seed + 104729 * stream + 8191 * i, 2147483646) + 1
    end do
    call random_seed(put=seed_values)
    deallocate(seed_values)
end subroutine seed_intrinsic_random

end module RuntimeConfig_mod
