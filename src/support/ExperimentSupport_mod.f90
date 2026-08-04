module ExperimentSupport_mod
    use mpi
    use GlobalSettings_mod, only: dp, ARL_Stats_type
    implicit none
    private

    public :: Init_Experiment_Common
    public :: Print_Banner
    public :: Write_OCARL_Header
    public :: Write_OCARL_Row

contains

subroutine Init_Experiment_Common(rank, process_count)
    integer, intent(out) :: rank, process_count
    integer :: ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, process_count, ierr)
end subroutine Init_Experiment_Common

subroutine Print_Banner(title, output_directory)
    character(len=*), intent(in) :: title, output_directory

    write(*, '(A)') "=========================================================="
    write(*, '(A)') " STARTING EXPERIMENT: " // trim(title)
    write(*, '(A)') " OUTPUT DIRECTORY: " // trim(output_directory)
    write(*, '(A)') "=========================================================="
end subroutine Print_Banner

subroutine Write_OCARL_Header(filename, prefix_columns)
    character(len=*), intent(in) :: filename, prefix_columns

    open(unit=30, file=filename, status='unknown', position='rewind')
    write(30, '(A)') trim(prefix_columns) // " &
        & OC_Mean OC_Std OC_Q1 OC_Med OC_Q3 OC_Max OC_Min &
        & Exp_Mean Exp_Std Exp_Q1 Exp_Med Exp_Q3 Exp_Max Exp_Min &
        & Num_Mean Num_Std Num_Q1 Num_Med Num_Q3 Num_Max Num_Min"
    close(30)
end subroutine Write_OCARL_Header

subroutine Write_OCARL_Row(filename, prefix_values, oc, exploration, exploitation)
    character(len=*), intent(in) :: filename
    real(dp), intent(in) :: prefix_values(:)
    type(ARL_Stats_type), intent(in) :: oc, exploration, exploitation
    integer :: i

    open(unit=30, file=filename, status='old', position='append')
    do i = 1, size(prefix_values)
        write(30, '(F8.4, 1X)', advance='no') prefix_values(i)
    end do
    write(30, '(3(1X, F12.4, F10.4, 5F10.2))') &
        oc%mean, oc%std_err, oc%q1, oc%median, oc%q3, oc%max_val, oc%min_val, &
        exploration%mean, exploration%std_err, exploration%q1, exploration%median, &
        exploration%q3, exploration%max_val, exploration%min_val, &
        exploitation%mean, exploitation%std_err, exploitation%q1, exploitation%median, &
        exploitation%q3, exploitation%max_val, exploitation%min_val
    close(30)
end subroutine Write_OCARL_Row

end module ExperimentSupport_mod
