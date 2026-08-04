program main_driver
    use Experiment_Runner_mod
    use RuntimeConfig_mod, only: initialize_evaluation, finalize_evaluation, results_root
    implicit none

    integer :: my_rank
    character(len=64) :: experiment

    call initialize_evaluation("kernel-circle", experiment, my_rank)

    select case (trim(experiment))
    case ("bspline")
        call Run_TOPR_Bspline_Test(results_root)
    case ("kernel-circle")
        call Run_TOPR_Kernel_CircleConv_Test(results_root)
    case default
        if (my_rank == 0) write(*, '(A)') "Unknown TOPR experiment: " // trim(experiment)
    end select
    call finalize_evaluation()
end program main_driver
