program main_driver
    use Experiment_Runner_mod
    use RuntimeConfig_mod, only: initialize_evaluation, finalize_evaluation, results_root
    implicit none

    integer :: my_rank
    character(len=64) :: experiment

    call initialize_evaluation("bspline", experiment, my_rank)

    if (trim(experiment) == "bspline") then
        call Run_TSBSS_Bspline_Test(results_root)
    else if (my_rank == 0) then
        write(*, '(A)') "Unknown TSBSS experiment: " // trim(experiment)
    end if
    call finalize_evaluation()
end program main_driver
