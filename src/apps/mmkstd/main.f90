program main_driver
    use Experiment_Runner_mod
    use RuntimeConfig_mod, only: initialize_evaluation, finalize_evaluation, results_root
    implicit none

    integer :: my_rank
    character(len=64) :: experiment

    call initialize_evaluation("circle-comparison", experiment, my_rank)

    select case (trim(experiment))
    case ("bspline")
        call Run_mMSTD_Bspline_Test(results_root)
    case ("kernel-circle")
        call Run_mMSTD_Kernel_CircleConv_Test(results_root)
    case ("circle-comparison")
        call Run_mMSTD_Circle_comparison(results_root)
    case ("circle")
        call Run_mMSTD_Circle_Test(results_root)
    case ("spatial-comparison")
        call Run_mMSTD_Spatial_comparison(results_root)
    case ("st-comparison")
        call Run_mMSTD_ST_comparison(results_root)
    case ("ellipse")
        call Run_mMSTD_Ellipse_Test(results_root)
    case ("crescent")
        call Run_mMSTD_Crescent_Test(results_root)
    case ("qnum")
        call Run_mMSTD_Circle_QnumTest(results_root)
    case default
        if (my_rank == 0) write(*, '(A)') "Unknown mMKSTD experiment: " // trim(experiment)
    end select

    call finalize_evaluation()
end program main_driver
