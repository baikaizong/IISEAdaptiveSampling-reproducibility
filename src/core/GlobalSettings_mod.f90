!======================================================================
!  Module: GlobalSettings
!  Purpose:
!     Centralized parameter control for simulation and sampling
!     - All parameters are SAVE variables
!     - set_params(): optional initialization of all parameters
!     - get_params(): return all parameters and print to file handle
!======================================================================
module GlobalSettings_mod
    implicit none
    save
    !==============================Global Parameters===================    
    integer, parameter :: dp = kind(1.0d0)          
    real(dp), parameter :: eps_perturb = 1.0d-12             ! small random perturbation to break ties
    real(dp), parameter :: eps_tol=1.0d-14                   ! eps_tol: threshold for treating very small entries as zero
    real(dp), parameter :: EPS_DEN = 1d-14
    real(dp), parameter :: EPS_END = 1d-12
    real(dp), parameter :: zero = 0.0_dp
    real(dp), parameter :: theta_coverage = 0.827d0, pi_value=3.14159265d0         
    

    
    !---------------------- Grid parameters ----------------------
    integer :: num_x = 100, num_y = 100
    integer:: num_disnearspatial=5000
    integer :: num_samplingnodes=10,num_samplingnodes_addone=10+1
    integer :: num_allnodes=100*100,num_allnodes_addone=100*100+1
    
    !---------------------- Simulation control -------------------
    integer :: simu = 5000, Maxarl = 2000, Time_window = 50
    integer :: IC_TestRuns=10000, down_gap=150, up_gap=50, IC_runs = 1000000
    !---------------------- In-control baseline ------------------
    real(dp) :: IcArl = 370.0d0, Icstd = 3.0d0, Exploration_IcArl=100.0d0
    real(dp) :: Df = 4.0d0
    
    !---------------------- MPI message tags ---------------------
    integer, parameter :: TAG_WORK_REQUEST = 1
    integer, parameter :: TAG_TASK_ASSIGNMENT = 2
    integer, parameter :: TAG_RESULT = 3
    integer, parameter :: TAG_TERMINATE = 4

    !---------------------- NOISE identifiers ---------------------
    integer, parameter :: NOISE_GAUSSIAN      = 1
    integer, parameter :: NOISE_EXPONENTIAL   = 2
    integer, parameter :: NOISE_CHISQUARE     = 3
    integer, parameter :: NOISE_CORR_GAUSSIAN = 4
    
    ! --- Background type constants ---
    integer, parameter :: BACK_NONE          = 0
    integer, parameter :: BACK_BSPLINE       = 1 
    integer, parameter :: BACK_KERNEL_RANDOM = 2 

    ! --- Anomaly Type Constants ---
    integer, parameter :: TYPE_CIRCLE          = 1
    integer, parameter :: TYPE_SPATIAL         = 2
    integer, parameter :: TYPE_ST              = 3
    integer, parameter :: TYPE_BSPLINE         = 4
    integer, parameter :: TYPE_RANDOM_POINTS_CONV     = 5
    integer, parameter :: TYPE_CIRCLE_CONV = TYPE_RANDOM_POINTS_CONV
    integer, parameter :: TYPE_ELLIPSE         = 6
    integer, parameter :: TYPE_CRESCENT        = 7
    
    integer, parameter :: num_disnear=200
    integer, parameter :: IC = 1, OC = 2  


    !---------------------- Noise parameters ---------------------
    real(dp) :: sigma_ground = 3.0d0, sigma_noise = 1.0d0, sigma_kernel=5.0d0
    real(dp) :: kernel_bandwidth=1.0d0
    logical :: initialized_general = .false.    

    !---------------------- Bspline parameters ------------------
    integer :: kx0_TS=5, ky0_TS=5, degx0_TS=2, degy0_TS=2
    integer :: kx1_TS=20, ky1_TS=20, degx1_TS=3, degy1_TS=3
    integer :: knot0_square_TS=5*5, knot1_square_TS=20*20
      
    !---------------------- CDS settings ---------------    
    integer ::  num_varnear_CDS = 50 
    
    !--------------------------------------------------------------
    !  AnomalyParams_type
    !  Unified structure to pass anomaly parameters, reducing
    !  argument list complexity in the main subroutine.
    !--------------------------------------------------------------
    type :: AnomalyParams_type
            integer  :: type_id              ! Anomaly type identifier
            real(dp) :: center_idx(2)        ! Index of the center node
            real(dp) :: radius               ! Initial radius
            real(dp) :: value                ! Shift value or Amplitude
            real(dp) :: area
            real(dp) :: ellipticity
            real(dp) :: theta                ! Spatial parameter
            real(dp) :: delta_r              ! Rate of change for radius (ST type)
            real(dp) :: delta_t              ! Rate of change for time/intensity
            integer  :: time_idx             ! Current time step index
            integer  :: bspline_idx          ! Column index for B-spline anomaly
        
            ! --- Random-point anomaly fields ---
            integer  :: num_points           ! Number of random points
            real(dp), allocatable :: points_coords(:,:) ! Coordinates of points (2, num_points)
        end type AnomalyParams_type  
    ! Structure to hold statistical results
    type :: ARL_Stats_type
        real(dp) :: mean
        real(dp) :: std_err     ! Standard Error of Mean
        real(dp) :: min_val
        real(dp) :: max_val
        real(dp) :: q1          ! 25% Quantile
        real(dp) :: median      ! 50% Quantile
        real(dp) :: q3          ! 75% Quantile
    end type ARL_Stats_type
    
contains
    !==================================================================
    !  Subroutine: set_params
    !  Purpose: Initialize or update parameters (all optional inputs)
    !==================================================================
    subroutine set_generalparams(num_x_in, num_y_in, num_samplingnodes_in, &
        num_disnearspatial_in, simu_in, Maxarl_in, Time_window_in, IC_runs_in, &
        IcArl_in, Icstd_in, Df_in, sigma_ground_in, sigma_noise_in)

        implicit none
        !------------- Optional Inputs -----------------
        integer,  intent(in), optional :: num_x_in, num_y_in,num_samplingnodes_in, num_disnearspatial_in
        integer,  intent(in), optional :: simu_in, Maxarl_in, Time_window_in
        integer,  intent(in), optional :: IC_runs_in
        real(dp),  intent(in), optional :: IcArl_in, Icstd_in, Df_in
        real(dp),  intent(in), optional :: sigma_ground_in, sigma_noise_in

        !------------- Assignments ----------------------
        if (present(num_x_in)) num_x = num_x_in
        if (present(num_y_in)) num_y = num_y_in
        if (present(num_samplingnodes_in)) num_samplingnodes = num_samplingnodes_in
        if (present(num_disnearspatial_in)) num_disnearspatial = num_disnearspatial_in
        if (present(simu_in))  simu  = simu_in
        if (present(Maxarl_in)) Maxarl = Maxarl_in
        if (present(Time_window_in)) Time_window = Time_window_in
        if (present(IC_runs_in)) IC_runs = IC_runs_in

        if (present(IcArl_in)) IcArl = IcArl_in
        if (present(Icstd_in)) Icstd = Icstd_in
        if (present(Df_in)) Df = Df_in

        if (present(sigma_ground_in)) sigma_ground = sigma_ground_in
        if (present(sigma_noise_in)) sigma_noise = sigma_noise_in

        !------------- Derived values -------------------
        num_samplingnodes_addone=num_samplingnodes+1
        num_allnodes = num_x * num_y
        num_allnodes_addone=num_allnodes+1
        initialized_general = .true.
    end subroutine set_generalparams


    !==================================================================
    !  Subroutine: get_params
    !  Purpose: Return and print all parameters to file handle
    !==================================================================
    subroutine get_generalparams(fid)
        implicit none
        integer, intent(in) :: fid  ! File handle

        write(fid,'(A)') "==============================================="
        write(fid,'(A)') "        GLOBAL PARAMETER CONFIGURATION"
        write(fid,'(A)') "==============================================="

        write(fid,'(A,I8)') "num_x               = ", num_x
        write(fid,'(A,I8)') "num_y               = ", num_y
        write(fid,'(A,I10)') "num_allnodes         = ", num_allnodes
        write(fid,'(A,I8)') "simu                = ", simu
        write(fid,'(A,I8)') "Maxarl              = ", Maxarl
        write(fid,'(A,I8)') "Time_window         = ", Time_window
        write(fid,'(A,I12)') "IC_runs             = ", IC_runs
        write(fid,'(A,F12.6)') "IcArl              = ", IcArl
        write(fid,'(A,F12.6)') "Icstd              = ", Icstd
        write(fid,'(A,F12.6)') "Df                 = ", Df
        write(fid,'(A,F12.6)') "sigma_ground        = ", sigma_ground
        write(fid,'(A,F12.6)') "sigma_noise         = ", sigma_noise

        write(fid,'(A)') "==============================================="
        write(fid,'(A)') "End of GLOBAL PARAMETER listing"
        write(fid,'(A)') "==============================================="
end subroutine get_generalparams

subroutine set_Bsplineparams( &
    kx0_TS_in, ky0_TS_in, degx0_TS_in, degy0_TS_in, &
    kx1_TS_in, ky1_TS_in, degx1_TS_in, degy1_TS_in)

    implicit none
    !------------------ Optional Inputs ------------------
    integer, intent(in), optional :: kx0_TS_in, ky0_TS_in, degx0_TS_in, degy0_TS_in
    integer, intent(in), optional :: kx1_TS_in, ky1_TS_in, degx1_TS_in, degy1_TS_in


    !------------------ Parameter Assignments ------------------
    if (present(kx0_TS_in))  kx0_TS  = kx0_TS_in
    if (present(ky0_TS_in))  ky0_TS  = ky0_TS_in
    if (present(degx0_TS_in)) degx0_TS = degx0_TS_in
    if (present(degy0_TS_in)) degy0_TS = degy0_TS_in

    if (present(kx1_TS_in))  kx1_TS  = kx1_TS_in
    if (present(ky1_TS_in))  ky1_TS  = ky1_TS_in
    if (present(degx1_TS_in)) degx1_TS = degx1_TS_in
    if (present(degy1_TS_in)) degy1_TS = degy1_TS_in

    !------------------ Derived values ------------------
    knot0_square_TS = kx0_TS * ky0_TS
    knot1_square_TS = kx1_TS * ky1_TS

end subroutine set_Bsplineparams

subroutine get_Bsplineparams(fid)
    implicit none
    integer, intent(in) :: fid  ! File handle for output

    write(fid,'(A)') "==============================================="
    write(fid,'(A)') "           Bspline PARAMETER CONFIGURATION"
    write(fid,'(A)') "==============================================="

    write(fid,'(A,I8)') "kx0_TS              = ", kx0_TS
    write(fid,'(A,I8)') "ky0_TS              = ", ky0_TS
    write(fid,'(A,I8)') "degx0_TS            = ", degx0_TS
    write(fid,'(A,I8)') "degy0_TS            = ", degy0_TS
    write(fid,'(A,I8)') "kx1_TS              = ", kx1_TS
    write(fid,'(A,I8)') "ky1_TS              = ", ky1_TS
    write(fid,'(A,I8)') "degx1_TS            = ", degx1_TS
    write(fid,'(A,I8)') "degy1_TS            = ", degy1_TS

    write(fid,'(A)') "==============================================="
    write(fid,'(A)') "End of Bspline PARAMETER listing"
    write(fid,'(A)') "==============================================="

end subroutine get_Bsplineparams
    
subroutine GeneratePoints(Nodesset)
        implicit none
        real(dp), intent(out) :: Nodesset(num_allnodes, 2)  ! Output: 2D array storing (x, y) coordinates
        integer :: i, j
        real(dp) :: x_start, y_start, x_step, y_step

        ! Ensure num_x and num_y are valid to prevent division errors
        if (num_x <= 0 .or. num_y <= 0) then
            print *, "Error: num_x and num_y must be greater than zero!"
            return
        end if

        ! Compute the step size for x and y directions
        x_step = 1.0 / num_x
        y_step = 1.0 / num_y

        ! Compute the starting positions (center of the first grid cell)
        x_start = x_step / 2.0
        y_start = y_step / 2.0

        ! Generate a uniform grid of points within the unit square (0,0) to (1,1)
        do i = 1, num_x
            do j = 1, num_y
                Nodesset((i - 1) * num_y + j, 1) = x_start + (i - 1) * x_step  ! X-coordinate
                Nodesset((i - 1) * num_y + j, 2) = y_start + (j - 1) * y_step  ! Y-coordinate
            end do
        end do
end subroutine GeneratePoints


end module GlobalSettings_mod
