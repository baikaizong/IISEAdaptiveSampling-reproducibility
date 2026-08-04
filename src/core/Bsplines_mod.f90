!==============================================================
!  Two-dimensional B-spline basis generator on the unit square.
!  - Cubic B-spline: set degx=3, degy=3
!  - Open uniform knots
!  - Optional first derivatives w.r.t x and y
!  - Sampling: "midpoints" (default) or "endpoints"
!==============================================================
module Bsplines_mod
    use GlobalSettings_mod
    implicit none
    public 
    contains

      !------------------------------------------------------------
      ! One-dimensional B-spline basis via Cox-de Boor recursion
      ! The right endpoint belongs to the final knot span.
      !------------------------------------------------------------
      recursive pure function Bspline(i, k, x, t) result(val)
        integer, intent(in) :: i, k
        real(dp), intent(in) :: x, t(:)
        real(dp) :: val, term1, term2, denom1, denom2
        integer :: nknots

        nknots = size(t)

        if (k == 1) then
           if ((x >= t(i) .and. x < t(i+1)) .or. &
               (i+1 == nknots .and. abs(x - t(i+1)) <= EPS_END)) then
              val = 1.0d0
           else
              val = 0.0d0
           end if
        else
           val = 0.0d0
           denom1 = t(i+k-1) - t(i)
           denom2 = t(i+k)   - t(i+1)

           if (abs(denom1) > EPS_DEN) then
              term1 = (x - t(i))/denom1 * Bspline(i, k-1, x, t)
              val = val + term1
           end if

           if (abs(denom2) > EPS_DEN) then
              term2 = (t(i+k) - x)/denom2 * Bspline(i+1, k-1, x, t)
              val = val + term2
           end if
        end if
      end function Bspline

      !------------------------------------------------------------
      ! 1D B-spline first derivative w.r.t x
      !------------------------------------------------------------
      recursive pure function Bspline_d(i, k, x, t) result(dval)
        integer, intent(in) :: i, k
        real(dp), intent(in) :: x, t(:)
        real(dp) :: dval, denom1, denom2

        if (k == 1) then
           dval = 0.0d0
        else
           dval  = 0.0d0
           denom1 = t(i+k-1) - t(i)
           denom2 = t(i+k)   - t(i+1)

           if (abs(denom1) > EPS_DEN) dval = dval + (k-1)/denom1 * Bspline(i,   k-1, x, t)
           if (abs(denom2) > EPS_DEN) dval = dval - (k-1)/denom2 * Bspline(i+1, k-1, x, t)
        end if
      end function Bspline_d

      !------------------------------------------------------------
      ! Open uniform knot vector on [0,1]
      ! Knots: [0,...,0,  t_{inner}..., 1,...,1] with "order" repeats at ends
      ! ncoef: number of basis functions (coefficients)
      ! size(knot) must be ncoef+order
      !------------------------------------------------------------
      subroutine open_uniform_knots(order, ncoef, knot)
        integer, intent(in)  :: order, ncoef
        real(dp), intent(out) :: knot(:)
        integer :: i, nknots, ninter
        real(dp) :: step

        nknots = ncoef + order
        if (size(knot) /= nknots) stop 'open_uniform_knots: knot size mismatch'
        if (order <= 0 .or. ncoef <= 0) stop 'open_uniform_knots: invalid order/ncoef'

        ! clamp ends
        do i = 1, order
           knot(i) = 0.0d0
           knot(nknots - order + i) = 1.0d0
        end do

        ! number of interior knots
        ninter = nknots - 2*order   ! = ncoef - order

        if (ninter > 0) then
           step = 1.0d0 / dble(ninter + 1)
           do i = 1, ninter
              knot(order + i) = dble(i)*step
           end do
        end if
      end subroutine open_uniform_knots

      !------------------------------------------------------------
      ! Generate 2D tensor-product B-spline basis on [0,1]^2
      !
      ! Inputs:
      !   nx, ny   : evaluation grid size
      !   kx, ky   : # of basis in x/y (coefficients per dimension)
      !   degx,degy: polynomial degrees (e.g., 3 for cubic)
      ! Optional:
      !   basis_dx, basis_dy : first derivatives w.r.t x or y
      !   sample_mode: "midpoints" (default) or "endpoints"
      !
      ! Outputs:
      !   basis_val (nx*ny, kx*ky), and optionally basis_dx, basis_dy
      !
      ! Column order note:
      !   col = (j-1)*kx + i   with i fastest (x-direction)
      !------------------------------------------------------------
      subroutine gen_bspline2D(nx, ny, kx, ky, degx, degy, &
                               basis_val, basis_dx, basis_dy, sample_mode)
        integer, intent(in) :: nx, ny, kx, ky, degx, degy
        real(dp), intent(out) :: basis_val(nx*ny, kx*ky)
        real(dp), intent(out), optional :: basis_dx(nx*ny, kx*ky), basis_dy(nx*ny, kx*ky)
        character(*), intent(in), optional :: sample_mode

        integer :: orderx, ordery, nknx, nkny
        real(dp), allocatable :: knotx(:), knoty(:)
        integer :: ix, iy, i, j, row, col
        real(dp) :: x, y, bx, by, dbx, dby
        logical :: want_dx, want_dy, use_endpoints

        ! ----------------- input validations -----------------
        if (nx <= 0 .or. ny <= 0) stop 'gen_bspline2D: invalid nx/ny'
        if (kx <= 0 .or. ky <= 0) stop 'gen_bspline2D: invalid kx/ky'
        if (degx < 0 .or. degy < 0) stop 'gen_bspline2D: invalid degx/degy'

        orderx = degx + 1
        ordery = degy + 1

        if (kx < orderx) stop 'gen_bspline2D: require kx >= degx+1'
        if (ky < ordery) stop 'gen_bspline2D: require ky >= degy+1'

        want_dx = present(basis_dx)
        want_dy = present(basis_dy)

        use_endpoints = .false.
        if (present(sample_mode)) then
           if (len_trim(sample_mode) > 0) then
              if (trim(sample_mode) == 'endpoints' .or. trim(sample_mode) == 'Endpoints') use_endpoints = .true.
           end if
        end if

        ! ----------------- knots -----------------
        nknx = kx + orderx
        nkny = ky + ordery
        allocate(knotx(nknx), knoty(nkny))
        call open_uniform_knots(orderx, kx, knotx)
        call open_uniform_knots(ordery, ky, knoty)

        ! ----------------- compute basis -----------------
        basis_val = 0.0d0
        if (want_dx) basis_dx = 0.0d0
        if (want_dy) basis_dy = 0.0d0

        do iy = 1, ny
           if (use_endpoints) then
              if (ny == 1) then
                 y = 0.0d0
              else
                 y = dble(iy-1) / dble(ny-1)
              end if
           else
              y = (0.5d0 + dble(iy-1)) / dble(ny)
           end if

           do ix = 1, nx
              if (use_endpoints) then
                 if (nx == 1) then
                    x = 0.0d0
                 else
                    x = dble(ix-1) / dble(nx-1)
                 end if
              else
                 x = (0.5d0 + dble(ix-1)) / dble(nx)
              end if

              row = (iy-1)*nx + ix
              col = 0

              do j = 1, ky
                 by  = Bspline( j, ordery, y, knoty)
                 if (want_dy) dby = Bspline_d(j, ordery, y, knoty)
                 do i = 1, kx
                    bx  = Bspline( i, orderx, x, knotx)
                    if (want_dx) dbx = Bspline_d(i, orderx, x, knotx)

                    col = col + 1  

                    basis_val(row, col) = bx * by

                    if (want_dx) basis_dx(row, col) = dbx * by
                    if (want_dy) basis_dy(row, col) = bx  * dby
                 end do
              end do
           end do
        end do

        deallocate(knotx, knoty)
    end subroutine gen_bspline2D
                           
end module Bsplines_mod
