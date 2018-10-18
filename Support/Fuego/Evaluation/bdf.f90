!
! BDF (backward differentiation formula) time-stepping routines.
!
! See
!
!   1. VODE: A variable-coefficient ODE solver; Brown, Byrne, and
!      Hindmarsh; SIAM J. Sci. Stat. Comput., vol. 10, no. 5, pp.
!      1035-1051, 1989.
!
!   2. An alternative implementation of variable step-size multistep
!      formulas for stiff ODES; Jackson and Sacks-Davis; ACM
!      Trans. Math. Soft., vol. 6, no. 3, pp. 295-318, 1980.
!
!   3. A polyalgorithm for the numerical solution of ODEs; Byrne and
!      Hindmarsh; ACM Trans. Math. Soft., vol. 1, no. 1, pp. 71-96,
!      1975.
!

module bdf
  implicit none

  integer, parameter  :: dp   = kind(1.d0)
  real(dp), parameter :: one  = 1.0_dp
  real(dp), parameter :: two  = 2.0_dp
  real(dp), parameter :: half = 0.5_dp

  integer, parameter :: bdf_max_iters = 666666666

  integer, parameter :: BDF_ERR_SUCCESS  = 0
  integer, parameter :: BDF_ERR_SOLVER   = 1
  integer, parameter :: BDF_ERR_MAXSTEPS = 2
  integer, parameter :: BDF_ERR_DTMIN    = 3

  character(len=64), parameter :: errors(0:3) = [ &
       'Success.                                                ', &
       'Newton solver failed to converge several times in a row.', &
       'Too many steps were taken.                              ', &
       'Minimum time-step reached several times in a row.       ' ]

  !
  ! bdf time-stepper
  !
  type :: bdf_ts

     ! parameters, set at build time by bdf_ts_build
     integer  :: neq                      ! number of equations (degrees of freedom) per point
     integer  :: npt                      ! number of points
     integer  :: max_order                ! maximum order (1 to 6)

     ! options, users are free to change these after bdf_ts_build
     integer  :: max_steps                ! maximum allowable number of steps
     integer  :: max_iters                ! maximum allowable number of newton iterations
     integer  :: verbose                  ! verbosity level
     real(dp) :: dt_min                   ! minimum allowable step-size
     real(dp) :: eta_min                  ! minimum allowable step-size shrink factor
     real(dp) :: eta_max                  ! maximum allowable step-size growth factor
     real(dp) :: eta_thresh               ! step-size growth threshold
     integer  :: max_j_age                ! maximum age of jacobian
     integer  :: max_p_age                ! maximum age of newton iteration matrix

     logical  :: debug
     integer  :: dump_unit

     real(dp), pointer :: rtol(:)         ! realtive tolerances
     real(dp), pointer :: atol(:)         ! absolute tolerances

     ! state (internal)
     real(dp) :: t                        ! current time
     real(dp) :: t1                       ! final time
     real(dp) :: dt                       ! current time step
     real(dp) :: dt_nwt                   ! dt used when building newton iteration matrix
     integer  :: k                        ! current order
     integer  :: n                        ! current step
     integer  :: j_age                    ! age of jacobian
     integer  :: p_age                    ! age of newton iteration matrix
     integer  :: k_age                    ! number of steps taken at current order
     real(dp) :: tq(-1:2)                 ! error coefficients (test quality)
     real(dp) :: tq2save
     logical  :: refactor

     real(dp), pointer :: J(:,:)          ! jacobian matrix
     real(dp), pointer :: P(:,:)          ! newton iteration matrix
     real(dp), pointer :: z(:,:,:)        ! nordsieck histroy array, indexed as (dof, p, n)
     real(dp), pointer :: z0(:,:,:)       ! nordsieck predictor array
     real(dp), pointer :: h(:)            ! time steps, h = [ h_n, h_{n-1}, ..., h_{n-k} ]
     real(dp), pointer :: l(:)            ! predictor/corrector update coefficients
     real(dp), pointer :: y(:,:)          ! current y
     real(dp), pointer :: yd(:,:)         ! current \dot{y}
     real(dp), pointer :: rhs(:,:)        ! solver rhs
     real(dp), pointer :: e(:,:)          ! accumulated correction
     real(dp), pointer :: e1(:,:)         ! accumulated correction, previous step
     real(dp), pointer :: ewt(:,:)        ! cached error weights
     real(dp), pointer :: b(:,:)          ! solver work space
     integer,  pointer :: ipvt(:)         ! pivots
     integer,  pointer :: A(:,:)          ! pascal matrix

     ! counters
     integer :: nfe                       ! number of function evaluations
     integer :: nje                       ! number of jacobian evaluations
     integer :: nlu                       ! number of factorizations
     integer :: nit                       ! number of non-linear solver iterations
     integer :: nse                       ! number of non-linear solver errors
     integer :: ncse                      ! number of consecutive non-linear solver errors
     integer :: ncit                      ! number of current non-linear solver iterations
     integer :: ncdtmin                   ! number of consecutive times we tried to shrink beyound the minimum time step

  end type bdf_ts

  private :: &
       rescale_timestep, decrease_order, increase_order, &
       alpha0, alphahat0, xi_j, xi_star_inv, ewts, norm, eye_r, eye_i, factorial

contains

  !
  ! Advance system from t0 to t1.
  !
  subroutine bdf_advance(ts, f, Jac, neq, npt, y0, t0, y1, t1, dt0, reset, reuse, ierr)
    type(bdf_ts), intent(inout) :: ts
    integer,      intent(in   ) :: neq, npt
    real(dp),     intent(in   ) :: y0(neq,npt), t0, t1, dt0
    real(dp),     intent(  out) :: y1(neq,npt)
    logical,      intent(in   ) :: reset, reuse
    integer,      intent(  out) :: ierr
    interface
       subroutine f(neq, npt, y, t, yd)
         import dp
         integer,  intent(in   ) :: neq, npt
         real(dp), intent(in   ) :: y(neq,npt), t
         real(dp), intent(  out) :: yd(neq,npt)
       end subroutine f
       subroutine Jac(neq, npt, y, t, J)
         import dp
         integer,  intent(in   ) :: neq, npt
         real(dp), intent(in   ) :: y(neq,npt), t
         real(dp), intent(  out) :: J(neq, neq)
       end subroutine Jac
    end interface

    integer  :: k
    logical  :: retry


    if (reset) then
       ts%t = t0
       call bdf_reset(ts, f, y0, dt0, reuse)
    end if

    ierr = BDF_ERR_SUCCESS

    ts%t1 = t1; ts%t = t0; ts%ncse = 0; ts%ncdtmin = 0;
    do k = 1, bdf_max_iters + 1
       if (ts%n > ts%max_steps .or. k > bdf_max_iters) then
          ierr = BDF_ERR_MAXSTEPS; return
       end if

       if (k == 1) &
            call bdf_dump(ts)

       call bdf_update(ts)                ! update various coeffs (l, tq) based on time-step history
       call bdf_predict(ts)               ! predict nordsieck array using pascal matrix
       call bdf_solve(ts, f, Jac)         ! solve for y_n based on predicted y and yd
       call bdf_check(ts, retry, ierr)    ! check for solver errors and test error estimate

       if (ierr /= BDF_ERR_SUCCESS) return
       if (retry) cycle

       call bdf_correct(ts)               ! new solution looks good, correct history and advance

       call bdf_dump(ts)
       if (ts%t >= t1) exit

       call bdf_adjust(ts)                ! adjust step-size/order
    end do

    if (ts%verbose > 0) &
         print '("BDF: n:",i6,", fe:",i6,", je: ",i3,", lu: ",i3,", it: ",i3,", se: ",i3,", dt: ",e15.8,", k: ",i2)', &
         ts%n, ts%nfe, ts%nje, ts%nlu, ts%nit, ts%nse, ts%dt, ts%k

    y1 = ts%z(:,:,0)

  end subroutine bdf_advance

  !
  ! Compute Nordsieck update coefficients l and error coefficients tq.
  !
  ! Regarding the l coefficients, see section 5, and in particular
  ! eqn. 5.2, of Jackson and Sacks-Davis (1980).
  !
  ! Regarding the error coefficients tq, these have been adapted from
  ! cvode.  The tq array is indexed as:
  !
  !  tq(-1) coeff. for order k-1 error est.
  !  tq(0)  coeff. for order k error est.
  !  tq(1)  coeff. for order k+1 error est.
  !  tq(2)  coeff. for order k+1 error est. (used for e_{n-1})
  !
  ! Note:
  !
  !   1. The input vector t = [ t_n, t_{n-1}, ... t_{n-k} ] where we
  !      are advancing from step n-1 to step n.
  !
  !   2. The step size h_n = t_n - t_{n-1}.
  !
  subroutine bdf_update(ts)
    type(bdf_ts), intent(inout) :: ts

    integer  :: j
    real(dp) :: a0, a0hat, a1, a2, a3, a4, a5, a6, xistar_inv, xi_inv, c

    ts%l  = 0
    ts%tq = 0

    ! compute l vector
    ts%l(0) = 1
    ts%l(1) = xi_j(ts%h, 1)
    if (ts%k > 1) then
       do j = 2, ts%k-1
          ts%l = ts%l + eoshift(ts%l, -1) / xi_j(ts%h, j)
       end do
       ts%l = ts%l + eoshift(ts%l, -1) * xi_star_inv(ts%k, ts%h)
    end if

    ! compute error coefficients (adapted from cvode)
    a0hat = alphahat0(ts%k, ts%h)
    a0    = alpha0(ts%k)

    xi_inv     = one
    xistar_inv = one
    if (ts%k > 1) then
       xi_inv     = one / xi_j(ts%h, ts%k)
       xistar_inv = xi_star_inv(ts%k, ts%h)
    end if

    a1 = one - a0hat + a0
    a2 = one + ts%k * a1
    ts%tq(0) = abs(a1 / (a0 * a2))
    ts%tq(2) = abs(a2 * xistar_inv / (ts%l(ts%k) * xi_inv))
    if (ts%k > 1) then
       c  = xistar_inv / ts%l(ts%k)
       a3 = a0 + one / ts%k
       a4 = a0hat + xi_inv
       ts%tq(-1) = abs(c * (one - a4 + a3) / a3)
    else
       ts%tq(-1) = one
    end if

    xi_inv = ts%h(0) / sum(ts%h(0:ts%k))
    a5 = a0 - one / (ts%k+1)
    a6 = a0hat - xi_inv
    ts%tq(1) = abs((one - a6 + a5) / a2 / (xi_inv * (ts%k+2) * a5))

    call ewts(ts)
  end subroutine bdf_update

  !
  ! Predict (apply Pascal matrix).
  !
  subroutine bdf_predict(ts)
    type(bdf_ts), intent(inout) :: ts
    integer :: i, j, m, p
    do i = 0, ts%k
       do p = 1, ts%npt
          ts%z0(:,p,i) = 0
          do j = i, ts%k
             do m = 1, ts%neq
                ts%z0(m,p,i) = ts%z0(m,p,i) + ts%A(i,j) * ts%z(m,p,j)
             end do
          end do
       end do
    end do
  end subroutine bdf_predict

  !
  ! Solve "y_n - dt f(y_n,t) = y - dt yd" for y_n where y and yd are
  ! predictors from the Nordsieck form.
  !
  ! Newton iteration is:
  !   solve:   P x = -c G(y(k)) for x
  !   update:  y(k+1) = y(k) + x
  ! where
  !   G(y) = y - dt * f(y,t) - rhs
  !
  subroutine bdf_solve(ts, f, Jac)
    type(bdf_ts), intent(inout) :: ts
    interface
       subroutine f(neq, npt, y, t, yd)
         import dp
         integer,  intent(in   ) :: neq, npt
         real(dp), intent(in   ) :: y(neq,npt), t
         real(dp), intent(  out) :: yd(neq,npt)
       end subroutine f
       subroutine Jac(neq, npt, y, t, J)
         import dp
         integer,  intent(in   ) :: neq, npt
         real(dp), intent(in   ) :: y(neq,npt), t
         real(dp), intent(  out) :: J(neq, neq)
       end subroutine Jac
    end interface

    include 'LinAlg.inc'

    integer  :: k, m, n, p, info
    real(dp) :: c, dt_adj, dt_rat, inv_l1
    logical  :: rebuild, iterating(ts%npt)

    inv_l1 = 1.0_dp / ts%l(1)
    do p = 1, ts%npt
       do m = 1, ts%neq
          ts%e(m,p)   = 0
          ts%rhs(m,p) = ts%z0(m,p,0) - ts%z0(m,p,1) * inv_l1
          ts%y(m,p)   = ts%z0(m,p,0)
       end do
    end do
    dt_adj = ts%dt / ts%l(1)

    dt_rat = dt_adj / ts%dt_nwt
    if (ts%p_age > ts%max_p_age) ts%refactor = .true.
    if (dt_rat < 0.7d0 .or. dt_rat > 1.429d0) ts%refactor = .true.

    iterating = .true.

    do k = 1, ts%max_iters

       ! build iteration matrix and factor
       if (ts%refactor) then
          rebuild = .true.
          if (ts%ncse == 0 .and. ts%j_age < ts%max_j_age) rebuild = .false.
          if (ts%ncse > 0  .and. (dt_rat < 0.2d0 .or. dt_rat > 5.d0)) rebuild = .false.

          if (rebuild) then
             call Jac(ts%neq, ts%npt, ts%y, ts%t, ts%J)
             ts%nje   = ts%nje + 1
             ts%j_age = 0
          end if

          call eye_r(ts%P)

          do m = 1, ts%neq
             do n = 1, ts%neq
                ts%P(n,m) = ts%P(n,m) - dt_adj * ts%J(n,m)
             end do
          end do

          call dgefa(ts%P, ts%neq, ts%neq, ts%ipvt, info)
          ! lapack      call dgetrf(neq, neq, ts%P, neq, ts%ipvt, info)
          ts%nlu    = ts%nlu + 1
          ts%dt_nwt = dt_adj
          ts%p_age  = 0
          ts%refactor  = .false.
       end if

       c = 2 * ts%dt_nwt / (dt_adj + ts%dt_nwt)

       call f(ts%neq, ts%npt, ts%y, ts%t, ts%yd)
       ts%nfe = ts%nfe + 1

       do p = 1, ts%npt
          if (.not. iterating(p)) cycle

          ! solve using factorized iteration matrix
          do m = 1, ts%neq
             ts%b(m,p) = c * (ts%rhs(m,p) - ts%y(m,p) + dt_adj * ts%yd(m,p))
          end do
          call dgesl(ts%P, ts%neq, ts%neq, ts%ipvt, ts%b(:,p), 0)
          ! lapack   call dgetrs ('N', neq, 1, ts%P, neq, ts%ipvt, ts%b, neq, info)
          ts%nit = ts%nit + 1

          do m = 1, ts%neq
             ts%e(m,p) = ts%e(m,p) + ts%b(m,p)
             ts%y(m,p) = ts%z0(m,p,0) + ts%e(m,p)
          end do
          if (norm(ts%b(:,p), ts%ewt(:,p)) < one) iterating(p) = .false.
       end do

       if (.not. any(iterating)) exit

    end do

    ts%ncit = k; ts%p_age = ts%p_age + 1; ts%j_age = ts%j_age + 1
  end subroutine bdf_solve

  !
  ! Check error estimates.
  !
  subroutine bdf_check(ts, retry, err)
    type(bdf_ts), intent(inout) :: ts
    logical,      intent(out)   :: retry
    integer,      intent(out)   :: err

    real(dp) :: error, eta
    integer  :: p

    retry = .false.; err = BDF_ERR_SUCCESS

    ! if solver failed many times, bail
    if (ts%ncit >= ts%max_iters .and. ts%ncse > 7) then
       err = BDF_ERR_SOLVER
       return
    end if

    ! if solver failed to converge, shrink dt and try again
    if (ts%ncit >= ts%max_iters) then
       ts%refactor = .true.; ts%nse = ts%nse + 1; ts%ncse = ts%ncse + 1
       call rescale_timestep(ts, 0.25d0)
       retry = .true.
       return
    end if
    ts%ncse = 0

    ! if local error is too large, shrink dt and try again
    do p = 1, ts%npt
       error = ts%tq(0) * norm(ts%e(:,p), ts%ewt(:,p))
       if (error > one) then
          eta = one / ( (6.d0 * error) ** (one / ts%k) + 1.d-6 )
          call rescale_timestep(ts, eta)
          retry = .true.
          if (ts%dt < ts%dt_min + epsilon(ts%dt_min)) ts%ncdtmin = ts%ncdtmin + 1
          if (ts%ncdtmin > 7) err = BDF_ERR_DTMIN
          return
       end if
    end do
    ts%ncdtmin = 0

  end subroutine bdf_check

  !
  ! Correct (apply l coeffs) and advance step.
  !
  subroutine bdf_correct(ts)
    type(bdf_ts), intent(inout) :: ts
    integer :: i, m, p

    do i = 0, ts%k
       do p = 1, ts%npt
          do m = 1, ts%neq
             ts%z(m,p,i) = ts%z0(m,p,i) + ts%e(m,p) * ts%l(i)
          end do
       end do
    end do

    ts%h     = eoshift(ts%h, -1)
    ts%h(0)  = ts%dt
    ts%t     = ts%t + ts%dt
    ts%n     = ts%n + 1
    ts%k_age = ts%k_age + 1
  end subroutine bdf_correct


  !
  ! Dump (for debugging)...
  !
  subroutine bdf_dump(ts)
    type(bdf_ts), intent(inout) :: ts
    integer :: i, m, p

    if (.not. ts%debug) return
    write(ts%dump_unit,*) ts%t, ts%z(:,:,0)
  end subroutine bdf_dump

  !
  ! Adjust step-size/order to maximize step-size.
  !
  subroutine bdf_adjust(ts)
    type(bdf_ts), intent(inout) :: ts

    real(dp) :: c, error, eta(-1:1), rescale, etamax(ts%npt), etaminmax, delta(ts%npt)
    integer  :: p

    rescale = 0

    do p = 1, ts%npt
       ! compute eta(k-1), eta(k), eta(k+1)
       eta = 0
       error  = ts%tq(0) * norm(ts%e(:,p), ts%ewt(:,p))
       eta(0) = one / ( (6.d0 * error) ** (one / ts%k) + 1.d-6 )
       if (ts%k_age > ts%k) then
          if (ts%k > 1) then
             error     = ts%tq(-1) * norm(ts%z(:,p,ts%k), ts%ewt(:,p))
             eta(-1) = one / ( (6.d0 * error) ** (one / ts%k) + 1.d-6 )
          end if
          if (ts%k < ts%max_order) then
             c = (ts%tq(2) / ts%tq2save) * (ts%h(0) / ts%h(2)) ** (ts%k+1)
             error  = ts%tq(1) * norm(ts%e(:,p) - c * ts%e1(:,p), ts%ewt(:,p))
             eta(1) = one / ( (10.d0 * error) ** (one / (ts%k+2)) + 1.d-6 )
          end if
          ts%k_age = 0
       end if

       ! choose which eta will maximize the time step
       etamax(p) = 0
       if (eta(-1) > etamax(p)) then
          etamax(p) = eta(-1)
          delta(p)  = -1
       end if
       if (eta(1) > etamax(p)) then
          etamax(p) = eta(1)
          delta(p)  = 1
       end if
       if (eta(0) > etamax(p)) then
          etamax(p) = eta(0)
          delta(p)  = 0
       end if
    end do

    p = minloc(etamax, dim=1)
    rescale = 0
    etaminmax = etamax(p)
    if (etaminmax > ts%eta_thresh) then
       if (delta(p) == -1) then
          call decrease_order(ts)
       else if (delta(p) == 1) then
          call increase_order(ts)
       end if
       rescale = etaminmax
    end if

    if (ts%t + ts%dt > ts%t1) then
       rescale = (ts%t1 - ts%t) / ts%dt
       call rescale_timestep(ts, rescale, .true.)
    else if (rescale /= 0) then
       call rescale_timestep(ts, rescale)
    end if

    ! save for next step (needed to compute eta(1))
    ts%e1 = ts%e
    ts%tq2save = ts%tq(2)

  end subroutine bdf_adjust

  !
  ! Reset counters, set order to one, init Nordsieck history array.
  !
  subroutine bdf_reset(ts, f, y0, dt, reuse)
    type(bdf_ts), intent(inout) :: ts
    real(dp),     intent(in   ) :: y0(ts%neq, ts%npt), dt
    logical,      intent(in   ) :: reuse
    interface
       subroutine f(neq, npt, y, t, yd)
         import dp
         integer,  intent(in   ) :: neq, npt
         real(dp), intent(in   ) :: y(neq,npt), t
         real(dp), intent(  out) :: yd(neq,npt)
       end subroutine f
    end interface

    ts%nfe = 0
    ts%nje = 0
    ts%nlu = 0
    ts%nit = 0
    ts%nse = 0

    ts%y  = y0
    ts%dt = dt
    ts%n  = 1
    ts%k  = 1

    ts%h        = ts%dt
    ts%dt_nwt   = ts%dt
    ts%refactor = .true.

    call f(ts%neq, ts%npt, ts%y, ts%t, ts%yd)
    ts%nfe = ts%nfe + 1

    ts%z(:,:,0) = ts%y
    ts%z(:,:,1) = ts%dt * ts%yd

    ts%k_age = 0
    if (.not. reuse) then
       ts%j_age = ts%max_j_age + 1
       ts%p_age = ts%max_p_age + 1
    else
       ts%j_age = 0
       ts%p_age = 0
    end if

  end subroutine bdf_reset

  !
  ! Rescale time-step.
  !
  ! This consists of:
  !   1. bound eta to honor eta_min, eta_max, and dt_min
  !   2. scale dt and adjust time array t accordingly
  !   3. rescale Nordsieck history array
  !
  subroutine rescale_timestep(ts, eta_in, force_in)
    type(bdf_ts), intent(inout)           :: ts
    real(dp),     intent(in   )           :: eta_in
    logical,      intent(in   ), optional :: force_in

    real(dp) :: eta
    integer  :: i
    logical  :: force

    force = .false.; if (present(force_in)) force = force_in

    if (force) then
       eta = eta_in
    else
       eta = max(eta_in, ts%dt_min / ts%dt, ts%eta_min)
       eta = min(eta, ts%eta_max)

       if (ts%t + eta*ts%dt > ts%t1) then
          eta = (ts%t1 - ts%t) / ts%dt
       end if
    end if

    ts%dt   = eta * ts%dt
    ts%h(0) = ts%dt

    do i = 1, ts%k
       ts%z(:,:,i) = eta**i * ts%z(:,:,i)
    end do
  end subroutine rescale_timestep

  !
  ! Decrease order.
  !
  subroutine decrease_order(ts)
    type(bdf_ts), intent(inout) :: ts
    integer  :: j
    real(dp) :: c(0:6)

    if (ts%k > 2) then
       c = 0
       c(2) = 1
       do j = 1, ts%k-2
          c = eoshift(c, -1) + c * xi_j(ts%h, j)
       end do

       do j = 2, ts%k-1
          ts%z(:,:,j) = ts%z(:,:,j) - c(j) * ts%z(:,:,ts%k)
       end do
    end if

    ts%z(:,:,ts%k) = 0
    ts%k = ts%k - 1
  end subroutine decrease_order

  !
  ! Increase order.
  !
  subroutine increase_order(ts)
    type(bdf_ts), intent(inout) :: ts
    integer  :: j
    real(dp) :: c(0:6)

    c = 0
    c(2) = 1
    do j = 1, ts%k-2
       c = eoshift(c, -1) + c * xi_j(ts%h, j)
    end do

    ts%z(:,:,ts%k+1) = 0
    do j = 2, ts%k+1
       ts%z(:,:,j) = ts%z(:,:,j) + c(j) * ts%e
    end do

    ts%k = ts%k + 1
  end subroutine increase_order

  !
  ! Return $\alpha_0$.
  !
  function alpha0(k) result(a0)
    integer,  intent(in) :: k
    real(dp) :: a0
    integer  :: j
    a0 = -1
    do j = 2, k
       a0 = a0 - one / j
    end do
  end function alpha0

  !
  ! Return $\hat{\alpha}_{n,0}$.
  !
  function alphahat0(k, h) result(a0)
    integer,  intent(in) :: k
    real(dp), intent(in) :: h(0:k)
    real(dp) :: a0
    integer  :: j
    a0 = -1
    do j = 2, k
       a0 = a0 - h(0) / sum(h(0:j-1))
    end do
  end function alphahat0

  !
  ! Return 1 / $\xi^*_k$.
  !
  ! Note that a lot of simplifications can be made to the formula for
  ! $\xi^*_k$ that appears in Jackson and Sacks-Davis.
  !
  function xi_star_inv(k, h) result(xii)
    integer,  intent(in) :: k
    real(dp), intent(in) :: h(0:)
    real(dp) :: xii, hs
    integer  :: j
    hs = 0.0_dp
    xii = -alpha0(k)
    do j = 0, k-2
       hs  = hs + h(j)
       xii = xii - h(0) / hs
    end do
  end function xi_star_inv

  !
  ! Return $\xi_j$.
  !
  function xi_j(h, j) result(xi)
    integer,  intent(in) :: j
    real(dp), intent(in) :: h(0:)
    real(dp) :: xi
    xi = sum(h(0:j-1)) / h(0)
  end function xi_j

  !
  ! Pre-compute error weights.
  !
  subroutine ewts(ts)
    type(bdf_ts), intent(inout) :: ts
    integer :: m, p
    do p = 1, ts%npt
       do m = 1, ts%neq
          ts%ewt(m,p) = one / (ts%rtol(m) * abs(ts%y(m,p)) + ts%atol(m))
       end do
    end do
  end subroutine ewts

  subroutine print_y(ts)
    type(bdf_ts), intent(in) :: ts
    integer :: p
    do p = 1, ts%npt
       print *, ts%y(:,p)
    end do
  end subroutine print_y

  !
  ! Compute weighted norm of y.
  !
  function norm(y, ewt) result(r)
    real(dp), intent(in) :: y(1:), ewt(1:)
    real(dp) :: r
    integer :: m, n
    n = size(y)
    r = 0.0_dp
    do m = 1, n
       r = r + (y(m)*ewt(m))**2
    end do
    r = sqrt(r/n)
  end function norm

  !
  ! Build/destroy BDF time-stepper.
  !
  subroutine bdf_ts_build(ts, neq, npt, rtol, atol, max_order)
    type(bdf_ts), intent(inout) :: ts
    integer,      intent(in   ) :: max_order, neq, npt
    real(dp),     intent(in   ) :: rtol(neq), atol(neq)

    integer :: k, U(max_order+1, max_order+1), Uk(max_order+1, max_order+1)

    allocate(ts%rtol(neq))
    allocate(ts%atol(neq))
    allocate(ts%z(neq, npt, 0:max_order))
    allocate(ts%z0(neq, npt, 0:max_order))
    allocate(ts%l(0:max_order))
    allocate(ts%h(0:max_order))
    allocate(ts%A(0:max_order, 0:max_order))
    allocate(ts%P(neq, neq))
    allocate(ts%J(neq, neq))
    allocate(ts%y(neq, npt))
    allocate(ts%yd(neq, npt))
    allocate(ts%rhs(neq, npt))
    allocate(ts%e(neq, npt))
    allocate(ts%e1(neq, npt))
    allocate(ts%ewt(neq, npt))
    allocate(ts%b(neq, npt))
    allocate(ts%ipvt(neq))

    ts%neq        = neq
    ts%npt        = npt
    ts%max_order  = max_order
    ts%max_steps  = 1000000
    ts%max_iters  = 10
    ts%verbose    = 0
    ts%dt_min     = epsilon(ts%dt_min)
    ts%eta_min    = 0.2_dp
    ts%eta_max    = 10.0_dp
    ts%eta_thresh = 1.50_dp
    ts%max_j_age  = 50
    ts%max_p_age  = 20

    ts%k = -1

    ts%rtol = rtol
    ts%atol = atol

    ts%J  = 0
    ts%P  = 0
    ts%yd = 0

    ts%j_age = 666666666
    ts%p_age = 666666666

    ts%debug = .false.

    ! build pascal matrix A using A = exp(U)
    U = 0
    do k = 1, max_order
       U(k,k+1) = k
    end do
    Uk = U
    call eye_i(ts%A)
    do k = 1, max_order+1
       ts%A  = ts%A + Uk / factorial(k)
       Uk = matmul(U, Uk)
    end do
  end subroutine bdf_ts_build

  subroutine bdf_ts_destroy(ts)
    type(bdf_ts), intent(inout) :: ts
    deallocate(ts%h,ts%l,ts%ewt,ts%rtol,ts%atol)
    deallocate(ts%y,ts%yd,ts%z,ts%z0,ts%A)
    deallocate(ts%P,ts%J,ts%rhs,ts%e,ts%e1,ts%b,ts%ipvt)
  end subroutine bdf_ts_destroy

  !
  ! Various misc. helper functions
  !
  subroutine eye_r(A)
    real(dp), intent(inout) :: A(:,:)
    integer :: i
    A = 0
    do i = 1, size(A, 1)
       A(i,i) = 1
    end do
  end subroutine eye_r
  subroutine eye_i(A)
    integer, intent(inout) :: A(:,:)
    integer :: i
    A = 0
    do i = 1, size(A, 1)
       A(i,i) = 1
    end do
  end subroutine eye_i
  recursive function factorial(n) result(r)
    integer, intent(in) :: n
    integer :: r
    if (n == 1) then
       r = 1
    else
       r = n * factorial(n-1)
    end if
  end function factorial

end module bdf
