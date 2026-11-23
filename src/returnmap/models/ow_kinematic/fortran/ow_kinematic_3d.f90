! =============================================================================
! J2 plasticity with Ohno-Wang type kinematic hardening -- UMAT for 3D solids
!
! A port of the hand-written return mapping and consistent tangent of the
! Python model OWKinematic3D (../model.py); the formulation and the solution
! method are in ../README.md.
!
!   f        = sqrt(3/2 xi:xi) - sigma_y0,        xi = dev(sigma) - alpha
!   d(alpha) = (2/3) C_k n dp - gamma alpha_bar alpha dp
!   alpha_bar = sqrt(3/2 alpha:alpha),            n = (3/2) xi / xi_bar
!
! Backward Euler.  With w = 1 / (1 + gamma dlambda alpha_bar_{n+1}), a number
! in (0, 1], a plastic step is closed in w:
!
!   eta     = s_trial - w alpha_n            q = sqrt(3/2 eta:eta)
!   n       = (3/2) eta / q                  c = n : alpha_n
!   dlambda = (q - sigma_y0) / (3 mu + C_k w)
!   b       = alpha_n + (2/3) C_k dlambda n  A = sqrt(3/2 b:b)
!   sigma   = sigma_trial - 2 mu dlambda n
!   alpha   = w b,     p = p_n + dlambda
!
! and w is the root of the one scalar equation
!
!   g(w) = w - 1 + gamma dlambda(w) w^2 A(w) = 0
!
! solved by a safeguarded Newton iteration from w = 1 inside the bracket
! (0, 1): each update narrows the bracket with the current point, then takes
! the Newton step w - g / dg if it lands strictly inside the bracket and the
! midpoint of the bracket otherwise.  A point with q <= sigma_y0 (no positive
! dlambda) counts as g < 0 with no Newton step.  Converged when
! |g| < 1e-12; a step that does not get there in 100 updates, or whose
! solution has dlambda < 0 or is not finite, has not converged.
!
! Subroutines:
!   ow_kinematic_3d_reduced -- the quantities of the reduction at one w
!   ow_kinematic_3d         -- f2py-callable core: one stress update
!                              (matches OWKinematic3D)
!   umat                    -- ABAQUS UMAT entry point (unpacks PROPS and
!                              STATEV, calls ow_kinematic_3d)
!
! Material properties (PROPS):
!   PROPS(1) = E        Young's modulus
!   PROPS(2) = nu       Poisson's ratio
!   PROPS(3) = sigma_y0 yield stress (constant)
!   PROPS(4) = C_k      kinematic hardening modulus
!   PROPS(5) = gamma    dynamic recovery coefficient (1 / stress)
!
! State variables (STATEV):
!   STATEV(1:6) = alpha backstress (deviatoric), Voigt order below
!   STATEV(7)   = ep    equivalent plastic strain p
!
! Voigt order: [11, 22, 33, 12, 13, 23]
!   stress, alpha : shear components are tensor components
!   strain        : shear components are engineering shear strains (2 eps)
!
! Build with f2py:
!   returnmap build --family returnmap.models.ow_kinematic
! =============================================================================


! -----------------------------------------------------------------------------
! ow_kinematic_3d_reduced -- everything the reduction expresses in w, at one w
!
! Parameters
! ----------
! mu, sigma_y0, C_k, gamma [in]  : shear modulus and material parameters
! w                        [in]  : the reduced unknown, in (0, 1]
! s_trial(6)               [in]  : deviator of the trial stress
! alpha_n(6)               [in]  : backstress at the start of the increment
! below                    [out] : 1 when q <= sigma_y0 (no positive dlambda:
!                                  w lies below the root), else 0.  When 1,
!                                  only q is set; g = -1 and dg = 0.
! q                        [out] : sqrt(3/2 eta:eta), eta = s_trial - w alpha_n
! n_flow(6)                [out] : n = (3/2) eta / q
! c_na                     [out] : c = n : alpha_n
! denom                    [out] : D = 3 mu + C_k w
! dlambda                  [out] : (q - sigma_y0) / D
! b(6)                     [out] : alpha_n + (2/3) C_k dlambda n
! a_norm                   [out] : A = sqrt(3/2 b:b)
! k1, k2                   [out] : C_k dlambda (c + C_k dlambda) / A and
!                                  C_k dlambda^2 / A (both 0 when A = 0)
! g                        [out] : w - 1 + gamma dlambda w^2 A
! dg                       [out] : dg/dw = 1 + 2 gamma w dlambda A
!                                  - gamma w^2 ((A + k1)(c + C_k dlambda) / D
!                                               + k2 (alpha_bar_n^2 - c^2) / q)
! -----------------------------------------------------------------------------
subroutine ow_kinematic_3d_reduced(mu, sigma_y0, C_k, gamma, w, s_trial, alpha_n, &
                                   below, q, n_flow, c_na, denom, dlambda, b, &
                                   a_norm, k1, k2, g, dg)
    implicit none
    double precision, intent(in)  :: mu, sigma_y0, C_k, gamma, w
    double precision, intent(in)  :: s_trial(6), alpha_n(6)
    integer,          intent(out) :: below
    double precision, intent(out) :: q, n_flow(6), c_na, denom, dlambda, b(6)
    double precision, intent(out) :: a_norm, k1, k2, g, dg

    double precision :: eta(6), wt(6), acc, a2
    integer          :: i

    ! -- weights of the double contraction of two stress-like Voigt vectors
    do i = 1, 3
        wt(i) = 1.0d0
    end do
    do i = 4, 6
        wt(i) = 2.0d0
    end do

    do i = 1, 6
        n_flow(i) = 0.0d0
        b(i) = 0.0d0
    end do
    c_na = 0.0d0
    denom = 0.0d0
    dlambda = 0.0d0
    a_norm = 0.0d0
    k1 = 0.0d0
    k2 = 0.0d0
    g = -1.0d0
    dg = 0.0d0

    ! -- eta = s_trial - w alpha_n,  q = sqrt(3/2 eta:eta)
    acc = 0.0d0
    do i = 1, 6
        eta(i) = s_trial(i) - w * alpha_n(i)
        acc = acc + wt(i) * eta(i)**2
    end do
    q = sqrt(1.5d0 * acc)

    ! -- written so that a NaN q also counts as "below"
    below = 1
    if (q > sigma_y0) below = 0
    if (below == 1) return

    ! -- n = (3/2) eta / q,  c = n : alpha_n,  alpha_bar_n^2
    c_na = 0.0d0
    a2 = 0.0d0
    do i = 1, 6
        n_flow(i) = 1.5d0 * eta(i) / q
        c_na = c_na + wt(i) * n_flow(i) * alpha_n(i)
        a2 = a2 + wt(i) * alpha_n(i)**2
    end do
    a2 = 1.5d0 * a2

    denom = 3.0d0 * mu + C_k * w
    dlambda = (q - sigma_y0) / denom

    ! -- b = alpha_n + (2/3) C_k dlambda n,  A = sqrt(3/2 b:b)
    acc = 0.0d0
    do i = 1, 6
        b(i) = alpha_n(i) + (2.0d0 / 3.0d0) * C_k * dlambda * n_flow(i)
        acc = acc + wt(i) * b(i)**2
    end do
    a_norm = sqrt(1.5d0 * acc)

    if (a_norm > 0.0d0) then
        k1 = C_k * dlambda * (c_na + C_k * dlambda) / a_norm
        k2 = C_k * dlambda**2 / a_norm
    else
        k1 = 0.0d0
        k2 = 0.0d0
    end if

    g = w - 1.0d0 + gamma * dlambda * w**2 * a_norm
    dg = 1.0d0 + 2.0d0 * gamma * w * dlambda * a_norm &
         - gamma * w**2 * ((a_norm + k1) * (c_na + C_k * dlambda) / denom &
                           + k2 * (a2 - c_na**2) / q)

end subroutine ow_kinematic_3d_reduced


! -----------------------------------------------------------------------------
! ow_kinematic_3d -- one stress update of J2 + Ohno-Wang type kinematic
!                    hardening (3D)
!
! Parameters
! ----------
! E, nu, sigma_y0, C_k, gamma [in]  : material parameters
! stress_in(6)                [in]  : stress at the start of the increment
! alpha_in(6)                 [in]  : backstress at the start
! ep_in                       [in]  : equivalent plastic strain at the start
! dstran(6)                   [in]  : strain increment (engineering shear)
! stress_out(6)               [out] : updated stress
! alpha_out(6)                [out] : updated backstress
! ep_out                      [out] : updated equivalent plastic strain
! ddsdde(6,6)                 [out] : consistent tangent
! converged                   [out] : 1 = converged, 0 = not converged
!
! When the step does not converge (the iteration limit is reached, dlambda < 0,
! or the solution is not finite), stress_out, alpha_out and ep_out are the
! values at the start of the increment and ddsdde is the elastic stiffness.
! -----------------------------------------------------------------------------
subroutine ow_kinematic_3d(E, nu, sigma_y0, C_k, gamma, &
                           stress_in, alpha_in, ep_in, dstran, &
                           stress_out, alpha_out, ep_out, ddsdde, converged)
    implicit none
    double precision, intent(in)  :: E, nu, sigma_y0, C_k, gamma
    double precision, intent(in)  :: stress_in(6), alpha_in(6), ep_in, dstran(6)
    double precision, intent(out) :: stress_out(6), alpha_out(6), ep_out
    double precision, intent(out) :: ddsdde(6,6)
    integer,          intent(out) :: converged

    ! -- iteration on g(w) = 0
    double precision, parameter :: NEWTON_TOL = 1.0d-12
    integer,          parameter :: NEWTON_MAX_ITER = 100

    ! -- local variables
    double precision :: mu, lam, cmat(6,6), stress_trial(6), s_trial(6)
    double precision :: p_mean, acc, q_trial, f_trial
    double precision :: w, w_lo, w_hi, w_new
    double precision :: q, n_flow(6), c_na, denom, dlambda, b(6)
    double precision :: a_norm, k1, k2, g, dg
    double precision :: h(6), gvec(6), h_n, h_a, g_n, g_a
    double precision :: t_dev, t_nn, t_hg, vol_c, idev_c
    integer          :: i, j, n_iter, below
    logical          :: finite

    ! -- elastic moduli
    mu  = E / (2.0d0 * (1.0d0 + nu))
    lam = E * nu / ((1.0d0 + nu) * (1.0d0 - 2.0d0 * nu))

    ! -- elastic stiffness (6x6 Voigt, engineering shear strain)
    do j = 1, 6
        do i = 1, 6
            cmat(i,j) = 0.0d0
        end do
    end do
    do i = 1, 3
        do j = 1, 3
            cmat(i,j) = lam
        end do
        cmat(i,i) = lam + 2.0d0 * mu
    end do
    do i = 4, 6
        cmat(i,i) = mu
    end do

    ! -- the output of a step that does not converge: the state at the start
    !    of the increment and the elastic stiffness
    do i = 1, 6
        stress_out(i) = stress_in(i)
        alpha_out(i) = alpha_in(i)
    end do
    ep_out = ep_in
    do j = 1, 6
        do i = 1, 6
            ddsdde(i,j) = cmat(i,j)
        end do
    end do
    converged = 0

    ! -- trial stress: stress_trial = stress_in + C * dstran
    do i = 1, 6
        stress_trial(i) = stress_in(i)
        do j = 1, 6
            stress_trial(i) = stress_trial(i) + cmat(i,j) * dstran(j)
        end do
    end do

    ! -- deviator of the trial stress
    p_mean = (stress_trial(1) + stress_trial(2) + stress_trial(3)) / 3.0d0
    do i = 1, 6
        s_trial(i) = stress_trial(i)
    end do
    do i = 1, 3
        s_trial(i) = s_trial(i) - p_mean
    end do

    ! -- yield function at the trial state: xi = s_trial - alpha_n
    acc = 0.0d0
    do i = 1, 3
        acc = acc + (s_trial(i) - alpha_in(i))**2
    end do
    do i = 4, 6
        acc = acc + 2.0d0 * (s_trial(i) - alpha_in(i))**2
    end do
    q_trial = sqrt(1.5d0 * acc)
    f_trial = q_trial - sigma_y0

    if (f_trial <= 0.0d0) then
        ! ----------------------------------------------------------------
        ! Elastic step: the trial stress is the stress, elastic tangent
        ! ----------------------------------------------------------------
        do i = 1, 6
            stress_out(i) = stress_trial(i)
        end do
        converged = 1
        return
    end if

    ! --------------------------------------------------------------------
    ! Plastic step: safeguarded Newton iteration on g(w) = 0, from w = 1
    ! inside the bracket (w_lo, w_hi) = (0, 1)
    ! --------------------------------------------------------------------
    w = 1.0d0
    w_lo = 0.0d0
    w_hi = 1.0d0
    call ow_kinematic_3d_reduced(mu, sigma_y0, C_k, gamma, w, s_trial, alpha_in, &
                                 below, q, n_flow, c_na, denom, dlambda, b, &
                                 a_norm, k1, k2, g, dg)
    n_iter = 0
    if (abs(g) < NEWTON_TOL) converged = 1

    do while (converged == 0 .and. n_iter < NEWTON_MAX_ITER)
        ! -- narrow the bracket with the current point
        if (g > 0.0d0) then
            w_hi = w
        else
            w_lo = w
        end if
        ! -- Newton step if it lands strictly inside the bracket, else bisect
        w_new = w_lo
        if (below == 0 .and. dg /= 0.0d0) w_new = w - g / dg
        if (.not. (w_new > w_lo .and. w_new < w_hi)) w_new = 0.5d0 * (w_lo + w_hi)
        w = w_new
        call ow_kinematic_3d_reduced(mu, sigma_y0, C_k, gamma, w, s_trial, alpha_in, &
                                     below, q, n_flow, c_na, denom, dlambda, b, &
                                     a_norm, k1, k2, g, dg)
        n_iter = n_iter + 1
        if (abs(g) < NEWTON_TOL) converged = 1
    end do

    ! -- not converged: the output set above stands
    if (converged == 0) return

    ! -- the solution must be admissible (dlambda >= 0) and finite
    finite = abs(dlambda) <= huge(dlambda)
    do i = 1, 6
        finite = finite .and. abs(stress_trial(i) - 2.0d0 * mu * dlambda * n_flow(i)) &
                              <= huge(dlambda)
        finite = finite .and. abs(w * b(i)) <= huge(dlambda)
    end do
    if (.not. finite .or. dlambda < 0.0d0) then
        converged = 0
        return
    end if

    ! -- stress, backstress and equivalent plastic strain
    do i = 1, 6
        stress_out(i) = stress_trial(i) - 2.0d0 * mu * dlambda * n_flow(i)
        alpha_out(i) = w * b(i)
    end do
    ep_out = ep_in + dlambda

    ! --------------------------------------------------------------------
    ! Consistent tangent
    !   D = C - (3 mu dlambda / q) I_dev C
    !       - 4 mu^2 (1 / (3 mu + C_k w) - dlambda / q) n (x) n
    !       - (2 mu / dg) h (x) G
    !   h = (2 mu (c + C_k dlambda) / (3 mu + C_k w) - 2 mu dlambda c / q) n
    !       + (3 mu dlambda / q) alpha_n
    !   G = gamma w^2 ((A + k1) / (3 mu + C_k w) - k2 c / q) n
    !       + (3/2) gamma w^2 (k2 / q) alpha_n
    ! --------------------------------------------------------------------
    h_n = 2.0d0 * mu * (c_na + C_k * dlambda) / denom - 2.0d0 * mu * dlambda * c_na / q
    h_a = 3.0d0 * mu * dlambda / q
    g_n = gamma * w**2 * ((a_norm + k1) / denom - k2 * c_na / q)
    g_a = 1.5d0 * gamma * w**2 * (k2 / q)
    do i = 1, 6
        h(i) = h_n * n_flow(i) + h_a * alpha_in(i)
        gvec(i) = g_n * n_flow(i) + g_a * alpha_in(i)
    end do

    t_dev = 3.0d0 * mu * dlambda / q
    t_nn  = 4.0d0 * mu**2 * (1.0d0 / denom - dlambda / q)
    t_hg  = 2.0d0 * mu / dg

    do j = 1, 6
        ! (I_dev C)(i,j): C(i,j) less the mean of the direct rows for i = 1..3
        vol_c = (cmat(1,j) + cmat(2,j) + cmat(3,j)) / 3.0d0
        do i = 1, 6
            if (i <= 3) then
                idev_c = cmat(i,j) - vol_c
            else
                idev_c = cmat(i,j)
            end if
            ddsdde(i,j) = cmat(i,j) - t_dev * idev_c &
                          - t_nn * n_flow(i) * n_flow(j) &
                          - t_hg * h(i) * gvec(j)
        end do
    end do

end subroutine ow_kinematic_3d


! -----------------------------------------------------------------------------
! umat -- ABAQUS UMAT entry point
!
! Unpacks PROPS and STATEV and delegates to ow_kinematic_3d.  For 3D solid
! elements (NTENS = 6).  When the return mapping does not converge, PNEWDT is
! set to 0.5 and STRESS and STATEV are left as they were; when it converges,
! PNEWDT is left as it was.
! -----------------------------------------------------------------------------
subroutine umat(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
                RPL, DDSDDT, DRPLDE, DRPLDT, &
                STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, &
                PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, &
                NSTATV, PROPS, NPROPS, COORDS, DROT, PNEWDT, &
                CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, &
                KSPT, KSTEP, KINC)
    implicit none
    character(len=80),    intent(in)    :: CMNAME
    integer,              intent(in)    :: NDI, NSHR, NTENS, NSTATV, NPROPS
    integer,              intent(in)    :: NOEL, NPT, LAYER, KSPT, KSTEP, KINC
    double precision,     intent(inout) :: STRESS(NTENS)
    double precision,     intent(inout) :: STATEV(NSTATV)
    double precision,     intent(out)   :: DDSDDE(NTENS, NTENS)
    double precision,     intent(out)   :: SSE, SPD, SCD, RPL, DRPLDT
    double precision,     intent(out)   :: DDSDDT(NTENS), DRPLDE(NTENS)
    double precision,     intent(in)    :: STRAN(NTENS), DSTRAN(NTENS)
    double precision,     intent(in)    :: TIME(2), DTIME, TEMP, DTEMP
    double precision,     intent(in)    :: PREDEF(1), DPRED(1)
    double precision,     intent(in)    :: PROPS(NPROPS), COORDS(3)
    double precision,     intent(in)    :: DROT(3,3), DFGRD0(3,3), DFGRD1(3,3)
    double precision,     intent(inout) :: PNEWDT
    double precision,     intent(in)    :: CELENT

    double precision :: stress_n(6), alpha_n(6), ep_n, dstran_n(6)
    double precision :: stress_out(6), alpha_out(6), ep_out
    integer          :: i, converged

    do i = 1, 6
        stress_n(i) = STRESS(i)
        alpha_n(i) = STATEV(i)
        dstran_n(i) = DSTRAN(i)
    end do
    ep_n = STATEV(7)

    call ow_kinematic_3d(PROPS(1), PROPS(2), PROPS(3), PROPS(4), PROPS(5), &
                         stress_n, alpha_n, ep_n, dstran_n, &
                         stress_out, alpha_out, ep_out, DDSDDE, converged)

    if (converged == 1) then
        do i = 1, 6
            STRESS(i) = stress_out(i)
            STATEV(i) = alpha_out(i)
        end do
        STATEV(7) = ep_out
    else
        ! -- not converged: ask for a smaller increment, keep STRESS and STATEV
        PNEWDT = 0.5d0
    end if

    ! -- zero unused output fields
    SSE = 0.0d0; SPD = 0.0d0; SCD = 0.0d0; RPL = 0.0d0; DRPLDT = 0.0d0
    do i = 1, NTENS
        DDSDDT(i) = 0.0d0
        DRPLDE(i) = 0.0d0
    end do

end subroutine umat
