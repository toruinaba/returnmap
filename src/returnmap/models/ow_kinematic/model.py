"""J2 plasticity with Ohno-Wang modified Armstrong-Frederick kinematic hardening.

This model is the canonical example of a *genuinely implicit* hardening law
in returnmap.  Unlike the standard Armstrong-Frederick model, the backward-Euler
discretization of the Ohno-Wang evolution equation cannot be solved in closed
form, so the model declares its states ``Implicit`` and writes
``state_residual``; they join σ and Δλ as Newton-Raphson unknowns.

Model parameters
----------------
E        : Young's modulus
nu       : Poisson's ratio
sigma_y0 : Initial yield stress (constant — no isotropic hardening)
C_k      : Kinematic hardening modulus
gamma    : Dynamic recovery parameter (Ohno-Wang nonlinearity), in 1/stress

Internal state
--------------
alpha : backstress tensor (shape (ntens,), deviatoric)
ep    : equivalent plastic strain (scalar, ≥ 0)

Yield function
--------------
f(state) = vonmises_norm(dev(σ) − α) − σ_y0

Ohno-Wang backstress evolution
-------------------------------
Standard Armstrong-Frederick has linear dynamic recovery:

    dα = (2/3) C_k dp n̂ − γ dp α

The Ohno-Wang variant weights the recovery by the current backstress norm:

    dα = (2/3) C_k dp n̂ − γ dp vonmises_norm(α) α

Backward-Euler discretisation as residual:

    R_α = α_{n+1} − α_n − (2/3) C_k Δλ n̂ + γ Δλ vonmises_norm(α_{n+1}) α_{n+1} = 0

Dimensional notes
-----------------
The three concrete classes share identical yield_function / state_residual
implementations via the OWKinematic parent.  Each stress-state dimension
provides dimension-specific dev / vonmises_norm that encapsulate
per-dimension arithmetic.

* **3D** (OWKinematic3D): α is a 6-component deviatoric tensor.
* **PS** (OWKinematicPS): α stores [α11, α22, α12]; vonmises_norm reconstructs
  α33 = −(α11+α22) without constructing a 6-component intermediate vector.
* **1D** (OWKinematic1D): α stores the deviatoric component α11_dev;
  vonmises_norm gives (3/2)|α11_dev|.

Notes
-----
* gamma=0 reduces to Prager's linear kinematic hardening (same limit as AF).
* σ, α and ep are all Newton-Raphson unknowns: ``x = [σ, Δλ, α, ep]``.
* ``OWKinematic3D`` adds a hand-written return mapping — one scalar unknown,
  ``w = 1 / (1 + γ Δλ ᾱ_{n+1})`` — and its consistent tangent, and has a UMAT
  (``fortran/ow_kinematic_3d.f90``).  ``README.md`` gives the reduction.
"""

import numpy as np

from returnmap.core.material import MaterialModel, ParamMeta
from returnmap.core.result import ReturnMappingResult
from returnmap.core.state import Implicit, NTENS, SCALAR
from returnmap.core.dimension import SOLID_3D, PLANE_STRESS, UNIAXIAL_1D


class OWKinematic(MaterialModel):
    """J2 + Ohno-Wang kinematic hardening — common physics across stress states.

    Pass ``dimension=`` to select the stress state, or use one of the
    pre-built concrete classes (:class:`OWKinematic3D`, :class:`OWKinematicPS`,
    :class:`OWKinematic1D`), which fix it.

    σ, α and ep are implicit: the Newton-Raphson unknowns are ``[σ, Δλ, α, ep]``.

    Parameters
    ----------
    dimension : StressDimension, optional
        Defaults to ``SOLID_3D``.
    E, nu, sigma_y0, C_k, gamma : float
        Material parameters.
    """

    stress = Implicit(shape=NTENS, doc="Cauchy stress")
    param_names = ["E", "nu", "sigma_y0", "C_k", "gamma"]
    param_meta = {
        "E": ParamMeta(doc="Young's modulus", unit="MPa",
                       default=210_000.0, bounds=(1e3, 1e6)),
        "nu": ParamMeta(doc="Poisson's ratio", default=0.3, bounds=(0.0, 0.49)),
        "sigma_y0": ParamMeta(doc="Initial yield stress", unit="MPa",
                              default=250.0, bounds=(1.0, 5e3)),
        "C_k": ParamMeta(doc="Kinematic hardening modulus (Ohno-Wang)", unit="MPa",
                         default=20_000.0, bounds=(0.0, 1e6)),
        # γ multiplies ‖α‖·α, so it carries 1/stress — unlike the AF γ.
        "gamma": ParamMeta(doc="Critical-state coefficient (Ohno-Wang)", unit="1/MPa",
                           default=100.0, bounds=(0.0, 1e4)),
    }
    alpha = Implicit(shape=NTENS, doc="backstress tensor (deviatoric)")
    ep = Implicit(shape=SCALAR, doc="equivalent plastic strain")

    def yield_function(self, state):
        s_xi = self.dev(state["stress"]) - state["alpha"]
        return self.vonmises_norm(s_xi) - self.sigma_y0

    def state_residual(self, new, dlambda, inc):
        alpha_new = new["alpha"]
        # α is stress-like, so the evolution law takes the physical-shear flow.
        n_hat = self.flow(new).stress_like
        a_norm = self.vonmises_norm(alpha_new)
        R_stress = self.default_stress_residual(new, dlambda, inc.trial["stress"])
        R_alpha = (
            alpha_new - inc.state_n["alpha"]
            - (2.0 / 3.0) * self.C_k * dlambda * n_hat
            + self.gamma * dlambda * a_norm * alpha_new
        )
        R_ep = new["ep"] - inc.state_n["ep"] - dlambda
        return [self.stress(R_stress), self.alpha(R_alpha), self.ep(R_ep)]


class OWKinematic3D(OWKinematic):
    """J2 + Ohno-Wang type kinematic hardening for 3D solids, solved by hand.

    Inherits the equations of :class:`OWKinematic` unchanged and adds the
    hand-written return mapping (``user_defined_return_mapping``) and its
    consistent tangent (``user_defined_tangent``).  Both rest on the identity
    ``C n = 2 mu n`` for a deviatoric direction, which holds for a full-rank
    stress state, so the class is fixed to ``SOLID_3D``.  The UMAT
    ``fortran/ow_kinematic_3d.f90`` is a port of these two methods.

    The reduction (``README.md``, section "手書きの解法 (3D)").  Write
    ``w = 1 / (1 + gamma dlambda alpha_bar_{n+1})``, a number in (0, 1].  The
    backward-Euler equations then give everything in closed form in w::

        eta     = s_trial - w alpha_n          q = sqrt(3/2 eta:eta)
        n       = (3/2) eta / q                c = n : alpha_n
        dlambda = (q - sigma_y0) / (3 mu + C_k w)
        b       = alpha_n + (2/3) C_k dlambda n,    A = sqrt(3/2 b:b)
        sigma   = sigma_trial - 2 mu dlambda n
        alpha   = w b,     alpha_bar = w A,     p = p_n + dlambda

    and the definition of w is the one scalar equation left::

        g(w) = w - 1 + gamma dlambda(w) w^2 A(w) = 0

    Parameters
    ----------
    E : float
        Young's modulus.
    nu : float
        Poisson's ratio.
    sigma_y0 : float
        Yield stress (constant).
    C_k : float
        Kinematic hardening modulus.
    gamma : float
        Dynamic recovery coefficient.
    """

    dimension = SOLID_3D

    #: Iteration on g(w) = 0: absolute tolerance on |g| (g is dimensionless,
    #: of order one) and the largest number of updates.
    NEWTON_TOL = 1.0e-12
    NEWTON_MAX_ITER = 100

    def _reduced(self, w, s_trial, alpha_n):
        """Everything the reduction expresses in w, at one value of w.

        Returns ``None`` when ``q <= sigma_y0``: the plastic multiplier that
        goes with this w is not positive, so w lies below the root.

        Otherwise a dict with ``q``, ``n``, ``c``, ``dlambda``, ``b``, ``A``,
        the denominator ``D = 3 mu + C_k w``, the two quotients
        ``k1 = C_k dlambda (c + C_k dlambda) / A`` and
        ``k2 = C_k dlambda^2 / A`` (both zero when ``A = 0``), the value
        ``g`` of the scalar equation and its derivative ``dg``::

            dg = 1 + 2 gamma w dlambda A
                 - gamma w^2 ((A + k1) (c + C_k dlambda) / D
                              + k2 (alpha_bar_n^2 - c^2) / q)
        """
        mu = self.E / (2.0 * (1.0 + self.nu))
        C_k = self.C_k
        gamma = self.gamma

        eta = s_trial - w * alpha_n
        q = np.sqrt(1.5 * _contract(eta, eta))
        if not q > self.sigma_y0:
            return None
        n = 1.5 * eta / q
        c = _contract(n, alpha_n)
        D = 3.0 * mu + C_k * w
        dlambda = (q - self.sigma_y0) / D
        b = alpha_n + (2.0 / 3.0) * C_k * dlambda * n
        A = np.sqrt(1.5 * _contract(b, b))
        if A > 0.0:
            k1 = C_k * dlambda * (c + C_k * dlambda) / A
            k2 = C_k * dlambda ** 2 / A
        else:
            k1 = 0.0
            k2 = 0.0
        a2 = 1.5 * _contract(alpha_n, alpha_n)

        g = w - 1.0 + gamma * dlambda * w ** 2 * A
        dg = (
            1.0 + 2.0 * gamma * w * dlambda * A
            - gamma * w ** 2 * (
                (A + k1) * (c + C_k * dlambda) / D + k2 * (a2 - c ** 2) / q
            )
        )
        return {
            "q": q, "n": n, "c": c, "D": D, "dlambda": dlambda, "b": b, "A": A,
            "k1": k1, "k2": k2, "g": g, "dg": dg,
        }

    def user_defined_return_mapping(self, inc):
        """Safeguarded Newton iteration on the scalar equation g(w) = 0.

        The root is bracketed from the start: ``g(1) >= 0`` on a plastic step
        and ``g(w) < 0`` wherever ``dlambda(w) <= 0``, in particular as
        ``w -> 0``.  From ``w = 1``, each update first narrows the bracket
        ``(w_lo, w_hi)`` with the current point (``g > 0`` replaces ``w_hi``,
        otherwise ``w_lo``), then takes the Newton step ``w - g / dg`` if it
        lands strictly inside the bracket and the midpoint of the bracket
        otherwise.  A point with ``q <= sigma_y0`` counts as ``g < 0`` with no
        Newton step.  The iteration stops when ``|g| < 1e-12``; a step that
        does not get there in 100 updates, or whose solution has
        ``dlambda < 0`` or is not finite, is reported as not converged, with
        the state at the start of the increment.

        With ``gamma = 0`` the start ``w = 1`` is the root and no update is
        made.

        Reads sigma_trial from ``inc.trial["stress"]``, and ``alpha`` and
        ``ep`` from ``inc.state_n``.

        Returns
        -------
        ReturnMappingResult
            ``residual_history`` holds ``|g|`` at the start and after every
            update (1.0 for a point with ``q <= sigma_y0``); ``n_iterations``
            is the number of updates.
        """
        mu = self.E / (2.0 * (1.0 + self.nu))
        stress_trial = np.asarray(inc.trial["stress"], dtype=float)
        alpha_n = np.asarray(inc.state_n["alpha"], dtype=float)
        ep_n = float(inc.state_n["ep"])
        s_trial = np.asarray(self.dev(stress_trial), dtype=float)

        w = 1.0
        w_lo = 0.0
        w_hi = 1.0
        r = self._reduced(w, s_trial, alpha_n)
        g = r["g"] if r is not None else -1.0
        residual_history = [abs(g)]
        n_iterations = 0
        converged = abs(g) < self.NEWTON_TOL
        while not converged and n_iterations < self.NEWTON_MAX_ITER:
            if g > 0.0:
                w_hi = w
            else:
                w_lo = w
            w_new = w_lo
            if r is not None and r["dg"] != 0.0:
                w_new = w - g / r["dg"]
            if not (w_lo < w_new < w_hi):
                w_new = 0.5 * (w_lo + w_hi)
            w = w_new
            r = self._reduced(w, s_trial, alpha_n)
            g = r["g"] if r is not None else -1.0
            residual_history.append(abs(g))
            n_iterations += 1
            converged = abs(g) < self.NEWTON_TOL

        if converged:
            # The solution must be admissible (dlambda >= 0) and finite.
            dlambda = r["dlambda"]
            stress_new = stress_trial - 2.0 * mu * dlambda * r["n"]
            alpha_new = w * r["b"]
            converged = bool(
                np.isfinite(dlambda)
                and np.all(np.isfinite(stress_new))
                and np.all(np.isfinite(alpha_new))
                and dlambda >= 0.0
            )

        if not converged:
            # Not a solution: hand back the state at the start of the increment.
            stress_n = np.asarray(inc.state_n["stress"], dtype=float)
            return ReturnMappingResult(
                stress=stress_n,
                state={"stress": stress_n, "alpha": alpha_n, "ep": np.array(ep_n)},
                dlambda=np.array(0.0),
                n_iterations=n_iterations,
                residual_history=residual_history,
                converged=False,
            )

        return ReturnMappingResult(
            stress=stress_new,
            state={
                "stress": stress_new,
                "alpha": alpha_new,
                "ep": np.array(ep_n + dlambda),
            },
            dlambda=np.array(dlambda),
            n_iterations=n_iterations,
            residual_history=residual_history,
            converged=True,
        )

    def user_defined_tangent(self, new, dlambda, inc):
        """Consistent tangent of the hand-written return mapping, in closed form::

            D = C - (3 mu dlambda / q) I_dev C
                - 4 mu^2 (1 / (3 mu + C_k w) - dlambda / q) n (x) n
                - (2 mu / dg) h (x) G

        The first three terms are the tangent at fixed w.  The last one is
        the change of w with the trial stress, ``dw = -(G : ds_trial) / dg``,
        carried to the stress by ``h = d(sigma)/dw``::

            h = (2 mu (c + C_k dlambda) / (3 mu + C_k w)
                 - 2 mu dlambda c / q) n + (3 mu dlambda / q) alpha_n
            G = gamma w^2 ((A + k1) / (3 mu + C_k w) - k2 c / q) n
                + (3/2) gamma w^2 (k2 / q) alpha_n

        n, alpha_n, h and G are stress-like Voigt vectors (tensor shear
        components); a column of D multiplies an engineering shear strain, and
        ``t : ds_trial = 2 mu t . d(eps)`` for a deviatoric t, which is where
        the plain outer products come from.  ``w`` is recovered from the
        converged state, ``w = 1 / (1 + gamma dlambda alpha_bar_{n+1})``, and
        every other quantity is the one of the return mapping at that w.

        A step the return mapping did not converge on (it reports
        ``dlambda = 0``, which no plastic solution has) gets the elastic
        stiffness.

        Returns
        -------
        ndarray, shape (6, 6)
        """
        mu = self.E / (2.0 * (1.0 + self.nu))
        gamma = self.gamma
        C = np.asarray(self.elastic_stiffness(inc.state_n), dtype=float)
        s_trial = np.asarray(self.dev(inc.trial["stress"]), dtype=float)
        alpha_n = np.asarray(inc.state_n["alpha"], dtype=float)
        alpha_new = np.asarray(new["alpha"], dtype=float)

        if not float(dlambda) > 0.0:
            return C

        alpha_bar = np.sqrt(1.5 * _contract(alpha_new, alpha_new))
        w = 1.0 / (1.0 + gamma * float(dlambda) * alpha_bar)
        r = self._reduced(w, s_trial, alpha_n)
        q, n, c, D = r["q"], r["n"], r["c"], r["D"]
        dl, A, k1, k2, dg = r["dlambda"], r["A"], r["k1"], r["k2"], r["dg"]

        h = (
            (2.0 * mu * (c + self.C_k * dl) / D - 2.0 * mu * dl * c / q) * n
            + (3.0 * mu * dl / q) * alpha_n
        )
        G = (
            gamma * w ** 2 * ((A + k1) / D - k2 * c / q) * n
            + 1.5 * gamma * w ** 2 * (k2 / q) * alpha_n
        )
        I_dev = np.asarray(self.I_dev(), dtype=float)

        return (
            C
            - (3.0 * mu * dl / q) * (I_dev @ C)
            - 4.0 * mu ** 2 * (1.0 / D - dl / q) * np.outer(n, n)
            - (2.0 * mu / dg) * np.outer(h, G)
        )


class OWKinematicPS(OWKinematic):
    """J2 + Ohno-Wang kinematic hardening for plane-stress elements.

    σ, α and ep are implicit: the Newton-Raphson unknowns are ``[σ, Δλ, α, ep]``.

    Parameters
    ----------
    E, nu, sigma_y0, C_k, gamma : float
        Material parameters.
    """

    dimension = PLANE_STRESS


class OWKinematic1D(OWKinematic):
    """J2 + Ohno-Wang kinematic hardening for uniaxial elements.

    σ, α and ep are implicit: the Newton-Raphson unknowns are ``[σ, Δλ, α, ep]``.

    The stored ``alpha`` is the deviatoric component α11_dev.  The effective
    backstress (yield-surface centre) is α_eff = (3/2)·α11_dev.

    Parameters
    ----------
    E, nu, sigma_y0, C_k, gamma : float
        Material parameters.
    """

    dimension = UNIAXIAL_1D


#: Weights of the double contraction a:b of two stress-like Voigt vectors
#: (tensor shear components, each standing for two entries of the tensor).
_CONTRACTION_WEIGHTS = np.array([1.0, 1.0, 1.0, 2.0, 2.0, 2.0])


def _contract(a, b):
    """Double contraction a:b of two stress-like 3D Voigt vectors."""
    return float(np.dot(a * _CONTRACTION_WEIGHTS, b))
