"""What the ow_kinematic family is verified on.

``MODELS`` are the models the family's tests run the common checks on
(:func:`returnmap.verification.check_model`): each class at the typical
parameters it documents, and at the typical parameters with ``gamma = 1`` —
a weak dynamic recovery, under which the backstress grows to the order of
the yield stress and its transient spans the loading histories (at the
typical ``gamma = 100`` it saturates within two or three steps).  The
standard loading histories are enough for this family, so no ``histories``
function is defined.  ``UMATS`` are the Fortran ports:
``returnmap build --family returnmap.models.ow_kinematic`` compiles them, and
``check_model`` compares a built one with the Python model step by step.
"""

from returnmap.models.ow_kinematic.model import OWKinematic3D, OWKinematicPS, OWKinematic1D
from returnmap.verification import Umat

MODELS = [
    OWKinematic3D.at_defaults(),
    OWKinematic3D.at_defaults(gamma=1.0),
    OWKinematicPS.at_defaults(),
    OWKinematicPS.at_defaults(gamma=1.0),
    OWKinematic1D.at_defaults(),
    OWKinematic1D.at_defaults(gamma=1.0),
]


def _umat_failure(ret):
    """The core returns ``(stress, alpha, ep, ddsdde, converged)``; 1 means converged."""
    return None if int(ret[4]) == 1 else "the UMAT's iteration on g(w) = 0 did not converge"


UMATS = [
    Umat(OWKinematic3D, module="ow_kinematic_3d",
         sources=["fortran/ow_kinematic_3d.f90"],
         integrator_options={"parse_umat_failure": _umat_failure}),
]
