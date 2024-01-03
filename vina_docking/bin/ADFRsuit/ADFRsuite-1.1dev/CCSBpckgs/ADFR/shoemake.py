import numpy as np
from math import sqrt, sin, cos, pi

TWOPI = pi * 2


def cube3_to_quaternion(xyz):
    """map (x,y,z) cube inside [1,2pi,2pi] to quaternion (w,x,y,z)
    if x, y and z are uniform in the range [0, 1/2pi], the
    quaternions will be evenly distributed
    http://planning.cs.uiuc.edu/node198.html
    https://doi.org/10.1016/B978-0-08-050755-2.50036-1
    Also, this is how quaternions are generated in the original autodock.
    """

    # expand input for readability
    u1, u2, u3 = xyz

    # quaterion = (w, x, y, z)
    w = sqrt(1 - u1) * sin(u2)
    x = sqrt(1 - u1) * cos(u2)
    y = sqrt(u1) * sin(u3)
    z = sqrt(u1) * cos(u3)

    return w, x, y, z


def quaternion_to_cube3(wxyz, wrap=True):
    """
    graphical solution:
      - y and z define a point within a 2D circle of radius 1.0, where
        sqrt(u1) is the distance from (0,0,0). Therefore, u1 is
        y**2 + z**2.
    """

    if abs(np.sum(np.array(wxyz) ** 2) - 1.0) > 1e-5:
        raise RuntimeError("input quaternion is not normalized")

    w, x, y, z = wxyz

    u1 = min(y**2 + z**2, 1.0)
    u1_ = -(w**2 + x**2) + 1

    if abs(u1 - u1_) > 1e-5:
        raise RuntimeError("expected u1 to be equal to u1_")

    ## if u1 == 0.:
    ##     print "WARNING: numerical singularity in quaternion_to_cube3"
    ##     print "WARNING: impossible to determine u3"

    ## if u1 == 1.:
    ##     print "WARNING: numerical singularity in quaternion_to_cube3"
    ##     print "WARNING: impossible to determine u2"

    # arctan2 returns values [-pi, +pi]
    u2 = np.arctan2(w, x)
    u3 = np.arctan2(y, z)

    if not wrap:
        return u1, u2, u3

    if u2 < 0.0:
        u2 += TWOPI
    if u3 < 0.0:
        u3 += TWOPI

    return 1.0, u2, u3
