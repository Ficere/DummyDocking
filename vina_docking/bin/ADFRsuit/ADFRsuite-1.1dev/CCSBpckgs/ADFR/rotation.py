import numpy as np


def normalize(v):
    return v / np.sqrt(np.sum(v**2))


def quaternion_rotation(pt, quaternion):
    # print 'INPUT:', pt, 'Q:', quaternion, 'OUT:', qv_mult(quaternion, pt)
    return qv_mult(quaternion, pt)


def quaternion_rotation_wrapper(coords, quaternion):
    return np.array([quaternion_rotation(coord, quaternion) for coord in coords])


def q_mult(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
    z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
    return w, x, y, z


def q_conjugate(q):
    w, x, y, z = q
    return (w, -x, -y, -z)


def qv_mult(q1, v1):
    q2 = (0.0,) + tuple(v1)
    return q_mult(q_mult(q1, q2), q_conjugate(q1))[1:]


def axisangle_to_q(v, theta):
    v = normalize(v)
    x, y, z = v
    theta /= 2.0
    w = np.cos(theta)
    x = x * np.sin(theta)
    y = y * np.sin(theta)
    z = z * np.sin(theta)
    return w, x, y, z
