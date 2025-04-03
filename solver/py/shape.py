import numpy as np
from numpy.typing import NDArray
from typing import Tuple

class Shapes:
    def __init__(self, n_local: int):
        if n_local == 2:
            self.shape_fun = shape_fun_edge2
            self.shape_fun_derivatives = shape_fun_edge2_derivatives
        elif n_local == 3:
            self.shape_fun = shape_fun_p1
            self.shape_fun_derivatives = shape_fun_p1_derivatives
        elif n_local == 4:
            self.shape_fun = shape_fun_q1
            self.shape_fun_derivatives = shape_fun_q1_derivatives


#
# -------------------- Shape Functions --------------------
#
# Each returns phi or derivatives. For 2D shapes (tri/quad),
# we have (xsi, eta). For edges, only (xsi).

def shape_fun_q1(xsi: float, eta: float) -> NDArray[np.float64]:
    """
    Bilinear shape functions for a 4-node Q1 quad
    on [-1,1] x [-1,1].
    phi[0] -> node ( x=+1, y=+1 )
    phi[1] -> node ( x=-1, y=+1 )
    phi[2] -> node ( x=-1, y=-1)
    phi[3] -> node ( x=+1, y=-1)
    """
    phi = np.zeros(4, dtype=np.float64)
    phi[0] = (1.0 + xsi)*(1.0 + eta)/4.0
    phi[1] = (1.0 - xsi)*(1.0 + eta)/4.0
    phi[2] = (1.0 - xsi)*(1.0 - eta)/4.0
    phi[3] = (1.0 + xsi)*(1.0 - eta)/4.0
    return phi

def shape_fun_q1_derivatives(xsi: float, eta: float) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Derivatives wrt xsi, eta for the Q1 shape functions.
    Returns (dphidxsi, dphideta), each shape (4,).
    """
    dphidxsi = np.zeros(4, dtype=np.float64)
    dphideta = np.zeros(4, dtype=np.float64)

    dphidxsi[0] =  (1.0 + eta)/4.0
    dphidxsi[1] = -(1.0 + eta)/4.0
    dphidxsi[2] = -(1.0 - eta)/4.0
    dphidxsi[3] =  (1.0 - eta)/4.0

    dphideta[0] =  (1.0 + xsi)/4.0
    dphideta[1] =  (1.0 - xsi)/4.0
    dphideta[2] = -(1.0 - xsi)/4.0
    dphideta[3] = -(1.0 + xsi)/4.0

    return dphidxsi, dphideta

def shape_fun_p1(xsi: float, eta: float) -> NDArray[np.float64]:
    """
    Linear shape functions for a 3-node triangle
    with reference coords: (0,0), (1,0), (0,1).
    phi[0] = 1 - xsi - eta
    phi[1] = xsi
    phi[2] = eta
    """
    phi = np.zeros(3, dtype=np.float64)
    phi[0] = 1.0 - xsi - eta
    phi[1] = xsi
    phi[2] = eta
    return phi

def shape_fun_p1_derivatives(xsi: float, eta: float) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Derivatives wrt xsi, eta for the P1 triangle shape functions.
    Returns (dphidxsi, dphideta), each shape (3,).
    """
    dphidxsi  = np.array([-1.0, 1.0,  0.0], dtype=np.float64)
    dphideta  = np.array([-1.0, 0.0,  1.0], dtype=np.float64)
    return dphidxsi, dphideta

def shape_fun_edge2(xsi: float) -> NDArray[np.float64]:
    """
    Linear shape functions for a 2-node edge on [-1,1].
    phi[0] = (1 - xsi)/2
    phi[1] = (1 + xsi)/2
    """
    phi = np.zeros(2, dtype=np.float64)
    phi[0] = (1.0 - xsi)*0.5
    phi[1] = (1.0 + xsi)*0.5
    return phi

def shape_fun_edge2_derivatives(xsi: float) -> NDArray[np.float64]:
    """
    Derivative wrt xsi for the 2-node edge shape functions.
    """
    dphidxsi = np.array([-0.5, 0.5], dtype=np.float64)
    return dphidxsi
