# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: base
#     language: python
#     name: python3
# ---

import numpy as np
import scipy as sp
import sympy as smp
from scipy.integrate import solve_ivp
import snoop
import matplotlib.pyplot as plt
import itertools
import mpmath
import plotly.graph_objects as go
import plotly.express as px

"""
Reference Paper:
https://doi.org/10.1046/j.1365-8711.1999.02754.x

"""


# %load_ext snoop

# %load_ext line_profiler

def get_norm(T,g):
    """
    T should have confuguration either 'u' or 'uu
    Returned Configuration: Scalar
    """
    if T.ndim==2:
        return np.einsum("ua,vb,ab,uv->",g,g,T,T)
    
    else:
        return np.einsum("ab,b,a->",g,T,T)


def levi_cevita_tensor(dim):
    """
    Define Levi Cevita Tensor
    Config: llll
    """

    arr = np.zeros(tuple([dim for _ in range(dim)]))
    for x in itertools.permutations(tuple(range(dim))):
        mat = np.zeros((dim, dim), dtype=np.int32)
        for i, j in zip(range(dim), x):
            mat[i, j] = 1
        arr[x] = int(np.linalg.det(mat))
    return arr


def antisymmetrize(arr):
    #* Function to antisymmetrize matrix
    arr = np.swapaxes(arr, -2, -1)
    anti_arr = arr - np.swapaxes(arr, -2, -1)
    anti_tensor=np.swapaxes(anti_arr,-2,-1)
    return anti_tensor


def symmetrize(arr):
    """
    Symmetrizes a 4-dimensional numpy array along last two axes.

    Args:
        arr (numpy.ndarray): The input array.

    Returns:
        numpy.ndarray: The symmetrized array.

    Example:

        arr = np.array([[[1, 2, 3, 4],
                         [5, 6, 7, 8],
                         [9, 10, 11, 12],
                         [13, 14, 15, 16]]])

        result = symmetrize(arr)
        print(result)
    """

    toput = np.diagonal(arr, 0, 1, 2)
    de_ = np.zeros((4, 4, 4),dtype=arr.dtype)


    idx = np.arange(de_.shape[0])
    de_[:, idx, idx] = toput

    return arr+np.transpose(arr,(0,2,1))-de_


def metric_tensor(r, theta,a,M=1):
    """
    Define the metric tensor function.

    Parameters:
        a (float): Parameter 'a' in the metric.
        r (float): Parameter 'r' in the metric.
        theta (float): Parameter 'theta' in the metric.

    Returns:
        np.ndarray: The metric tensor at the given values of a, r, and theta.
        Configuration: ll
    """
    g = np.array(
        [
            [
                (-(a**2) + 2 * M * r - r**2 + a**2 * np.sin(theta) ** 2)
                / (r**2 + a**2 * np.cos(theta) ** 2),
                0,
                0,
                -(
                    (2 * a * M * r * np.sin(theta) ** 2)
                    / (r**2 + a**2 * np.cos(theta) ** 2)
                ),
            ],
            [
                0,
                (r**2 + a**2 * np.cos(theta) ** 2) / (a**2 - 2 * M * r + r**2),
                0,
                0,
            ],
            [0, 0, r**2 + a**2 * np.cos(theta) ** 2, 0],
            [
                -(
                    (2 * a * M * r * np.sin(theta) ** 2)
                    / (r**2 + a**2 * np.cos(theta) ** 2)
                ),
                0,
                0,
                (
                    np.sin(theta) ** 2
                    * (
                        (a**2 + r**2) ** 2
                        - a**2 * (a**2 - 2 * M * r + r**2) * np.sin(theta) ** 2
                    )
                )
                / (r**2 + a**2 * np.cos(theta) ** 2),
            ],
        ]
    )

    return g


def kerr_christoffel(r, theta, a, M=1):

    """
    Function to give the christoffel symbols for kerr metric.
    The christoffel symbols are given as \Gamma ^i _{jk}

    From Reference Paper, Appendix
    Config: ull
    """
    
    cs = np.zeros((4, 4, 4))

    # Definitions
    Delta = a**2 - 2 * M * r + r**2
    scA = (r**2 + a**2) ** 2 - Delta * a**2 * np.sin(theta) ** 2
    omega_k = 2 * M * a * r / scA
    Sigma = r**2 + a**2 * np.cos(theta) ** 2

    cs[3, 0, 1] = M * (2 * r**2 - Sigma) / (Delta * Sigma**2) * a
    cs[0, 0, 1] = cs[3, 0, 1] * (r**2 + a**2) / a

    cs[3, 0, 2] = -2 * M * a * r / (np.tan(theta) * Sigma**2)
    cs[0, 0, 2] = a * np.sin(theta) ** 2 * cs[3, 0, 2]

    cs[0, 1, 3] = (
        -M
        * a
        * (2 * r**2 * (r**2 + a**2) + Sigma * (r**2 - a**2))
        * np.sin(theta) ** 2
        / (Delta * Sigma**2)
    )

    cs[3, 1, 3] = (
        r * Sigma * (Sigma - 2 * M * r)
        - M * a**2 * (2 * r**2 - Sigma) * np.sin(theta) ** 2
    ) / (Delta * Sigma**2)

    cs[0, 2, 3] = M * a**3 * r * np.sin(theta) ** 2 * np.sin(2 * theta) / Sigma**2

    cs[3, 2, 3] = (scA - Sigma * a**2 * np.sin(theta) ** 2) / (
        np.tan(theta) * Sigma**2
    )

    cs[1, 0, 0] = M * Delta * (2 * r**2 - Sigma) / Sigma**3

    cs[1, 0, 3] = -cs[1, 0, 0] * a * np.sin(theta) ** 2

    cs[2, 0, 0] = -M * a * r * np.sin(2 * theta) / Sigma**3 * a
    cs[2, 0, 3] = -cs[2, 0, 0] * (r**2 + a**2) / a

    cs[1, 1, 1] = r / Sigma + (M - r) / Delta

    cs[1, 2, 2] = -r * Delta / Sigma
    cs[2, 1, 2] = -cs[1, 2, 2] / Delta

    cs[1, 1, 2] = cs[2, 2, 2] = -(a**2) * np.sin(2 * theta) / (2 * Sigma)
    cs[2, 1, 1] = -cs[1, 1, 2] / Delta

    cs[1, 3, 3] = (
        -Delta
        * (r * Sigma**2 - M * a**2 * (2 * r**2 - Sigma) * np.sin(theta) ** 2)
        * np.sin(theta) ** 2
        / Sigma**3
    )

    cs[2, 3, 3] = (
        -(Delta * Sigma**2 + 2 * M * r * (r**2 + a**2) ** 2)
        * np.sin(2 * theta)
        / (2 * Sigma**3)
    )

    return symmetrize(arr=cs)
    # return cs


def RRC(r, theta, a, M):
    Delta = r**2 + a**2 - 2 * M * r
    Sigma = r**2 + a**2 * np.cos(theta) ** 2

    Sigma_32 = np.power(Sigma, 3 / 2)
    Delta_12 = np.power(Delta, 1 / 2)

    omega = np.zeros((4, 4, 4))

    w1 = (
        r * a**2 * np.sin(theta) ** 2 - M * r**2 + M * a**2 * np.cos(theta) ** 2
    ) / (Sigma_32 * Delta_12)
    w2 = a * r * np.sin(theta) / Sigma_32
    w3 = r * Delta_12 / Sigma_32
    w4 = a**2 * np.cos(theta) * np.sin(theta) / Sigma_32
    w5 = a * np.cos(theta) * Delta_12 / Sigma_32
    w6 = (r**2 + a**2) * np.cos(theta) / (Sigma_32 * np.sin(theta))

    omega[0, 1, 0] = omega[0, 0, 0] = w1

    omega[3, 1, 0] = omega[3, 0, 1] = omega[1, 3, 0] = omega[1, 0, 3] = omega[
        0, 3, 1
    ] = w2
    omega[0, 1, 3] = -w2

    omega[2, 2, 1] = omega[3, 3, 1] = w3
    omega[2, 1, 2] = omega[3, 1, 3] = -w3

    omega[0, 2, 0] = omega[0, 0, 2] = omega[1, 2, 1] = w4
    omega[1, 1, 2] = -w4

    omega[3, 2, 0] = omega[3, 0, 2] = omega[0, 3, 2] = w5
    omega[2, 3, 0] = omega[2, 0, 3] = omega[0, 2, 3] = -w5

    omega[3, 3, 2] = w6
    omega[3, 2, 3] = -w6

    return omega


def kerr_riemann_tensor(r, theta, a, M=1, config="ulll"):
    """
    Define variables

    Components of the Riemann tensor for Kerr Metric
    From Reference Paper, Appendix
    The Configuration is ulll
    """

    rijkl = np.zeros((4, 4, 4, 4))

    X = r**2 - 3 * a**2 * np.cos(theta) ** 2
    Y = 3 * r**2 - a**2 * np.cos(theta) ** 2

    # Definitions
    Delta = a**2 - 2 * M * r + r**2
    scA = (r**2 + a**2) ** 2 - Delta * a**2 * np.sin(theta) ** 2
    omega_k = 2 * M * a * r / scA
    Sigma = r**2 + a**2 * np.cos(theta) ** 2

    rijkl[0, 0, 0, 3] = 2 * M**2 * a * r**2 * X * np.sin(theta) ** 2 / Sigma**4
    rijkl[3, 3, 0, 3] = -rijkl[0, 0, 0, 3]
    rijkl[0, 3, 0, 3] = -rijkl[0, 0, 0, 3] / omega_k
    rijkl[3, 0, 0, 3] = -rijkl[0, 0, 0, 3] / (
        2 * M * a * r / (Delta - a**2 * np.sin(theta) ** 2)
    )

    rijkl[0, 0, 1, 2] = -(
        M**2 * a**2 * r * Y * np.sin(2 * theta) / (Delta * Sigma**3)
    )
    rijkl[3, 3, 1, 2] = -rijkl[0, 0, 1, 2]
    rijkl[0, 3, 1, 2] = -rijkl[0, 0, 1, 2] / omega_k
    rijkl[3, 0, 1, 2] = -rijkl[0, 0, 1, 2] / (
        2 * M * a * r / (Delta - a**2 * np.sin(theta) ** 2)
    )

    rijkl[3, 2, 2, 3] = -(
        M * r * X * (2 * (r**2 + a**2) + a**2 * np.sin(theta) ** 2) / Sigma**3
    )

    rijkl[0, 1, 0, 1] = -rijkl[3, 2, 2, 3] / Delta

    rijkl[0, 2, 0, 2] = -(
        M * r * X * ((r**2 + a**2) + 2 * a**2 * np.sin(theta) ** 2) / Sigma**3
    )

    rijkl[3, 1, 1, 3] = -rijkl[0, 2, 0, 2] / Delta

    rijkl[0, 1, 0, 2] = rijkl[3, 2, 1, 3] = (
        -M
        * a**2
        / (Delta * Sigma**3)
        * Y
        * (3 * (r**2 + a**2) - 2 * M * r)
        * np.sin(theta)
        * np.cos(theta)
    )

    rijkl[0, 2, 0, 1] = rijkl[3, 2, 1, 3] = (
        -M
        * a**2
        / (Delta * Sigma**3)
        * Y
        * (3 * (r**2 + a**2) - 4 * M * r)
        * np.sin(theta)
        * np.cos(theta)
    )

    rijkl[3, 2, 0, 2] = -3 * M * a * r * X / Sigma**3
    rijkl[3, 1, 0, 1] = -rijkl[3, 2, 0, 2] / Delta

    rijkl[0, 2, 2, 3] = rijkl[3, 2, 0, 2] * np.sin(theta) ** 2 * (r**2 + a**2)
    rijkl[0, 1, 1, 3] = -rijkl[0, 2, 2, 3] / Delta

    rijkl[1, 0, 0, 2] = (
        -3 * M * a**2 * Delta / Sigma**4 * Y * np.sin(theta) * np.cos(theta)
    )
    rijkl[2, 0, 0, 1] = rijkl[1, 0, 0, 2] / Delta

    rijkl[1, 0, 1, 3] = (
        M
        * a
        * r
        / Sigma**4
        * X
        * np.sin(theta) ** 2
        * (3 * (r**2 + a**2) - 4 * M * r)
    )
    rijkl[1, 3, 0, 1] = -rijkl[1, 0, 1, 3]

    rijkl[2, 0, 2, 3] = -(
        M
        * a
        * r
        / Sigma**4
        * X
        * np.sin(theta) ** 2
        * (3 * (r**2 + a**2) - 2 * M * r)
    )

    rijkl[2, 3, 0, 2] = -rijkl[2, 0, 2, 3]

    rijkl[1, 0, 2, 3] = (
        -M
        * a
        * Delta
        / Sigma**4
        * Y
        * np.sin(theta)
        * np.cos(theta)
        * (2 * (r**2 + a**2) + a**2 * np.sin(theta) ** 2)
    )
    rijkl[2, 3, 0, 1] = -rijkl[1, 0, 2, 3] / Delta

    rijkl[1, 3, 0, 2] = (
        M
        * a
        * Delta
        / Sigma**4
        * Y
        * np.sin(theta)
        * np.cos(theta)
        * ((r**2 + a**2) + 2 * a**2 * np.sin(theta) ** 2)
    )

    rijkl[2, 0, 1, 3] = -rijkl[1, 3, 0, 2] / Delta

    rijkl[1, 2, 0, 3] = Delta**2 * rijkl[0, 0, 1, 2] / (2 * M * a * r)
    rijkl[2, 1, 0, 3] = -rijkl[1, 2, 0, 3] / Delta

    rijkl[1, 3, 2, 3] = -(r**2 + a**2) * np.sin(theta) ** 2 * rijkl[1, 0, 0, 2]
    rijkl[2, 3, 1, 3] = rijkl[1, 3, 2, 3] / Delta

    rijkl[1, 2, 1, 2] = -M * r * X / Sigma**2
    rijkl[2, 1, 1, 2] = -rijkl[1, 2, 1, 2] / Delta

    rijkl[0, 1, 2, 3] = (
        -M
        * a
        * Y
        * (2 * (r**2 + a**2) ** 2 + Delta * a**2 * np.sin(theta) ** 2)
        * np.sin(theta)
        * np.cos(theta)
        / (Delta * Sigma**3)
    )
    rijkl[0, 2, 1, 3] = (
        -M
        * a
        * Y
        * ((r**2 + a**2) ** 2 + 2 * Delta * a**2 * np.sin(theta) ** 2)
        * np.sin(theta)
        * np.cos(theta)
        / (Delta * Sigma**3)
    )

    rijkl[3, 1, 0, 2] = (
        -M
        * a
        * Y
        * (Delta + 2 * a**2 * np.sin(theta) ** 2)
        / (np.tan(theta) * Delta * Sigma**3)
    )
    rijkl[3, 2, 0, 1] = (
        -M
        * a
        * Y
        * (2 * Delta + a**2 * np.sin(theta) ** 2)
        / (np.tan(theta) * Delta * Sigma**3)
    )

    rijkl[1, 0, 0, 1] = (
        M * r * X * (2 * Delta + a**2 * np.sin(theta) ** 2) / Sigma**4
    )

    rijkl[2, 0, 0, 2] = -(
        M * r * X * (Delta + 2 * a**2 * np.sin(theta) ** 2) / Sigma**4
    )

    rijkl[1, 3, 1, 3] = (
        -M
        * r
        * X
        * ((r**2 + a**2) ** 2 + 2 * Delta * a**2 * np.sin(theta) ** 2)
        * np.sin(theta) ** 2
        / Sigma**4
    )
    rijkl[2, 3, 2, 3] = (
        M
        * r
        * X
        * (2 * (r**2 + a**2) ** 2 + Delta * a**2 * np.sin(theta) ** 2)
        * np.sin(theta) ** 2
        / Sigma**4
    )

    if config == "ulll":
        return antisymmetrize(arr=rijkl)
    elif config == "llll":
        rulll = antisymmetrize(arr=rijkl)
        return np.einsum(
            "ij,jklm->iklm", metric_tensor(r=r, theta=theta, a=a, M=M), rulll
        )

    else:
        # #!Config: lluu
        gkinv = np.linalg.inv(metric_tensor(r=r, theta=theta, a=a, M=M))

        rllll = antisymmetrize(
            arr=np.einsum(
                "ij,jklm->iklm", metric_tensor(r=r, theta=theta, a=a, M=M), rijkl
            )
        )
        return np.einsum("kjlm,al,bm->kjab", rllll, gkinv, gkinv)


def metric_D(r, theta, a, M=1.0):
    dgkt = np.zeros((4, 4))
    dgkp = np.zeros((4, 4))

    dgkt[1, 0] = (M * (a**2 - 2 * r**2 + a**2 * np.cos(2 * theta))) / (
        r**2 + a**2 * np.cos(theta) ** 2
    ) ** 2

    dgkt[1, 3] = -(
        (2 * a * M * (-(r**2) + a**2 * np.cos(theta) ** 2) * np.sin(theta) ** 2)
        / (r**2 + a**2 * np.cos(theta) ** 2) ** 2
    )

    dgkt[2, 0] = (2 * a**2 * M * r * np.sin(2 * theta)) / (
        r**2 + a**2 * np.cos(theta) ** 2
    ) ** 2

    dgkt[2, 3] = -(
        (2 * a * M * r * (a**2 + r**2) * np.sin(2 * theta))
        / (r**2 + a**2 * np.cos(theta) ** 2) ** 2
    )

    dgkp[1, 0] = -(
        (2 * a * M * (-(r**2) + a**2 * np.cos(theta) ** 2) * np.sin(theta) ** 2)
        / (r**2 + a**2 * np.cos(theta) ** 2) ** 2
    )

    dgkp[1, 3] = (
        2
        * np.sin(theta) ** 2
        * (
            np.cos(theta) ** 2
            * (
                2 * a**2 * r * (a**2 + r**2)
                + a**4 * (M - r) * np.sin(theta) ** 2
            )
            + r * (-(a**4) + r**4 + a**2 * (a**2 - M * r) * np.sin(theta) ** 2)
        )
    ) / (r**2 + a**2 * np.cos(theta) ** 2) ** 2

    # Compute the expressions
    dgkp[2,0] = -(
        (2 * a * M * r * (a**2 + r**2) * np.sin(2 * theta))
        / (r**2 + a**2 * np.cos(theta) ** 2) ** 2
    )
    dgkp[2,3] = (
        (
            3 * a**6
            + 10 * a**4 * M * r
            + 11 * a**4 * r**2
            + 16 * a**2 * M * r**3
            + 16 * a**2 * r**4
            + 8 * r**6
            + 4* a**2* (a**2 + 2 * r**2)* (a**2 + r * (-2 * M + r))
            * np.cos(2 * theta)
            + a**4 * (a**2 - 2 * M * r + r**2) * np.cos(4 * theta)
        )* np.sin(2 * theta)
    ) / (8 * (r**2 + a**2 * np.cos(theta) ** 2) ** 2)

    return dgkt.T, dgkp.T


def spin_matrix(sa,pb,gk,epsilon):
    """
    Resultant Config: aa
    Input Config: u,u,uu,llll
    In Tetrad Basis or Coordinate basis
    Output Config: uu
    """
    gkinv=np.linalg.inv(gk)
    return np.einsum("ab,cd,bdkl,l,k->ac",gkinv,gkinv,epsilon,sa,pb)



def get_dual(r,theta,a,epsilon,config='lluu',M=1):
    """
    Eq 5.28 Format,Output: 'llll'
    """
    rijkl=kerr_riemann_tensor(r=r,theta=theta,a=a,M=M,config=config)
    return np.einsum("ijkl,klmn->ijmn",rijkl,epsilon)*0.5



def get_dual_dual(r,theta,a,gk,epsilon,M=1):
    """
    Eq 5.31 Format,Output: 'llll'
    """
    gkinv=np.linalg.inv(gk)
    drijkl=get_dual(r,theta,a,epsilon,M=M)
    # #!Config: llll

    sdrijkl=np.einsum('ij,kl,jlps->ikps',gkinv,gkinv,drijkl)
    # #!Config: uull

    return np.einsum("uvab,abps->uvps",epsilon,sdrijkl)*0.5


def momentas(rsol,a,M):
    """
    Function returns the momentas given the current 4-position vectors in covariant fashion
    """
    
    guv=metric_tensor(r=rsol[1],theta=rsol[2],a=a,M=M)
    p=rsol[4:8]
    return np.einsum('ij,j->i',guv,p)


def to_cartesian(rsol):
    r, t, p = rsol[:, 1], rsol[:, 2], rsol[:, 3]

    x = r * np.sin(t) * np.cos(p)
    y = r * np.sin(t) * np.sin(p)
    z = r * np.cos(t)

    return np.concatenate((x[:, None], y[:, None], z[:, None]), axis=1)


# +
#-------------New Code---------

# +
M = 1.0
m = 1e-6

S = -1.0*m*M
E = 0.9328 * m
Jz = 2.8 * m * M

r0 = 6.0 * M
theta0 = np.pi / 2
phi0 = 0.0
a0 = 0.8 * M

S1 = 0.0 
S3 = 0.0 

P2=0.0
# -

gk=metric_tensor(r=r0,theta=theta0,a=a0,M=M)

gkinv=np.linalg.inv(gk)


def rhs22(vmu):
    return vmu


# @snoop
def rhs23(y, vmu, sp, ps, gk, epsilon, a, M=1.0):
    _, r, theta, _ = y
    gkinv = np.linalg.inv(gk)
    drijkl = get_dual(r=r, theta=theta, a=a, epsilon=epsilon, M=M)
    sdrijkl = np.einsum("ua,avps->uvps", gkinv, drijkl)
    cs = kerr_christoffel(r=r, theta=theta, a=a, M=M)

    term1 = np.einsum("uvps,v,p,s->u", sdrijkl, vmu, sp, ps) / m

    # term2 = np.einsum("usp,p,s->u", cs, ps, ps)/m
    term2 = np.einsum("usp,p,s->u", cs, vmu, ps)

    return term1 - term2


def rhs24(y, vmu, sp, ps, epsilon, a,M=1.):
    _, r, theta, _ = y
    drijkl = get_dual(r=r, theta=theta, a=a, epsilon=epsilon,M=M)
    cs = kerr_christoffel(r=r, theta=theta, a=a,M=M)

    term1 = np.einsum("u,vpsg,v,p,s,g->u", ps, drijkl, sp, vmu, sp, ps) / (m**3)

    term2 = np.einsum("usp,p,s->u", cs, vmu, sp)

    return term1 - term2


def rhs27(y, umu, su, gk, eps, a,m, config="l", M=1):
    _, r, theta, _ = y
    gkinv = np.linalg.inv(gk)
    ddrijkl = get_dual_dual(r=r, theta=theta, a=a, gk=gk, epsilon=eps, M=M)

    if config == "l":
        # #!Umu is u_mu
        term2 = np.einsum("uvps,v,p,sj,j->u", ddrijkl, su, su,gkinv, umu) / (m**2)

        return umu + term2

    else:
        # #!Umu is u^mu
        sddrijkl = np.einsum("ab,bvps->avps", gkinv, ddrijkl)
        term2 = np.einsum("uvps,v,p,s->u", sddrijkl, su, su, umu) / (m**2)

        return umu + term2


def mpd(t, y,a,m,M):
    """
    Solve the MPD Equations Given the Initial Conditions
    """

    dy=np.zeros_like(y)
    ps=y[4:8]
    su=y[8:]

    gk=metric_tensor(r=y[1],theta=y[2],a=a,M=M)

    eps=levi_cevita_tensor(dim=4)

    vmu=rhs27(y=y[:4],umu=ps/m,su=su,gk=gk,eps=eps,a=a,m=m,config='u',M=M)

    dy[:4]=vmu

    dy[4:8]=rhs23(y=y[:4],vmu=vmu,sp=su,ps=ps,gk=gk,epsilon=eps,a=a,M=M)
    dy[8:]=rhs24(y=y[:4],vmu=vmu,sp=su,ps=ps,epsilon=eps,a=a,M=M)

    return dy
    


p0,p3,p1=smp.symbols("p^0,p^3,p^1",real=True)
# #!In tetrad basis p^i

pi = np.array([p0, p1, P2, p3])
# #!In tetrad basis p^i

# ## The equations are in Tetrad Basis

# Definitions
Delta = a0**2 - 2 * M * r0 + r0**2
scA = (r0**2 + a0**2) ** 2 - Delta * a0**2 * np.sin(theta0) ** 2
omega_k = 2 * M * a0 * r0 / scA
Sigma = r0**2 + a0**2 * np.cos(theta0) ** 2

Suin = np.array([0.0, S1, -S, S3])
# #!In tetrad basis, S^p

epsl=levi_cevita_tensor(dim=4)
# #!Config: llll, Tetrad Basis

eta = np.array([[-1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]).astype(float)
# #!Config: uu or ll

Sab=spin_matrix(sa=Suin,pb=pi,gk=eta,epsilon=epsl)
# #!Config: uu, In tetrad Basis

Sab

eq210a=np.sqrt(Delta/Sigma)*pi[0]+a0/np.sqrt(Sigma) *np.sin(theta0)*pi[3]+ M/Sigma**2 *(r0**2-a0**2*np.cos(theta0)**2)*Sab[1,0]+ 2*M*r0*a0*np.cos(theta0)/Sigma**2 *Sab[2,3] - E

eq210a

eq210b=a0*np.sin(theta0)**2*np.sqrt(Delta/Sigma)*pi[0]+(r0**2+a0**2)*np.sin(theta0)/np.sqrt(Sigma)*pi[3]+a0*((r0-M)*Sigma+2*M*r0**2)/Sigma**2 *np.sin(theta0)**2*Sab[1,0]+a0*np.sqrt(Delta)*np.sin(theta0)*np.cos(theta0)*Sab[2,0]/Sigma +\
    ((r0**2+a0**2)**2-a0**2*Delta*np.sin(theta0)**2)*np.cos(theta0)*Sab[2,3]/Sigma**2 + r0*np.sqrt(Delta)*np.sin(theta0)*Sab[1,3]/Sigma- Jz

eq210b

eq533=get_norm(T=pi,g=eta)+m**2

eq533

smp.nonlinsolve([eq533,eq210a,eq210b],[p0,p1,p3])

eq210a,eq210b,eq533

smp.nsolve([eq533, eq210a, eq210b], [p0, p1, p3],[1e-5,1e-5,1e-5])

smp.solve([eq210a, eq210b, eq533], [p0, p1, p3])[0]

# +
from logging import raiseExceptions


try:
    p0,p1,p3=smp.solve([eq210a,eq210b,eq533],[p0,p1,p3])[0]
except:
    raiseExceptions("Program failed to find the solution")
# -

get_norm(T=np.array([p0,p1,0.,p3]),g=eta)

p0,p1,p3


# ## Convert from Tetrad to Coordinate bases

# @snoop
def co_tat(r, theta, a, vec, method="to"):
    """
    This is a helper function to convert from tetrad basis to
    coordinate basis and vica-verse.

    The function needs the information on r, theta and a.

    The method takes argument 'to' or 'from' taking coordinate_tetrad or
    coordinate_from_tetrad

    The function takes an array vec. This is a upper indexed four-vector in
    either of basis.

    In tetrad basis, the latin indexed are lowered by minkowski metric while
    the greek index are lowered by the spacetime metric.

    So, if a vec with configuration 'u' is passed in latin index fashion, the covarint
    vector in coordinate basis is e^i_mu vec^mu= vec^i

    Next, to convert from tetrad basis to coordinate basis, we need the inverse of tetrad basis vectors
    It is given by e_i^mu. Hence, if a vec with configuration 'u' is passed as vec^i, the resultant coordinate basis
    vector is e_i^mu vec^i=vec^mu
    """
    # Definitions
    Delta = a**2 - 2 * M * r + r**2
    scA = (r**2 + a**2) ** 2 - Delta * a**2 * np.sin(theta) ** 2
    omega_k = 2 * M * a * r / scA
    Sigma = r**2 + a**2 * np.cos(theta) ** 2

    e0 = np.array(
        [
            [
                np.sqrt(Delta / Sigma),
                0.0,
                0.0,
                -a * np.sin(theta) ** 2 * np.sqrt(Delta / Sigma),
            ]
        ]
    ).T

    e1 = np.array([[0.0, np.sqrt(Sigma / Delta), 0.0, 0.0]]).T
    e2 = np.array([[0.0, 0.0, np.sqrt(Sigma), 0.0]]).T
    e3 = np.array(
        [
            [
                -a * np.sin(theta) / np.sqrt(Sigma),
                0.0,
                0.0,
                (r**2 + a**2) * np.sin(theta) / np.sqrt(Sigma),
            ]
        ]
    ).T

    tet = np.hstack((e0, e1, e2, e3))

    # #! Now for tetrad to Coordinate basis
    tet_inv = np.linalg.inv(tet)

    # coo_vec = np.einsum("ij,j->i", tet_inv, vec)

    if method == "to":
        if vec.ndim == 1:
            return np.einsum("ij,i->j", tet, vec)
        else:
            return np.einsum("ij,kl,ik->jl", tet, tet, vec)

    elif method == "from":
        if vec.ndim == 1:
            return np.einsum("ij,i->j", tet_inv, vec)
        else:
            return np.einsum("ij,kl,ik->jl", tet_inv, tet_inv, vec)


pi = np.array([p0, p1, 0.0, p3]).astype(float)

pi

p_coordu=co_tat(r=r0, theta=theta0, a=a0, vec=pi,method="from")

get_norm(T=p_coordu,g=gk)

sinl=np.array([0.,0.,-S,0.])
# #!Spin s^i in tetrad basis

s_coordu=co_tat(r=r0,theta=theta0,a=a0,vec=sinl,method='from')
# #!Spin s^mu in coordinate basis

get_norm(T=s_coordu,g=gk)

# +
vmu0=rhs27(y=[0.,r0,theta0,phi0],umu=p_coordu/m,su=s_coordu,gk=gk,eps=levi_cevita_tensor(dim=4),a=a0,m=m,config='u',M=M)

# #! Get v^mu according to equation 2.27  in paper
# -

vmu0

N=1/(1-(M*S**2)*(1+3*p3**2/m**2)/r0**3)

# +
v_tet=np.array([N * (1 - M * S**2 / r0**3) * p0 / m, N * (1 - M * S**2 / r0**3) * p1 / m,0., N * (
    1 + 2*M * S**2 / r0**3
) * p3 / m]).astype(float)

# #!Eq: 2.15-17 in paper
# -

v_tet

co_tat(r=r0, theta=theta0, a=a0,vec=vmu0,method="to") #!v^i


def get_constants(rsol, eta, gk, a=a0,M=1):
    _, r, theta, _ = rsol[:4]
    pmuu = rsol[4:8]
    smuu = rsol[8:]

    # pmul = np.einsum("ij,j->i", gk, pmuu) #!'l', Coordinate
    # smul = np.einsum("ij,j->i", gk, smuu) #! 'l', Coordinate

    pi = co_tat(r=r, theta=theta, a=a, vec=pmuu, method="to") #! 'u', tetrad
    si = co_tat(r=r, theta=theta, a=a, vec=smuu, method="to") #! 'u', tetrad

    # pi = np.einsum("ij,j->i", eta, pil) #! 'u', tetrad
    # si = np.einsum("ij,j->i", eta, sil) #! 'u', tetrad

    epsilon=levi_cevita_tensor(dim=4) #!'llll', tetrad
    Sab=spin_matrix(sa=si,pb=pi,gk=eta,epsilon=epsilon) 

    # Definitions
    Delta = a**2 - 2 * M * r + r**2
    scA = (r**2 + a**2) ** 2 - Delta * a**2 * np.sin(theta) ** 2
    omega_k = 2 * M * a * r / scA
    Sigma = r**2 + a**2 * np.cos(theta) ** 2

    E = (
        np.sqrt(Delta / Sigma) * pi[0]
        + a / np.sqrt(Sigma) * np.sin(theta) * pi[3]
        + M / Sigma**2 * (r**2 - a**2 * np.cos(theta) ** 2) * Sab[1, 0]
        + 2 * M * r * a * np.cos(theta) / Sigma**2 * Sab[2, 3]
        
    )

    Jz = (
        a * np.sin(theta) ** 2 * np.sqrt(Delta / Sigma) * pi[0]
        + (r**2 + a**2) * np.sin(theta) / np.sqrt(Sigma) * pi[3]
        + a* ((r - M) * Sigma + 2 * M * r**2)/ Sigma**2
        * np.sin(theta) ** 2
        * Sab[1, 0]
        + a * np.sqrt(Delta) * np.sin(theta) * np.cos(theta) * Sab[2, 0] / Sigma
        + ((r**2 + a**2) ** 2 - a**2 * Delta * np.sin(theta) ** 2)
        * np.cos(theta) * Sab[2, 3] / Sigma**2
        + r * np.sqrt(Delta) * np.sin(theta) * Sab[1, 3] / Sigma
    )

    mu=get_norm(T=pi,g=eta)
    Smag=get_norm(T=si,g=eta)

    return E,Jz,mu,Smag


# +
#In Coordinate basis

# +

# #!The initial conditions
# #!We pass r0^mu,p^mu and S^mu
y0=np.array([0.,r0,theta0,phi0,p_coordu[0],p_coordu[1],p_coordu[2],p_coordu[3],s_coordu[0],s_coordu[1],s_coordu[2],s_coordu[3]]).astype(float)
# -

y0

mpd(0.,y0,a0,m,M)


# ## Solution of Diffferential Equation

# +
# @snoop(depth=2)
def integrate_BulirschStoeir(F, x, y, xStop, tol):
    def midpoint(F, x, y, xStop, nSteps):
        h = (xStop - x) / nSteps

        y0 = y
        y1 = y0 + h * F(x, y0, a=a0, m=m, M=M)
        for i in range(nSteps - 1):
            x = x + h
            y2 = y0 + 2.0 * h * F(x, y1, a=a0, m=m, M=M)
            y0 = y1
            y1 = y2

        return 0.5 * (y1 + y0 + h * F(x, y2, a=a0,m=m,M=M))

    def richardson(r, k):
        for j in range(k - 1, 0, -1):
            const = (k / (k - 1.0)) ** (2.0 * (k - j))
            r[j] = (const * r[j + 1] - r[j]) / (const - 1.0)
        return

    kMax = 51  ## Indentation problem from here
    n = len(y)
    r = np.zeros((kMax, n), dtype=float)

    # Start with two integration steps
    nSteps = 2
    r[1] = midpoint(F, x, y, xStop, nSteps)
    r_old = r[1].copy()

    # Increase the number of integration points by 2 and refine result by Richardson extrapolation
    for k in range(2, kMax):
        nSteps = 2 * k
        r[k] = midpoint(F, x, y, xStop, nSteps)
        richardson(r, k)

        # Compute RMS change in solution
        e = np.sqrt(sum((r[1] - r_old) ** 2) / n)

        # Check for convergence
        if e < tol:
            return r[1]
        r_old = r[1].copy()
    print("Midpoint method did not converge")


# Bulirsch-Stoer Algorithm:-

""" X, Y = bulStoer(F, x, y, xStop, H, tol=1.0e-6).
    Simplified Bulirsch-Stoer method for solving the
    initial value problem {y}’ = {F(x,{y})}, where {y} = {y[0],y[1],...y[n-1]}
    x, y = initial conditions
    xStop = terminal value of x
    H = increment of x at which results are stored
    F = user-supplied function that returns the array F(x,y) = {y’[0],y’[1],...,y’[n-1]} """

# from numpy import array


def bulStoer(F, x, y, xStop, H, tol=1.0e-6):
    X = []
    Y = []
    X.append(x)
    Y.append(y)
    while x < xStop:
        H = min(H, xStop - x)
        y = integrate_BulirschStoeir(F=F, x=x, y=y, xStop=x + H, tol=tol)  # Midpoint method
        x = x + H
        X.append(x)
        Y.append(y)
    return np.array(X), np.array(Y)


def printSoln(X, Y, freq):
    def printHead(n):
        print("\n        x  ", end=" ")  ## end=" " here
        for i in range(n):
            print("      y[", i, "] ", end=" ")  ## end=" " here
        print()

    def printLine(x, y, n):
        print("{:13.4e}".format(x), end=" ")  ## end=" " here
        for i in range(n):
            print("{:13.4e}".format(y[i]), end=" ")  ## end=" " here
        print()

    m = Y.shape[0]
    try:
        n = Y.shape[1]
    except Exception:
        n = 1
    if freq == 0:
        freq = m
    printHead(n)
    for i in range(0, m, freq):
        printLine(X[i], Y[i], n)
    if i != m - 1:
        printLine(X[m - 1], Y[m - 1], n)


# +
# rsol2=bulStoer(mpd,0.0,y0,20.,0.01,1e-13)[1]
# -

sol = solve_ivp(
    fun=mpd,
    t_span=[0.0, 20.0],
    y0=y0,
    method="LSODA",
    t_eval=np.arange(0.0, 20.0, 1e-3),
    args=(a0,m,M),
    rtol=1e-10,
    atol=1e-8,
)

rsol=sol.y.T

# +
#Check for the constants
# -

E,Jz

i = -1
get_constants(rsol=rsol[i], eta=eta, gk=gk, a=a0, M=M)[0] - E, get_constants(rsol=rsol[i], eta=eta, gk=gk, a=a0, M=M)[
    1
] - Jz, get_constants(rsol=rsol[i], eta=eta, gk=gk, a=a0, M=M)[2] + m**2, get_constants(rsol=
    rsol[i], eta=eta, gk=gk, a=a0, M=M
)[
    3
] - S**2

crsol=to_cartesian(rsol=rsol)

np.average(crsol[:,2])

# +
fig = go.Figure(
    data=go.Scatter3d(
        x=crsol[:, 0],
        y=crsol[:, 1],
        z=crsol[:, 2],
        marker=dict(
            size=4,
            color=crsol[:, 2],
            colorscale="Viridis",
        ),
        line=dict(color="darkblue", width=2),
    )
)

fig.update_layout(
    title=f"Particle motion with the initial conditions E={E:.3E}, Jz={Jz:.3E} and S={S:.3E}",
    width=800,
    height=700,
    autosize=False,
    scene=dict(
        camera=dict(
            up=dict(x=0, y=0, z=1),
            eye=dict(
                x=0,
                y=1.0707,
                z=1,
            ),
        ),
        aspectratio=dict(x=1, y=1, z=0.7),
        aspectmode="manual",
    ),
)

# fig.update_traces(line_colorbar_exponentformat="E", selector=dict(type='parcoords'))

fig.show()
# -

# ## Let us now generalize the model by adding one more constraints

# +
M = 1.0
m = 1e-6

S = 1.0 * m * M
E = 0.9328* m
Jz = 2.8 * m * M

r0 = 6.0 * M
theta0 = np.pi / 2
phi0 = 0.0
a0 = 0.8*M
# -

P1=-1e-3*m

S1=1e-5*m
S2=1e-3*m

p0,p2,p3=smp.symbols("p0 p2 p3",real=True)
# #!Tetrad Components, 'u'

s0,s3=smp.symbols("s0 s3",real=True)
# #!Tetrad Components, 'u'

# +
gpi=np.array([p0,P1,p2,p3])
gsi=np.array([s0,S1,S2,s3])

# #!Config='u'
# -

# Definitions
Delta = a0**2 - 2 * M * r0 + r0**2
scA = (r0**2 + a0**2) ** 2 - Delta * a0**2 * np.sin(theta0) ** 2
omega_k = 2 * M * a0 * r0 / scA
Sigma = r0**2 + a0**2 * np.cos(theta0) ** 2

Sab = spin_matrix(sa=gsi, pb=gpi, gk=eta, epsilon=epsl)
# #!Config: uu, In tetrad Basis

Sab

eq210a = (
    np.sqrt(Delta / Sigma) * gpi[0]
    + a0 / np.sqrt(Sigma) * np.sin(theta0) * gpi[3]
    + M / Sigma**2 * (r0**2 - a0**2 * np.cos(theta0) ** 2) * Sab[1, 0]
    + 2 * M * r0 * a0 * np.cos(theta0) / Sigma**2 * Sab[2, 3]
    - E
)

eq210b = (
    a0 * np.sin(theta0) ** 2 * np.sqrt(Delta / Sigma) * gpi[0]
    + (r0**2 + a0**2) * np.sin(theta0) / np.sqrt(Sigma) * gpi[3]
    + a0
    * ((r0 - M) * Sigma + 2 * M * r0**2)
    / Sigma**2
    * np.sin(theta0) ** 2
    * Sab[1, 0]
    + a0 * np.sqrt(Delta) * np.sin(theta0) * np.cos(theta0) * Sab[2, 0] / Sigma
    + ((r0**2 + a0**2) ** 2 - a0**2 * Delta * np.sin(theta0) ** 2)
    * np.cos(theta0)
    * Sab[2, 3]
    / Sigma**2
    + r0 * np.sqrt(Delta) * np.sin(theta0) * Sab[1, 3] / Sigma
    - Jz
)

eq210b

eq210a = eq210a.xreplace(
    dict([(n, 0) for n in eq210a.atoms(smp.Float) if abs(n) < 1e-15])
)

eq210b=eq210b.xreplace(dict([(n, 0) for n in eq210b.atoms(smp.Float) if abs(n) < 1e-15]))

eq533 = get_norm(T=gpi, g=eta) + m**2

eq229=np.einsum("ij,j,i->",eta,gpi,gsi)

eq228=np.einsum("ij,j,i->",eta,gsi,gsi)-S**2

eq533

eq210a

eq228

eq229

eq210b


# +
# @snoop(depth=1,watch='x')
def solve_eq(p_0, p_2, p_3, s_0, s_3):
    f = smp.lambdify([p0, p2, p3, s0, s3], [eq210a, eq210b, eq533, eq228, eq229])
    return f(p_0, p_2, p_3, s_0, s_3)


# f = smp.lambdify([p0, p2, p3, s0, s3], [eq210a, eq210b, eq533, eq228, eq229])

# +
# @snoop(depth=1,watch='x')
def solve_eq_sp(p):
    f = smp.lambdify([p0, p2, p3, s0, s3], [eq210a, eq210b, eq533, eq228, eq229])
    return f(p[0], p[1], p[2], p[3], p[4])


# f = smp.lambdify([p0, p2, p3, s0, s3], [eq210a, eq210b, eq533, eq228, eq229])
# -

def jac_nl_symb():
    X = smp.Matrix([eq210a, eq210b, eq533, eq228, eq229])
    Y = smp.Matrix([p0, p2, p3, s0, s3])

    return X.jacobian(Y)


jac=jac_nl_symb()


def jac_nl(p_0, p_2, p_3, s_0, s_3):
    f = smp.lambdify([p0, p2, p3, s0, s3], jac)

    return f(p_0, p_2, p_3, s_0, s_3)



def jac_nl_sp(p):
    f = smp.lambdify([p0, p2, p3, s0, s3], jac)

    return f(p[0], p[1], p[2], p[3], p[4])


solve_eq_sp([0.0, 0.0, 0.0, 1e-6, 1e-6])

roots=mpmath.findroot(solve_eq, (1e-6, 1e-6, 1e-6, 1e-6, 1e-6),solver='mdnewton' ,tol=1e-27, J=jac_nl,verify=True)

roots=np.array(roots.tolist(),dtype=np.float64)[:,0]

roots

solve_eq(*roots)

# +
# roots_ls = sp.optimize.least_squares(solve_eq_sp, x0=[1e-6, 1e-6, 1e-6, 1e-6, 1e-6],jac=jac_nl_sp,ftol=1e-15,loss='cauchy')
# -

roots_sp = sp.optimize.root(solve_eq_sp, [1e-6, 1e-6, 1e-6, 1e-6, 1e-6],jac=jac_nl_sp, tol=1e-26, options={'maxfev':10000,'xtol':1e-20})

solve_eq_sp(roots_sp.x)

p0, p2, p3, s0, s3=roots

ngpi = np.array([p0, P1, p2, p3]).astype(float)
ngsi = np.array([s0, S1, S2, s3]).astype(float)

ngpi,ngsi

get_norm(ngpi,eta)

p_coordu = co_tat(r=r0, theta=theta0, a=a0, vec=ngpi, method="from")
#! "u"

get_norm(T=p_coordu,g=gk)

get_norm(ngsi,eta)

s_coordu = co_tat(r=r0, theta=theta0, a=a0, vec=ngsi, method="from")
#! "u"

get_norm(s_coordu, gk)

# +
vmu0 = rhs27(
    y=[0.0, r0, theta0, phi0],
    umu=p_coordu / m,
    su=s_coordu,
    gk=gk,
    eps=levi_cevita_tensor(dim=4),
    a=a0,
    m=m,
    config="u",
    M=M,
)

# #! Get v^mu according to equation 2.27  in paper
# -

vmu0

# #!The initial conditions
# #!We pass r0^mu,p^mu and S^mu
y0 = np.array(
    [
        0.0,
        r0,
        theta0,
        phi0,
        p_coordu[0],
        p_coordu[1],
        p_coordu[2],
        p_coordu[3],
        s_coordu[0],
        s_coordu[1],
        s_coordu[2],
        s_coordu[3],
    ]
).astype(float)

y0


def poincare(t,y,*args):
    # return y[3]%np.pi -1.
    if y[0]>10:
        return y[2]-np.pi/2
    else:
        return y[1]


sol = solve_ivp(
    fun=mpd,
    t_span=[0.0, 100000.0],
    y0=y0,
    method="LSODA",
    t_eval=np.arange(0.0, 50000.0, 1.0),
    events=poincare,
    args=(a0, m, M),
    rtol=1e-10,
    atol=1e-8,
)

poi_points=sol.y_events[0]

rsol = sol.y.T

i = -1
get_constants(rsol=rsol[i], eta=eta, gk=gk, a=a0, M=M)[0] - E, get_constants(
    rsol=rsol[i], eta=eta, gk=gk, a=a0, M=M
)[1] - Jz, get_constants(rsol=rsol[i], eta=eta, gk=gk, a=a0, M=M)[
    2
] + m**2, get_constants(
    rsol=rsol[i], eta=eta, gk=gk, a=a0, M=M
)[
    3
] - S**2

crsol = to_cartesian(rsol=rsol)

np.average(crsol[:, 2])

crsol.shape

# +
fig = go.Figure(
    data=go.Scatter3d(
        x=crsol[::2, 0],
        y=crsol[::2, 1],
        z=crsol[::2, 2],
        marker=dict(
            size=4,
            color=crsol[::2, 2],
            colorscale="Viridis",
        ),
        line=dict(color="darkblue", width=2),
    )
)

fig.update_layout(
    title=f"Particle motion with the initial conditions E={E:.3E}, Jz={Jz:.3E} and S={S:.3E}",
    width=800,
    height=700,
    autosize=False,
    scene=dict(
        camera=dict(
            up=dict(x=0, y=0, z=1),
            eye=dict(
                x=0,
                y=1.0707,
                z=1,
            ),
        ),
        aspectratio=dict(x=1, y=1, z=0.7),
        aspectmode="manual",
    ),
)

# fig.update_traces(line_colorbar_exponentformat="E", selector=dict(type='parcoords'))

fig.show()

# +
# fig.write_html("chaoticfig.html")

# +
#Poicare Plot
# -

plt.plot(poi_points[:,2])

# +
# fig = go.Figure(
#     data=go.Scatter(
#         x=poi_points[:, 1],
#         y=poi_points[:, 5],
#         mode="markers",
#         marker_color="rgba(0, 0, 0, 1)",
#         marker_size=1,
#     )
# )
# fig.update_layout(
#     autosize=False,
#     width=500,
#     height=500,
#     margin=go.layout.Margin(
#         l=10.2,  # left margin
#         r=10.2,  # right margin
#         b=10.2,  # bottom margin
#         t=10.2,  # top margin
#     ),
# )

# fig.show()
# -

plt.figure(dpi=1200)
plt.scatter(poi_points[:,1],poi_points[:,5],s=0.1,c='black')
plt.xlim(5.90,6.65)
# plt.ylim(-4,-3)

def Smetric_tensor():
    """
    Define the metric tensor function.

    Parameters:
        a (float): Parameter 'a' in the metric.
        r (float): Parameter 'r' in the metric.
        theta (float): Parameter 'theta' in the metric.

    Returns:
        np.ndarray: The metric tensor at the given values of a, r, and theta.
        Configuration: ll
    """

    r, theta, a, M= smp.symbols("r theta a M")
    g = np.array(
        [
            [
                (-(a**2) + 2 * M * r - r**2 + a**2 * smp.sin(theta) ** 2)
                / (r**2 + a**2 * smp.cos(theta) ** 2),
                0,
                0,
                -(
                    (2 * a * M * r * smp.sin(theta) ** 2)
                    / (r**2 + a**2 * smp.cos(theta) ** 2)
                ),
            ],
            [
                0,
                (r**2 + a**2 * smp.cos(theta) ** 2) / (a**2 - 2 * M * r + r**2),
                0,
                0,
            ],
            [0, 0, r**2 + a**2 * smp.cos(theta) ** 2, 0],
            [
                -(
                    (2 * a * M * r * smp.sin(theta) ** 2)
                    / (r**2 + a**2 * smp.cos(theta) ** 2)
                ),
                0,
                0,
                (
                    smp.sin(theta) ** 2
                    * (
                        (a**2 + r**2) ** 2
                        - a**2 * (a**2 - 2 * M * r + r**2) * smp.sin(theta) ** 2
                    )
                )
                / (r**2 + a**2 * smp.cos(theta) ** 2),
            ],
        ]
    )

    return g


print(Smetric_tensor())


def Skerr_christoffel():
    """
    Function to give the christoffel symbols for kerr metric.
    The christoffel symbols are given as \Gamma ^i _{jk}

    From Reference Paper, Appendix
    Config: ull
    """
    r, theta, a, M = smp.symbols("r theta a M")

    cs = np.zeros((4, 4, 4),dtype=object)

    # Definitions
    Delta = a**2 - 2 * M * r + r**2
    scA = (r**2 + a**2) ** 2 - Delta * a**2 * smp.sin(theta) ** 2
    omega_k = 2 * M * a * r / scA
    Sigma = r**2 + a**2 * smp.cos(theta) ** 2

    cs[3, 0, 1] = M * (2 * r**2 - Sigma) / (Delta * Sigma**2) * a
    cs[0, 0, 1] = cs[3, 0, 1] * (r**2 + a**2) / a

    cs[3, 0, 2] = -2 * M * a * r / (smp.tan(theta) * Sigma**2)
    cs[0, 0, 2] = a * smp.sin(theta) ** 2 * cs[3, 0, 2]

    cs[0, 1, 3] = (
        -M
        * a
        * (2 * r**2 * (r**2 + a**2) + Sigma * (r**2 - a**2))
        * smp.sin(theta) ** 2
        / (Delta * Sigma**2)
    )

    cs[3, 1, 3] = (
        r * Sigma * (Sigma - 2 * M * r)
        - M * a**2 * (2 * r**2 - Sigma) * smp.sin(theta) ** 2
    ) / (Delta * Sigma**2)

    cs[0, 2, 3] = M * a**3 * r * smp.sin(theta) ** 2 * smp.sin(2 * theta) / Sigma**2

    cs[3, 2, 3] = (scA - Sigma * a**2 * smp.sin(theta) ** 2) / (
        smp.tan(theta) * Sigma**2
    )

    cs[1, 0, 0] = M * Delta * (2 * r**2 - Sigma) / Sigma**3

    cs[1, 0, 3] = -cs[1, 0, 0] * a * smp.sin(theta) ** 2

    cs[2, 0, 0] = -M * a * r * smp.sin(2 * theta) / Sigma**3 * a
    cs[2, 0, 3] = -cs[2, 0, 0] * (r**2 + a**2) / a

    cs[1, 1, 1] = r / Sigma + (M - r) / Delta

    cs[1, 2, 2] = -r * Delta / Sigma
    cs[2, 1, 2] = -cs[1, 2, 2] / Delta

    cs[1, 1, 2] = cs[2, 2, 2] = -(a**2) * smp.sin(2 * theta) / (2 * Sigma)
    cs[2, 1, 1] = -cs[1, 1, 2] / Delta

    cs[1, 3, 3] = (
        -Delta
        * (r * Sigma**2 - M * a**2 * (2 * r**2 - Sigma) * smp.sin(theta) ** 2)
        * smp.sin(theta) ** 2
        / Sigma**3
    )

    cs[2, 3, 3] = (
        -(Delta * Sigma**2 + 2 * M * r * (r**2 + a**2) ** 2)
        * smp.sin(2 * theta)
        / (2 * Sigma**3)
    )

    return symmetrize(arr=cs)
    # return cs


def Skerr_riemann_tensor(config="ulll"):
    """
    Define variables

    Components of the Riemann tensor for Kerr Metric
    From Reference Paper, Appendix
    The Configuration is ulll
    """

    r, theta, a, M = smp.symbols("r theta a M")

    rijkl = np.zeros((4, 4, 4, 4), dtype=object)

    X = r**2 - 3 * a**2 * smp.cos(theta) ** 2
    Y = 3 * r**2 - a**2 * smp.cos(theta) ** 2

    # Definitions
    Delta = a**2 - 2 * M * r + r**2
    scA = (r**2 + a**2) ** 2 - Delta * a**2 * smp.sin(theta) ** 2
    omega_k = 2 * M * a * r / scA
    Sigma = r**2 + a**2 * smp.cos(theta) ** 2

    rijkl[0, 0, 0, 3] = 2 * M**2 * a * r**2 * X * smp.sin(theta) ** 2 / Sigma**4
    rijkl[3, 3, 0, 3] = -rijkl[0, 0, 0, 3]
    rijkl[0, 3, 0, 3] = -rijkl[0, 0, 0, 3] / omega_k
    rijkl[3, 0, 0, 3] = -rijkl[0, 0, 0, 3] / (
        2 * M * a * r / (Delta - a**2 * smp.sin(theta) ** 2)
    )

    rijkl[0, 0, 1, 2] = -(
        M**2 * a**2 * r * Y * smp.sin(2 * theta) / (Delta * Sigma**3)
    )
    rijkl[3, 3, 1, 2] = -rijkl[0, 0, 1, 2]
    rijkl[0, 3, 1, 2] = -rijkl[0, 0, 1, 2] / omega_k
    rijkl[3, 0, 1, 2] = -rijkl[0, 0, 1, 2] / (
        2 * M * a * r / (Delta - a**2 * smp.sin(theta) ** 2)
    )

    rijkl[3, 2, 2, 3] = -(
        M * r * X * (2 * (r**2 + a**2) + a**2 * smp.sin(theta) ** 2) / Sigma**3
    )

    rijkl[0, 1, 0, 1] = -rijkl[3, 2, 2, 3] / Delta

    rijkl[0, 2, 0, 2] = -(
        M * r * X * ((r**2 + a**2) + 2 * a**2 * smp.sin(theta) ** 2) / Sigma**3
    )

    rijkl[3, 1, 1, 3] = -rijkl[0, 2, 0, 2] / Delta

    rijkl[0, 1, 0, 2] = rijkl[3, 2, 1, 3] = (
        -M
        * a**2
        / (Delta * Sigma**3)
        * Y
        * (3 * (r**2 + a**2) - 2 * M * r)
        * smp.sin(theta)
        * smp.cos(theta)
    )

    rijkl[0, 2, 0, 1] = rijkl[3, 2, 1, 3] = (
        -M
        * a**2
        / (Delta * Sigma**3)
        * Y
        * (3 * (r**2 + a**2) - 4 * M * r)
        * smp.sin(theta)
        * smp.cos(theta)
    )

    rijkl[3, 2, 0, 2] = -3 * M * a * r * X / Sigma**3
    rijkl[3, 1, 0, 1] = -rijkl[3, 2, 0, 2] / Delta

    rijkl[0, 2, 2, 3] = rijkl[3, 2, 0, 2] * smp.sin(theta) ** 2 * (r**2 + a**2)
    rijkl[0, 1, 1, 3] = -rijkl[0, 2, 2, 3] / Delta

    rijkl[1, 0, 0, 2] = (
        -3 * M * a**2 * Delta / Sigma**4 * Y * smp.sin(theta) * smp.cos(theta)
    )
    rijkl[2, 0, 0, 1] = rijkl[1, 0, 0, 2] / Delta

    rijkl[1, 0, 1, 3] = (
        M
        * a
        * r
        / Sigma**4
        * X
        * smp.sin(theta) ** 2
        * (3 * (r**2 + a**2) - 4 * M * r)
    )
    rijkl[1, 3, 0, 1] = -rijkl[1, 0, 1, 3]

    rijkl[2, 0, 2, 3] = -(
        M
        * a
        * r
        / Sigma**4
        * X
        * smp.sin(theta) ** 2
        * (3 * (r**2 + a**2) - 2 * M * r)
    )

    rijkl[2, 3, 0, 2] = -rijkl[2, 0, 2, 3]

    rijkl[1, 0, 2, 3] = (
        -M
        * a
        * Delta
        / Sigma**4
        * Y
        * smp.sin(theta)
        * smp.cos(theta)
        * (2 * (r**2 + a**2) + a**2 * smp.sin(theta) ** 2)
    )
    rijkl[2, 3, 0, 1] = -rijkl[1, 0, 2, 3] / Delta

    rijkl[1, 3, 0, 2] = (
        M
        * a
        * Delta
        / Sigma**4
        * Y
        * smp.sin(theta)
        * smp.cos(theta)
        * ((r**2 + a**2) + 2 * a**2 * smp.sin(theta) ** 2)
    )

    rijkl[2, 0, 1, 3] = -rijkl[1, 3, 0, 2] / Delta

    rijkl[1, 2, 0, 3] = Delta**2 * rijkl[0, 0, 1, 2] / (2 * M * a * r)
    rijkl[2, 1, 0, 3] = -rijkl[1, 2, 0, 3] / Delta

    rijkl[1, 3, 2, 3] = -(r**2 + a**2) * smp.sin(theta) ** 2 * rijkl[1, 0, 0, 2]
    rijkl[2, 3, 1, 3] = rijkl[1, 3, 2, 3] / Delta

    rijkl[1, 2, 1, 2] = -M * r * X / Sigma**2
    rijkl[2, 1, 1, 2] = -rijkl[1, 2, 1, 2] / Delta

    rijkl[0, 1, 2, 3] = (
        -M
        * a
        * Y
        * (2 * (r**2 + a**2) ** 2 + Delta * a**2 * smp.sin(theta) ** 2)
        * smp.sin(theta)
        * smp.cos(theta)
        / (Delta * Sigma**3)
    )
    rijkl[0, 2, 1, 3] = (
        -M
        * a
        * Y
        * ((r**2 + a**2) ** 2 + 2 * Delta * a**2 * smp.sin(theta) ** 2)
        * smp.sin(theta)
        * smp.cos(theta)
        / (Delta * Sigma**3)
    )

    rijkl[3, 1, 0, 2] = (
        -M
        * a
        * Y
        * (Delta + 2 * a**2 * smp.sin(theta) ** 2)
        / (smp.tan(theta) * Delta * Sigma**3)
    )
    rijkl[3, 2, 0, 1] = (
        -M
        * a
        * Y
        * (2 * Delta + a**2 * smp.sin(theta) ** 2)
        / (smp.tan(theta) * Delta * Sigma**3)
    )

    rijkl[1, 0, 0, 1] = (
        M * r * X * (2 * Delta + a**2 * smp.sin(theta) ** 2) / Sigma**4
    )

    rijkl[2, 0, 0, 2] = -(
        M * r * X * (Delta + 2 * a**2 * smp.sin(theta) ** 2) / Sigma**4
    )

    rijkl[1, 3, 1, 3] = (
        -M
        * r
        * X
        * ((r**2 + a**2) ** 2 + 2 * Delta * a**2 * smp.sin(theta) ** 2)
        * smp.sin(theta) ** 2
        / Sigma**4
    )
    rijkl[2, 3, 2, 3] = (
        M
        * r
        * X
        * (2 * (r**2 + a**2) ** 2 + Delta * a**2 * smp.sin(theta) ** 2)
        * smp.sin(theta) ** 2
        / Sigma**4
    )

    if config == "ulll":
        return antisymmetrize(arr=rijkl)
    elif config == "llll":
        rulll = antisymmetrize(arr=rijkl)
        return np.einsum(
            "ij,jklm->iklm", metric_tensor(r=r, theta=theta, a=a, M=M), rulll
        )

    else:
        # #!Config: lluu
        gkinv = np.linalg.inv(metric_tensor(r=r, theta=theta, a=a, M=M))

        rllll = antisymmetrize(
            arr=np.einsum(
                "ij,jklm->iklm", metric_tensor(r=r, theta=theta, a=a, M=M), rijkl
            )
        )
        return np.einsum("kjlm,al,bm->kjab", rllll, gkinv, gkinv)

# !jupytext --output Einstenpy.py Einsteinpy.ipynb
