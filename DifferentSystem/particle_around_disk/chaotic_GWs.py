# * Imports

import numpy as np 
import matplotlib.pyplot as plt #* For Plotting
import scipy as sp
from scipy import integrate #* For Solving Differential Equations
import snoop  #* For Debugging Purpose

from numbalsoda import lsoda_sig, lsoda #Faster Method
import numba as nb #For Speeding Up
import scipy.integrate as integrate
from scipy.interpolate import InterpolatedUnivariateSpline as spline


#! Parameters

mu = 1.0  #! This Value was not given in paper. I solved it
#! using the constraint on the total energy
M = 1.0
L = 1.0
z0 = 0.5
alpha = 0.01  # TODO Important Parameter.
ics = np.array(
    [1.2, 0.0, 0.0, 0.0, 0.76, L]
)  # * Initial conditions w,z,phi, pw,pz,pphi from paper
t = np.arange(0.0, 100000.0, 1)  # Time Steps for evolution

tau = 100  # Renormalization time
tf = 1000000.0  # Final time
steps = 10000  # *Number of time steps between each time interval
# Theory Parameters

Theta = np.pi / 2
Psi = 0.0
R = 1e10

tm = np.linspace(0.0, 500000.0, 5000000)  # Time for integrations
tm2 = np.linspace(0.0, 10000.0, 100000)

## Noise level in the Science Requirement Document
sqSnoise_SciRD = [3e-15, 15e-12]

## Speed of light
c = 299792458.0

## Year in seconds
year = 31558149.763545603

## Astronomical Unit in seconds
AUs = 499.00478383615643

## Armlength in seconds
L_m = 2.5e9
L = L_m / c

lower = -1e-21
upper = 1e-21

#
@nb.njit
def logcosh(x):
    # s always has real part >= 0
    s = np.sign(x) * x
    p = np.exp(-2 * s)
    return s + np.log1p(p) - np.log(2)
# @nb.njit


# @snoop
@nb.njit
def H(mu: np.float64,ics: np.ndarray,alpha: np.int32) -> np.ndarray:
    """

    Args:
        mu (float): Mass of the Test Particle
        ics (np.ndarray): The Phase Space Coordinates of the particle
                        Shape(n,4)

    Returns:
        Energy of the particel (np.ndarray): Shape (n,)
    """
    if ics.ndim==1:
        w,z,phi,pw,pz,pphi=ics[0],ics[1],ics[2],ics[3],ics[4],ics[5]
    else:
        w,z,phi,pw,pz,pphi=ics[:,0],ics[:,1],ics[:,2],ics[:,3],ics[:,4],ics[:,5]
    L= pphi

    return mu*(pw**2/(2*mu) + pz**2/(2*mu) + L**2/(2* mu**2*w**2) - M/np.sqrt(w**2 + z**2) + alpha*z0*logcosh(z/z0))


@nb.njit
def evolve(t: float, x: np.ndarray, alpha: int) -> np.ndarray:
    """
    Function for solving Differential Equation
    Args:
        t (float): Time Step
        x (np.ndarray): Shape (n,)

    Returns:
        Output (np.ndarray): Shape (n,)
    """
    w, z, phi, pw, pz, pphi = x
    dw = pw
    dz = pz
    dphi = L / (mu * w**2)

    dpw = -mu * (-(L**2) / (mu**2 * w**3) + M * w / ((w**2 + z**2) ** (3 / 2)))
    dpz = -mu * (M * z / ((w**2 + z**2) ** (3 / 2)) + alpha * np.tanh(z / z0))
    dpphi = 0.0

    return np.array((dw, dz, dphi, dpw, dpz, dpphi))


# @nb.cfunc('')
@nb.njit
def poincare(t: float, y: np.ndarray, args: np.int32) -> float:
    """Creating Poincare Section

    Args:
        t (float): time Step
        y (np.ndarray): Shape (n,4)
        args (np.int32) : Not needed. But used for compatibility for
                            below methods

    Returns:
        float: The value where the function gets zero (z value in this case)
    """
    return y[1]


def chaotic_lyp(t, Y, alpha=alpha):
    """
    Function for solving the system of equation simultaneously with
    the Jacobian

    """

    # First system of equations

    w, z, phi, pw, pz, pphi = Y[:6]
    dw = pw
    dz = pz
    dphi = L / (mu * w**2)

    dpw = -mu * (-(L**2) / ((mu**2) * w**3) + M * w / ((w**2 + z**2) ** (3 / 2)))
    dpz = -mu * (M * z / ((w**2 + z**2) ** (3 / 2)) + alpha * np.tanh(z / z0))
    dpphi = 0.0

    dhdwdw = -mu * (
        3 * L**2 / ((mu**2) * (w**4))
        - 3 * M * w**2 / ((w**2 + z**2) ** (5 / 2))
        + M / ((w**2 + z**2) ** (3 / 2))
    )
    dhdwdz = 3 * M * mu * w * z / ((w**2 + z**2) ** (5 / 2))
    dhdzdz = -mu * (
        -3 * M * z**2 / ((w**2 + z**2) ** (5 / 2))
        + M / ((w**2 + z**2) ** (3 / 2))
        + alpha * (1 / (np.cosh(z / z0)) ** 2) / z0
    )

    # Now calculate the jacobian
    J = np.array(
        [
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            [-2 * L / (mu * w**3), 0.0, 0.0, 0.0, 0.0, 1 / (mu * w**2)],
            [dhdwdw, dhdwdz, 0.0, 0.0, 0.0, 2 * L / (mu * w**3)],
            [dhdwdz, dhdzdz, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ]
    )

    dY = Y[6:].reshape((6, 1))
    dY_dt = np.dot(J, dY)

    return np.concatenate(([dw, dz, dphi, dpw, dpz, dpphi], dY_dt.flatten()))


#! Faster implementaton of above function
@nb.cfunc(lsoda_sig)
def faster_chaotic_system(t, u, du, p):
    """
    The same system as above but faster than above method. This uses Numba
    and Numbalsoda for faster evaluations
    """
    u_ = nb.carray(u, (12,))
    p_ = nb.carray(p, (1,))
    alpha = p_[0]
    w, z, phi, pw, pz, pphi = u_[:6]
    dw = pw
    dz = pz
    dphi = L / (mu * w**2)

    dpw = -mu * (-(L**2) / ((mu**2) * w**3) + M * w / ((w**2 + z**2) ** (3 / 2)))
    dpz = -mu * (M * z / ((w**2 + z**2) ** (3 / 2)) + alpha * np.tanh(z / z0))
    dpphi = 0.0

    dhdwdw = -mu * (
        3 * L**2 / ((mu**2) * w**4)
        - 3 * M * w**2 / ((w**2 + z**2) ** (5 / 2))
        + M / ((w**2 + z**2) ** (3 / 2))
    )
    dhdwdz = 3 * M * mu * w * z / ((w**2 + z**2) ** (5 / 2))
    dhdzdz = -mu * (
        -3 * M * z**2 / ((w**2 + z**2) ** (5 / 2))
        + M / ((w**2 + z**2) ** (3 / 2))
        + alpha * (1 / (np.cosh(z / z0)) ** 2) / z0
    )

    # Now calculate the jacobian

    J = np.array(
        [
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            [-2 * L / (mu * w**3), 0.0, 0.0, 0.0, 0.0, 1 / (mu * w**2)],
            [dhdwdw, dhdwdz, 0.0, 0.0, 0.0, 2 * L / (mu * w**3)],
            [dhdwdz, dhdzdz, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ]
    )

    dY = u_[6:].reshape((6, 1))
    dY_dt = np.dot(J, dY)
    dydtf = dY_dt.flatten()
    du_ = np.zeros((12,), dtype=nb.float64)
    du_[:6] = np.array([dw, dz, dphi, dpw, dpz, dpphi])
    du_[6:] = dydtf

    for i in range(len(du_)):
        du[i] = du_[i]


@nb.njit
def evolve_lyap(ics0, tau, tf, alpha, steps, sanity_check=0):

    ti = 0.0  # * Initial Time
    timestep = np.arange(ti, tf, tau)  # * Time interval
    size = timestep.shape[0]

    X1t = np.zeros((size - 1))  # *Lyapunaov Exponents
    sum = 0.0
    ics = ics0.copy()  # *We do not want to change our initial conditions
    h = np.zeros(size - 1)  # *Check for energy
    sol = np.zeros((size, 6))
    sol[0, :] = ics[:6]
    # traj=np.zeros((steps,6),dtype=np.float64)

    for i in range(0, size - 1):
        timeint = np.linspace(
            timestep[i], timestep[i + 1], steps
        )  #!Our time interval with 10000 points between them

        usol_, _ = lsoda(
            pfunctptr,
            ics,
            t_eval=timeint,
            data=np.array((alpha)),
            rtol=1e-12,
            atol=1e-12,
            mxstep=100000,
        )  # * Solver
        # if i == 0:
        #     traj = np.zeros_like(usol_)
        #     traj = usol_

        usol = usol_[-1]

        # Uncomment following line to use scipy solver and comment above two lines
        # usol=integrate.solve_ivp(chaotic_lyp_paper,[0.,tf],ics0,'LSODA',t_eval=timeint,atol=1e-10,rtol=1e-10).y[:,-1]

        xk = usol[:6]  # Orbit at t= tau
        wk = usol[6:]  # Deviations at t=tau

        # *---------------------The following algorithm is from the paper

        alphak = np.linalg.norm(wk)
        sum = sum + np.log(alphak)
        X1t[i] = sum / timestep[i + 1]
        wk0 = wk / alphak

        # *----------------------Lyapunaov Exponents calculated above
        #! Now for another time step. Set the current solution as the initial conditions for next time steps
        ics[:6] = xk
        ics[6:] = wk0

        #! ------------------------Sanity Check--------------------------
        # if sanity_check!=0:
        h[i] = H(mu, ics, alpha)  #! Current Energy
        sol[i + 1, :] = ics[:6]
        #! Comment above line if you do not want to perform the sanity check

    # if sanity_check!=0:
    #     return X1t,h,sol, traj

    return X1t, sol, h


@nb.cfunc(lsoda_sig)
def faster_evolve(t, u, du, p):
    """
    Function for solving Differential Equation
    Args:
        t (float): Time Step
        x (np.ndarray): Shape (n,)

    Returns:
        Output (np.ndarray): Shape (n,)
    """
    # w,z, phi,pw,pz,pphi=x
    p_ = nb.carray(p, (1,))
    alpha = p_[0]
    du[0] = u[3]
    du[1] = u[4]
    du[2] = L / (mu * u[0] ** 2)

    du[3] = -mu * (
        -(L**2) / (mu**2 * u[0] ** 3) + M * u[0] / ((u[0] ** 2 + u[1] ** 2) ** (3 / 2))
    )
    du[4] = -mu * (
        M * u[1] / ((u[0] ** 2 + u[1] ** 2) ** (3 / 2)) + alpha * np.tanh(u[1] / z0)
    )
    du[5] = 0.0


@nb.njit(
    nb.float64[:](
        nb.float64,
        nb.float64[:, :],
        nb.float64,
        nb.float64,
        nb.float64,
        nb.float64,
        nb.float64,
    ),
    parallel=True,
)
def hplus(t, dH, R, mu, Theta, Psi, alpha):
    """

    Args:
        t ([float]): Time, not needed in this function
        dH ([ndarray]): Solution of differential Equation at various time steps, shape (n,6)
        R ([float]): Distance of observer and the system
        mu ([float]): Mass of the particle
        Theta ([float]): Location of the observer
        Psi ([float]): Location of the observer
        alpha ([float]): Parameter of the theory

    Returns:
        Expression of hplus ([ndarray]) at various time: Shape (n,6)
    """
    # Compute the hplus expression using NumPy

    r, zc, phi, pr, pzc, pphi = (
        dH[:, 0],
        dH[:, 1],
        dH[:, 2],
        dH[:, 3],
        dH[:, 4],
        dH[:, 5],
    )
    dr_dt = pr / mu
    dzc_dt = pzc / mu
    dPhi_dt = L / (mu**2 * r**2)
    d2r_dt2 = -(-(L**2) / (mu**2 * r**3) + M * r / ((r**2 + zc**2) ** (3 / 2)))
    d2zc_dt2 = -(M * zc / ((r**2 + zc**2) ** (3 / 2)) + alpha * np.tanh(zc / z0))
    d2Phi_dt2 = 0.0

    expression = (
        (1 / (2 * R))
        * mu
        * (
            -2
            * np.sin(Theta) ** 2
            * (dr_dt**2 - 2 * dzc_dt**2 + r * d2r_dt2 - 2 * zc * d2zc_dt2)
            + (3 + np.cos(2 * Theta))
            * np.sin(2 * phi)
            * (
                np.sin(2 * Psi) * (dr_dt**2 + r * (-2 * r * dPhi_dt**2 + d2r_dt2))
                - np.cos(2 * Psi) * r * (4 * dr_dt * dPhi_dt + r * d2Phi_dt2)
            )
            + (3 + np.cos(2 * Theta))
            * np.cos(2 * phi)
            * (
                np.cos(2 * Psi) * (dr_dt**2 + r * (-2 * r * dPhi_dt**2 + d2r_dt2))
                + r * np.sin(2 * Psi) * (4 * dr_dt * dPhi_dt + r * d2Phi_dt2)
            )
            + 2
            * np.sin(2 * Theta)
            * np.sin(phi)
            * (
                -np.sin(Psi)
                * (
                    2 * dr_dt * dzc_dt
                    + zc * d2r_dt2
                    + r * (-zc * dPhi_dt**2 + d2zc_dt2)
                )
                + np.cos(Psi)
                * (2 * (zc * dr_dt + r * dzc_dt) * dPhi_dt + r * zc * d2Phi_dt2)
            )
            - 2
            * np.cos(phi)
            * np.sin(2 * Theta)
            * (
                np.cos(Psi)
                * (
                    2 * dr_dt * dzc_dt
                    + zc * d2r_dt2
                    + r * (-zc * dPhi_dt**2 + d2zc_dt2)
                )
                + np.sin(Psi)
                * (2 * (zc * dr_dt + r * dzc_dt) * dPhi_dt + r * zc * d2Phi_dt2)
            )
        )
    )


    return expression
# @snoop
@nb.njit(
    nb.float64[:](
        nb.float64,
        nb.float64[:, :],
        nb.float64,
        nb.float64,
        nb.float64,
        nb.float64,
        nb.float64,
    ),
    parallel=True,
)
# Define the Python function 'hcross'
def hcross(t, dH, R, mu, Theta, Psi, alpha):
    """

    Args:
        t ([float]): Time, not needed in this function
        dH ([ndarray]): Solution of differential Equation at various time steps, shape (n,6)
        R ([float]): Distance of observer and the system
        mu ([float]): Mass of the particle
        Theta ([float]): Location of the observer
        Psi ([float]): Location of the observer
        alpha ([float]): Parameter of the theory

    Returns:
        Expression of hcross ([ndarray]) at various time: Shape (n,6)
    """
    # Compute the hcross expression using NumPy

    r, zc, phi, pr, pzc, pphi = (
        dH[:, 0],
        dH[:, 1],
        dH[:, 2],
        dH[:, 3],
        dH[:, 4],
        dH[:, 5],
    )
    dr_dt = pr
    dzc_dt = pzc
    dPhi_dt = L / (mu * r**2)
    d2r_dt2 = -mu * (-(L**2) / (mu**2 * r**3) + M * r / ((r**2 + zc**2) ** (3 / 2)))
    d2zc_dt2 = -mu * (M * zc / ((r**2 + zc**2) ** (3 / 2)) + alpha * np.tanh(zc / z0))
    d2Phi_dt2 = 0.0

    # Compute the expression using NumPy
    expression = (
        (1 / R)
        * 2
        * mu
        * (
            -np.cos(Theta) * np.sin(2 * (Psi - phi)) * dr_dt**2
            + 2
            * dr_dt
            * (
                np.sin(Theta) * np.sin(Psi - phi) * dzc_dt
                + (
                    2 * np.cos(Theta) * np.cos(2 * (Psi - phi)) * r
                    - np.cos(Psi - phi) * np.sin(Theta) * zc
                )
                * dPhi_dt
            )
            + np.cos(Theta)
            * r
            * (
                np.sin(2 * (Psi - phi)) * (2 * r * dPhi_dt**2 - d2r_dt2)
                + np.cos(2 * (Psi - phi)) * r * d2Phi_dt2
            )
            + np.sin(Theta)
            * (
                np.sin(Psi - phi) * zc * d2r_dt2
                - r
                * (
                    2 * np.cos(Psi - phi) * dzc_dt * dPhi_dt
                    + np.sin(Psi - phi) * (zc * dPhi_dt**2 - d2zc_dt2)
                    + np.cos(Psi - phi) * zc * d2Phi_dt2
                )
            )
        )
    )

    return expression


def de_dw(t, hp, hc):
    """
    Calculate dE/dw

    Args:
        t (ndarray): Time array between which to calculate the FFT, Shape: (n,)
        hp (ndarray): hplus array, Shape: (n,)
        hc (ndarray): hcross array, Shape: (n,)

    Returns:
        dEdw (ndarray): dEdw, shape(n//2,)
        df (ndarray): Fourier Transform Frequencies, shape(n//2,)
    """
    # window=np.hanning(t.size)
    fhp = sp.fft.fft(hp, norm="ortho")
    fhc = sp.fft.fft(hc, norm="ortho")
    N = len(fhc)

    dt = t[1] - t[0]
    # df=sp.fft.fftfreq(len(fhp),dt)[1:N//2]
    df = (np.arange(len(t)) / (dt * N))[: N // 2]
    # const=(3e8)**3/(4*6.674e-11)

    E = (np.abs(fhp) ** 2 + np.abs(fhc) ** 2)[: N // 2]
    dEdw = (E * df**2) * R**2

    return dEdw, df


def PSD_Noise_components(fr, sqSnoise):
    [sqSacc_level, sqSoms_level] = sqSnoise
    # sqSacc_level: Amplitude level of acceleration noise [3e-15]
    # sqSoms_level: Amplitude level of OMS noise [15e-12]

    ### Acceleration noise
    Sa_a = sqSacc_level**2 * (1.0 + (0.4e-3 / fr) ** 2) * (1.0 + (fr / 8e-3) ** 4)
    Sa_d = Sa_a * (2.0 * np.pi * fr) ** (-4.0)
    Sa_nu = Sa_d * (2.0 * np.pi * fr / c) ** 2

    ### Optical Metrology System
    Soms_d = sqSoms_level**2 * (1.0 + (2.0e-3 / fr) ** 4)
    Soms_nu = Soms_d * (2.0 * np.pi * fr / c) ** 2

    return [Sa_nu, Soms_nu]


def PSD_Noise_X20(fr, sqSnoise):
    #! This gives S_n for X_{2.0} as described in Eq 120 in LISA SciRD
    [Sa_nu, Soms_nu] = PSD_Noise_components(fr, sqSnoise)
    phiL = 2 * np.pi * fr * L
    return (
        64
        * (np.sin(phiL)) ** 2
        * (np.sin(2 * phiL)) ** 2
        * (Soms_nu + Sa_nu * (3 + np.cos(2 * phiL)))
    )


def GenGal(fr, Tobs, A, f2, alp, a1, b1, ak, bk):

    f1 = 10 ** (a1 * np.log10(Tobs) + b1)
    fk = 10 ** (ak * np.log10(Tobs) + bk)

    Sgal = (
        0.5
        * A
        * fr ** (-7.0 / 3.0)
        * np.exp(-((fr / f1) ** alp))
        * (1.0 + np.tanh((fk - fr) / f2))
    )

    return Sgal


def colour(f, log_fr, log_psds):
    with np.errstate(divide="ignore"):
        return np.exp(np.interp(np.log(f), log_fr, log_psds))


def normalize(
    x, newRange=(0, 1)
):  # x is an array. Default range is between zero and one
    xmin, xmax = np.min(x), np.max(x)  # get max and min from input array
    norm = (x - xmin) / (xmax - xmin)  # scale between zero and one

    if newRange == (0, 1):
        return norm  # wanted range is the same as norm
    elif newRange != (0, 1):
        return (
            norm * (newRange[1] - newRange[0]) + newRange[0]
        )  # scale to a different range.
    # add other conditions here. For example, an error message


#####################Code Execution starts###############

print("Executing the main code now....\n")

sol = integrate.solve_ivp(
    evolve,
    (0.0, t[-1]),
    ics,
    "LSODA",
    t_eval=t,
    events=poincare,
    dense_output=True,
    args=(alpha,),
    atol=1e-11,
    rtol=1e-11,
)

mask = sol.y_events[0][:, 4] > 0  # * Get those values where v_z is greater than zero


# * Poincare Plot

plt.scatter(sol.y_events[0][mask, 0], sol.y_events[0][mask, 3], s=0.2)
plt.title(r"Poincare Plot for z=0")
plt.xlabel(r"w")
plt.ylabel(r"$p_w$")
plt.show()


energy = H(mu, sol.y.T, alpha=alpha)  # Get the energy of given solutions


plt.plot(energy-energy[0])
plt.title("Energy vs time")
plt.xlabel("Time")
plt.ylabel("Energy")
plt.show()

# * The value of energy is within -0.2 (+/- 39e^(-8))

print("Calculating Lyapunov Exponents\n")

pfunctptr = faster_chaotic_system.address


u0 = np.copy(ics)  # Initial condition for x(0)
w0 = np.random.normal(0.0, 1.0, (6))
w0 = w0 / np.linalg.norm(w0)  # Initial unitary deviation vector w(0)

ics0 = np.array((u0, w0)).flatten()  # Combined Initial vector

_, sol_, traj_ = evolve_lyap(
    ics0, 10.0, 100.0, alpha, steps, 0
)  #!Run it for the first time so that numba understands the call signature


lyaps1 = evolve_lyap(ics0, tau, tf, alpha, steps)  # * Lyapunaov Exponents.

#! If you get a warning the t+h_=t on the next step, ignore it.


plt.plot(np.arange(0.0, tf, tau)[:-1], lyaps1[2])
plt.title("Energy difference when calculating Lyapunov Exponent")
plt.show()

plt.figure()
plt.plot(
    np.arange(0.0, tf, tau)[:-2],
    lyaps1[0][:-1],
    label="Maximum Lyapunov Exponents",
    linewidth=2,
    linestyle="--",
)
ax = plt.gca()
ax.set_xscale("log", base=10)
ax.set_yscale("log", base=10)
plt.xlabel("Time")
plt.ylabel("Lyapunov Exponent")
plt.legend()
# plt.xlim(1e4,1e6)
plt.grid(True)
plt.title(r"Time Evolution of Lyapunov Exponents for $\alpha$={0}".format(alpha))
plt.show()


# * ----------------------------------------------------

print("Calculating GWs...\n")

funcptr = faster_evolve.address


# Solve the system
usol10, success = lsoda(
    funcptr, ics, tm2, data=np.array((alpha)), atol=1e-12, rtol=1e-12
)

hp10 = hplus(0, usol10, R, mu, Theta, Psi, alpha)  # hplus
hc10 = hcross(0, usol10, R, mu, Theta, Psi, alpha)  # hcross


fig = plt.figure()
# fig.set_size_inches(10.5, 10.5)
plt.plot(tm2, R * hp10, "r--")
plt.title(r"r$h_+$, $\alpha$={0}".format(alpha))
plt.ylabel(r"r$h_+$")
plt.xlabel("t")
plt.show()

fig = plt.figure()
# fig.set_size_inches(10.5, 10.5)
plt.plot(tm2, R * hc10, "r--")
plt.title(r"$h_x, \alpha$={0}".format(alpha))
plt.ylabel(r"r$h_+$")
plt.xlabel("t")
plt.show()

#!----------------------------------------------------------------

E10, f10 = de_dw(tm2, hp10, hc10)
# * For alpha

plt.plot(np.log10(f10[1:]), np.log10(E10[1:]))
plt.title(r"dE/d$\omega$ vs f for $\alpha=$ {0}".format(alpha))
plt.ylabel(r"$Log_{10}(\frac{dE}{dw}$)")
plt.xlabel(r"$Log_{10}f$")
plt.show()


#! Noise + Signal
print("Adding noise to Signal....\n")

noise = np.random.normal(0.0, 1.0, len(hp10))  # * Normally distributed Random Noise

nhp10 = hp10 + noise  # *Adding noise to signals
nhc10 = hc10 + noise


fig = plt.figure()
# fig.set_size_inches(10.5, 10.5)
plt.plot(tm2, R * nhc10, "r--")
plt.title(r"$h_x+ noise, \alpha=${0}".format(alpha))
plt.ylabel(r"r$h_+$")
plt.xlabel("t")
plt.show()

fig = plt.figure()
# fig.set_size_inches(10.5, 10.5)
plt.plot(tm2, R * nhp10, "r--")
plt.title(r"$h_++$noise, $\alpha$={0}".format(alpha))
plt.ylabel(r"r$h_+$")
plt.xlabel("t")
plt.show()


#! Get the energy spectra for noise with signals
NE10, nf10 = de_dw(tm2, nhp10, nhc10)

fig, ax = plt.subplots(nrows=1, ncols=1)
ax.plot(np.log10(f10[1:]), np.log10(E10[1:]))
ax.set_title("Signal case")
plt.ylabel(r"$Log_{10}(\frac{dE}{dw}$)")
plt.xlabel(r"$Log_{10}f$")
plt.tight_layout()
plt.savefig("psd_signal.png", dpi=300)

fig, ax = plt.subplots(nrows=1, ncols=1)
# fig.set_size_inches(18.5, 10.5)
ax.plot(np.log10(nf10[1:]), np.log10(NE10[1:]))

# fig.suptitle(r"dE/d$\omega$ vs f for $\alpha=10.$")
ax.set_title("Signal + Noise case")

plt.ylabel(r"$Log_{10}(\frac{dE}{dw}$)")
plt.xlabel(r"$Log_{10}f$")
plt.tight_layout()
plt.savefig("psd_signal_noise.png", dpi=300)

freqs, AvFXp2_Raw = np.load("./AvFXp2_Raw.npy")

## Minimal frequency for plotting
fMin = 1e-5

## Minimal frequency for plotting
fMax = 1

fr_ = 10 ** np.linspace(np.log10(fMin), np.log10(fMax), 100000)
AvFXp2_ = np.interp(fr_, freqs, AvFXp2_Raw)
[Sa_nu_, Soms_nu_] = PSD_Noise_components(fr_, sqSnoise_SciRD)
phiL = 2 * np.pi * fr_ * L
S_hX_SA_ = (Soms_nu_ + Sa_nu_ * (3 + np.cos(2 * phiL))) / (phiL**2 * AvFXp2_ / 4**2)
# S_h_SA = S_hX_SA / 2

log_fr_ = np.log(fr_)
log_psds_ = np.log(S_hX_SA_ + S_gal_)

Nhp10 = normalize(hp10, (lower, upper))
Nhc10 = normalize(hc10, (lower, upper))

n_, dt_ = 100000, 1.0 / 16.0
f_ = np.fft.rfftfreq(n_, dt_)
rng = np.random.default_rng()
x_f = rng.normal(0, 0.5, len(f_)) + 1j * rng.normal(0, 0.5, len(f_))
x_f *= np.sqrt(n_ / dt_)

x_f[0] = np.abs(x_f[0])
if len(f) % 2 == 0:
    x_f[-1] = np.abs(x_f[-1])
x_f *= np.sqrt(colour(f_, log_fr_, log_psds_))
noise = np.fft.irfft(x_f)

lnoise = noise.copy()

lnoise[10000:20000] = noise[10000:20000] + Nhp10[10000:20000]

template_fft = abs(np.fft.fft(Nhp10)) ** 2 + abs(np.fft.fft(Nhc10)) ** 2

dt = t[1] - t[0]
fr = np.fft.fftfreq(len(hp10), 1 / dt)

ind = np.argwhere(fr < 0.0)[0, 0]
pfr = fr[1:ind]
template_fft = template_fft[1:ind]


SAv= np.interp(pfr, freqs,AvFXp2_Raw)

SRXSAX20 = (
    (2 * np.pi * pfr * (L_m / c)) ** 2
    * (np.sin(2 * np.pi * pfr * L)) ** 2
    * (2 * np.sin(2 * 2 * np.pi * pfr * L)) ** 2
    * SAv
)

[Sa_nu, Soms_nu] = PSD_Noise_components(pfr, sqSnoise_SciRD)
phiL = 2 * np.pi * pfr * L
S_hX_SA = (Soms_nu + Sa_nu * (3 + np.cos(2 * phiL))) / (phiL**2 * SAv / 4**2)

noiseratio = (template_fft) / (S_hX_SA)

SNR = np.sqrt(np.trapz(noiseratio, pfr) * 4)
print(SNR)

print("Finished execution\n")
