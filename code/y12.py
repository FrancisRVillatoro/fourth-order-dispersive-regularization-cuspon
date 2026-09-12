import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def rhs0(t, s):
    Y, P, K, R = s
    return [P, K, R, 2.0*Y*K]

c_star = 4.459498847498201
solU = solve_ivp(rhs0, [0, 5.35], [-1.0, 0.0, c_star, 0.0],
                 method='DOP853', rtol=1e-13, atol=1e-16,
                 dense_output=True, max_step=0.01)

Pinf = 4.054004070          # from cutoff study, +-~1e-8
kap  = Pinf**(-1.0/3.0)

def Y0f(X):   # scaled layer profile and second derivative
    s = solU.sol(kap*X)
    return kap**2*s[0], kap**4*s[2]

B0 = kap**2*(-3.3834633950)   # B* from early-cutoff value

# ---- pass 1: Y1 ----
def rhs1(t, s):
    Y0, K0 = Y0f(t)
    out = np.empty(12)
    for j, f in zip(range(3), (1.0, 0.0, 0.0)):
        y, p, k, r = s[4*j:4*j+4]
        out[4*j:4*j+4] = (p, k, r, 2*Y0*k + 2*K0*y + f)
    return out

def solve_Y1(Xmax):
    ic = np.zeros(12); ic[4] = 1.0; ic[10] = 1.0   # phi_p, phi_A(0)=1, phi_B''(0)=1
    sol = solve_ivp(rhs1, [0, Xmax], ic, method='DOP853',
                    rtol=1e-12, atol=1e-14, dense_output=True, max_step=0.02)
    z = Xmax + B0
    tgtK = -1/(2*z) - 1/(2*z**4)
    tgtP = -0.5*np.log(z) + 1/(6*z**3)
    e = sol.sol(Xmax)
    Amat = np.array([[e[6],  e[10]],
                     [e[5],  e[9]]])
    rvec = np.array([tgtK - e[2], tgtP - e[1]])
    A, B = np.linalg.solve(Amat, rvec)
    def Y1(X, d=0):
        e = sol.sol(X)
        return e[0+d] + A*e[4+d] + B*e[8+d]
    C1 = Y1(Xmax) + 0.5*z*np.log(z) - 0.5*z + 1/(12*z**2)
    return A, B, C1, sol, (A, B)

for Xmax in (6.0, 7.0, 8.0, 8.4):
    A, B, C1, *_ = solve_Y1(Xmax)
    print(f"Xmax={Xmax}:  Y1(0)={A:.9f}  Y1''(0)={B:.9f}  C1={C1:.9f}")

Xmax = 8.0
A1c, B1c, C1, sol1, coef = solve_Y1(Xmax)

# ---- pass 2: Y2 with forcing 2*Y1*Y1'' ----
def rhs2(t, s):
    Y0, K0 = Y0f(t)
    e = sol1.sol(t)
    Y1v = e[0] + coef[0]*e[4] + coef[1]*e[8]
    K1v = e[2] + coef[0]*e[6] + coef[1]*e[10]
    F = 2.0*Y1v*K1v
    out = np.empty(12)
    for j, f in zip(range(3), (F, 0.0, 0.0)):
        y, p, k, r = s[4*j:4*j+4]
        out[4*j:4*j+4] = (p, k, r, 2*Y0*k + 2*K0*y + f)
    return out

def solve_Y2(Xm):
    ic = np.zeros(12); ic[4] = 1.0; ic[10] = 1.0
    sol = solve_ivp(rhs2, [0, Xm], ic, method='DOP853',
                    rtol=1e-12, atol=1e-14, dense_output=True, max_step=0.02)
    z = Xm + B0; lz = np.log(z)
    tgtK = -(lz - 1)/(4*z) + C1/(2*z**2) + (-0.5*lz + 5/6)/z**4
    tgtP = -lz**2/8 + lz/4 - 0.25 - C1/(2*z) + lz/(6*z**3) - 2/(9*z**3)
    e = sol.sol(Xm)
    Amat = np.array([[e[6],  e[10]],
                     [e[5],  e[9]]])
    rvec = np.array([tgtK - e[2], tgtP - e[1]])
    A, B = np.linalg.solve(Amat, rvec)
    return A, B, sol, (A, B)

for Xm in (6.0, 7.0, 8.0):
    A2, B2, *_ = solve_Y2(Xm)
    print(f"Xm={Xm}:  Y2(0)={A2:.8f}  Y2''(0)={B2:.8f}")

A2c, B2c, sol2, coef2 = solve_Y2(8.0)

print()
print(f"Y1(0)   = {A1c:.9f}    ms -1.1065076820  ->  A1 = {-A1c:.9f}")
print(f"Y1''(0) = {B1c:.9f}    ms  0.3022408770")
print(f"Y2(0)   = {A2c:.8f}    ms -0.0131925717  ->  A2 = {-A2c:.8f}")
print(f"Y2''(0) = {B2c:.8f}    ms  0.6573620821")

# ---- interface coefficients ----
X0 = brentq(lambda x: Y0f(x)[0], 0.5, 2.0, xtol=1e-14)
h = 1e-6
Y0p  = (Y0f(X0+h)[0]-Y0f(X0-h)[0])/(2*h)
Y0pp = Y0f(X0)[1]
e = sol1.sol(X0)
Y1v  = e[0] + coef[0]*e[4] + coef[1]*e[8]
Y1p  = e[1] + coef[0]*e[5] + coef[1]*e[9]
e2 = sol2.sol(X0)
Y2v = e2[0] + coef2[0]*e2[4] + coef2[1]*e2[8]

X1 = -Y1v/Y0p
X2 = -(Y2v + X1*Y1p + 0.5*X1**2*Y0pp)/Y0p
print()
print(f"X0 = {X0:.10f}   ms  1.1090768620")
print(f"X1 = {X1:.9f}   ms  1.4615818195")
print(f"X2 = {X2:.8f}   ms -1.6180901746")
