import numpy as np
from scipy.integrate import solve_ivp

# Universal inner problem: Y'''' = 2 Y Y''.
# Shooting normalization of the manuscript: Y(0)=-1, Y'(0)=0, K(0)=c, K'(0)=0.

def rhs(t, s):
    Y, P, K, R = s
    return [P, K, R, 2.0*Y*K]

def evK(t, s):  return s[2]
evK.terminal = True; evK.direction = -1.0
def evR(t, s):  return s[3]
evR.terminal = True; evR.direction = 1.0

def classify(c, Xmax=40.0):
    sol = solve_ivp(rhs, [0, Xmax], [-1.0, 0.0, c, 0.0],
                    events=[evK, evR], method='DOP853',
                    rtol=1e-13, atol=1e-14, max_step=0.05)
    tK = sol.t_events[0]; tR = sol.t_events[1]
    if tK.size and (not tR.size or tK[0] < tR[0]):
        return 1  # type I: K hits 0 with K'<0
    if tR.size:
        return 2  # type II: K'=0 with K>0
    return 0

# sanity: manuscript lemmas
print("c=0.3 ->", classify(0.3), " (expect 1)")
print("c=45  ->", classify(45.0), " (expect 2)")

lo, hi = 0.3, 45.0
assert classify(lo) == 1 and classify(hi) == 2
for it in range(70):
    mid = 0.5*(lo+hi)
    t = classify(mid)
    if t == 1: lo = mid
    elif t == 2: hi = mid
    else:
        # undecided: treat as separatrix-side; refine bracket around it
        # decide by longer integration
        t2 = classify(mid, Xmax=80.0)
        if t2 == 1: lo = mid
        elif t2 == 2: hi = mid
        else:
            print("undecided at", mid); break
c_star = 0.5*(lo+hi)
print("c* =", repr(c_star), " bracket width =", hi-lo)

# Integrate only through the numerically reliable early tail. Because the
# separatrix is exponentially ill-conditioned under forward integration, later
# values are not used. The cutoff below is reached near X=5.
sol = solve_ivp(rhs, [0, 5.35], [-1.0, 0.0, c_star, 0.0],
                method='DOP853', rtol=1e-13, atol=1e-15,
                dense_output=True, max_step=0.02)
Y, P, K, R = sol.y
t = sol.t
# Choose an early-tail cutoff where K~1e-8, before forward-integration
# contamination of the exponentially delicate separatrix.
mask = (K > 0) & (R < 0) & (K < 1e-8)
if not mask.any():
    raise RuntimeError('reliable early-tail cutoff K<1e-8 was not reached')
icut = np.where(mask)[0][0]
Xc, Yc, Pc, Kc, Rc = t[icut], Y[icut], P[icut], K[icut], R[icut]
q = 2.0*Yc
# tail integrals for K ~ e^{-Theta}, Theta' = sqrt(q):
# int_X^inf K ds ~ K/sqrt(q) * (1 + O(q^{-3/2})); int_X^inf (P_inf - P) ~ K/q
Pinf = Pc + Kc/np.sqrt(q)
Bstar = Yc - Pinf*Xc - Kc/q
print("cutoff X =", Xc, " K =", Kc)
print("P_inf =", repr(Pinf))
print("B*    =", repr(Bstar))

kappa = Pinf**(-1.0/3.0)
Y0_0   = -kappa**2
K0_0   = kappa**4 * c_star
B0     = kappa**2 * Bstar
print()
print("Y0(0)   =", repr(Y0_0),  "   ms: -0.3933180528")
print("Y0''(0) =", repr(K0_0),  "   ms:  0.6898804165")
print("B0      =", repr(B0),    "   ms: -1.3307771994")

# zero of Y0: zero of Y at X0_raw, then X0 = X0_raw/kappa
from scipy.optimize import brentq
i0 = np.where(Y > 0)[0][0]
X0_raw = brentq(lambda x: sol.sol(x)[0], t[i0-1], t[i0], xtol=1e-14)
X0 = X0_raw/kappa
print("X0      =", repr(X0),    "   ms:  1.1090768620")
print()
print("A0 = -Y0(0) =", repr(-Y0_0), "  ms A0: 0.3933180528")
