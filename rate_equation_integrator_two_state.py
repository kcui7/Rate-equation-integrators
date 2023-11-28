import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


class RateEquationIntegratorTwoState(object):
    """
    This class is used to numerically solve the rate equations for the following system:
    Two molecules can convert to each other, 
        A <==> B
    This reaction is characterized by the forward rate constant kAB, backward rate constant kBA, and equilibrium constant K.
    Both A and B can react with another molecule M, forming the same product P, these reactions are IRREVERSIBLE
        A + R --> P  k1
        B + R --> P  k2

    Input parameters are:
    kAB (float): rate constant for reaction A --> B
    kBA (float): rate constant for reaction B --> A
    k1 (float): rate constant for reaction A + R --> P
    k2 (float): rate constant for reaction A + R --> P
    cA0 (float): initial concentration of A (in mol/L)
    cB0 (float): initial concentration of B (in mol/L)
    cR0 (float): initial concentration of R (in mol/L)
    """
    def __init__(self, kAB, kBA, k1, k2, cA0, cB0, cR0):

        self.k_consume = np.array([[-k1, 0],
                                   [0, -k2]])

        self.k_eq = np.array([[-kAB, kBA],
                              [kAB, -kBA]])

        self.c0 = np.array([cA0, cB0])
        self.cR0 = cR0

    def evolve(self, t_tot, dt=0.001):
        """
        t_tot (float): total time to evolve the system in reduced time unit. The recommeneded reduced unit is microsecond (us)
        dt (float): time step to evolve the system in reduced time unit. default value = 0.001
        """
        nstep = int(t_tot/dt)

        c_old = self.c0
        cR_old = self.cR0
        ct = []
        cRt = []
        for i in range(nstep):
            ct.append(c_old)
            cRt.append(cR_old)
            k = self.k_consume*cR_old + self.k_eq

            c_new = c_old + np.matmul(k,c_old)*dt
            cR_new = cR_old + np.sum(np.matmul(self.k_consume,c_old))*cR_old*dt
            c_old = c_new
            cR_old = cR_new

        t = np.arange(int(t_tot/dt))*dt
        ct = np.transpose(np.array(ct))
        cRt = np.array(cRt)

        return t, ct, cRt



#=========================================================
# Example 1.1: equilibrium assumption for WEE oxidation
#=========================================================

# parameters for the WEE system
KOH = 7e4
pH = 5
cH = 10**(-pH)   # in M
cOH = 1e-14/cH   # in M

kp = 1e4*cOH   # in us^-1
km = kp/cOH/KOH   # in us^-1
kET1EPT1 = 0.71   # in M^-1 us^-1
kEPT2 = 4e3   # in M^-1 us^-1

PS0 = 2e-5   # in M
cWEE = 5e-3   # in M


c0 = cWEE*np.array([1/(1+KOH*cOH), KOH*cOH/(1+KOH*cOH)])
print(c0)

system = RateEquationIntegratorTwoState(kAB=kp,
                                        kBA=km,
                                        k1=kET1EPT1,
                                        k2=kEPT2,
                                        cA0=c0[0],
                                        cB0=c0[1],
                                        cR0=PS0,
                                        )

t, ct, PSt = system.evolve(t_tot=200)

plt.plot(t, ct[0]/ct[0][0], 'k-', label=r'${\rm TrpNH\cdot\!\cdot\!\cdot H_2O}$', lw=1.5)
plt.plot(t, ct[1]/ct[1][0], 'r-', label=r'${\rm TrpNH\cdot\!\cdot\!\cdot OH^-}$', lw=1.5)
plt.legend(fontsize=17)

plt.xlim(0,200)
plt.ylim(0,1.1)
plt.xlabel(r't / $\rm \mu s$', fontsize=17)
plt.ylabel(r'${[{\rm Red}_i]}/{[{\rm Red}_i]_{t=0}}$', fontsize=17)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()

plt.show()



#=========================================================
# Example 1.2: kinetics from different models
#=========================================================

def single_exp(t, k):
    return np.exp(-k*t)

# apparent rate constant calculated using the kinetic model with equilibrium assumption
k_app = 1/(1+KOH*cOH)*kET1EPT1 + KOH*cOH/(1+KOH*cOH)*kEPT2

# fit the time evolution of PS to single exponential decay
k_fit = curve_fit(single_exp, t, PSt/PS0)[0][0]/cWEE

print('Apparent rate constant obtained from kinetic model: %f M^-1 us^-1'%k_app)
print('Apparent rate constant obtained from numerical integration + single exp fit: %f M^-1 us^-1'%(k_fit))

plt.plot(t, PSt/PS0, 'k-', label='numerical integration', lw=1.5)
plt.plot(t, np.exp(-k_fit*cWEE*t), 'r-', label='single exp fit', lw=1.5)
plt.plot(t, np.exp(-k_app*cWEE*t), 'b-', label='kinetic model', lw=1.5)

plt.legend(fontsize=17)
plt.xlim(0,200)
plt.ylim(0,1.1)
plt.xlabel(r't / $\rm \mu s$', fontsize=17)
plt.ylabel(r'${[{\rm PS^{\!+\!}}]}/{[{\rm PS^{\!+\!}}]_{t=0}}$', fontsize=17)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()

plt.show()
