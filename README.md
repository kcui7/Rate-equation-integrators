# Rate-equation-generators
This repository collects various kinetic models which numerically integrate the rate equations. These kinetic models are highly system-specific. Please refer to the following list and the header of the corresponding codes for the applicable scenarios. 

## Code list
### rate_equation_integrator_two_state.py
This code is used to numerically solve the rate equations for the following system:

Two molecules can convert to each other, 

    A <==> B
    
This reaction is characterized by the forward rate constant `kAB`, backward rate constant `kBA`, and equilibrium constant K.
Both A and B can react with another molecule R, forming the same product P, these reactions are IRREVERSIBLE:

    A + R --> P  k1
    B + R --> P  k2

Input parameters are:

1. `kAB` (float): rate constant for reaction A --> B

2. `kBA` (float): rate constant for reaction B --> A

3. `k1` (float): rate constant for reaction A + R --> P

4. `k2` (float): rate constant for reaction A + R --> P

5. `cA0` (float): initial concentration of A (in mol/L)

6. `cB0` (float): initial concentration of B (in mol/L)
 
7. `cR0` (float): initial concentration of R (in mol/L)

To create an instance and integrate `t_tot` time with time step `dt`, write
```python
import RateEquationIntegratorTwoState

system = RateEquationIntegratorTwoState(kAB=kAB,
                                        kBA=kBA,
                                        k1=k1,
                                        k2=k2,
                                        cA0=cA0,
                                        cB0=cB0,
                                        cR0=cR0,
                                        )

t, cABt, cRt = system.evolve(t_tot=200, dt=1e-3)

```
Due to the implementation of the integrator, `cABt` combines the concentrations cAt and cBt, 
```python
cAt = cABt[0]
cBt = cABt[1]
```
Please refer to **Ref. 1** for the detailed derivation and application of this model. An example is included in `rate_equation_integrator_two_state.py`

### rate_equation_integrator_water_autoprotolysis.py
This class is used to numerically solve the rate equations for the following system:

Two molecules can convert to each other, 

        A + OH- <==> B
        
This reaction is characterized by the forward rate constant `kAB`, backward rate constant `kBA`, and equilibrium constant K. OH- comes from water autoprotolysis,

        H2O <==> H+ + OH-
        
The forward rate constant `kauto = 2.6e-5*55.5 s^-1` and backward rate constant `kneu = 1.3e11 M^-1s^-1`. 
Both A and B can react with another molecule R, forming the same product P, these reactions are IRREVERSIBLE,

        A + R --> P  k1
        B + R --> P  k2

Input parameters are:

1. `kAB` (float): rate constant for reaction A + OH- --> B

2. `kBA` (float): rate constant for reaction B --> A + OH-

3. `kauto` (float): rate constant for water autoprotolysis

4. `kneu` (float): rate constant for neutralization

5. `k1` (float): rate constant for reaction A + R --> P

6. `k2` (float): rate constant for reaction A + R --> P

7. `cA0` (float): initial concentration of A (in mol/L)

8. `cB0` (float): initial concentration of B (in mol/L)

9. `cOH0` (float): initial concentration of OH- (in mol/L)

10. `cR0` (float): initial concentration of R (in mol/L)

To create an instance and integrate `t_tot` time with time step `dt`, write
```python
import RateEquationIntegratorTwoState

system = RateEquationIntegratorTwoState(kAB=kAB,
                                        kBA=kBA,
                                        kauto=kauto,
                                        kneu=kneu,
                                        k1=k1,
                                        k2=k2,
                                        cA0=cA0,
                                        cB0=cB0,
                                        cH0=cH0,
                                        cOH0=cOH0,
                                        cR0=cR0,
                                        )

t, cABt, cRt, cHt, cOHt = system.evolve(t_tot=200, dt=1e-3)
```
Please refer to **Ref. 1** for the detailed derivation and application of this model. An example is included in `rate_equation_integrator_water_autoprotolysis.py`

## Citation
If you find these kinetic models helpful, please cite the corresponding papers: 
1. Cui, K.; Soudackov, A. V.; Hammes-Schiffer, S. Modeling the weak pH Dependence of Proton-Coupled Electron Transfer for Tryptophan Derivatives, *J. Phys. Chem. Lett.* In Press. DOI:[10.1021/acs.jpclett.3c02282](https://doi.org/10.1021/acs.jpclett.3c02282)
