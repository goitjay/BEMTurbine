# BEMTurbine

This is a simple blade element momentum (BEM) analysis tool for horizontal axis wind turbines. 

The tool consists of two main python scripts:

1. optimum_rotor.py: Optimizes the blade shape (twist and chord) to the 'Betz
optimum' blade. The input parameters for this script should be defined in
setup/optimum_rotor.inp.

2. BEM_analysis.py: Estimates axial induction factor and angular induction factor
iteratively using blade element momentum (BEM) theory. It can perform two computations.
  (a)estimate thrust, torque and power coefficients as a function of tip-speed ratio. The script also
generates axial induction factor, angular induction factor, angle of attack, Cl
and Cd for each element for the given tip speed ratio. The input parameters for
this script should be defined in setup/BEM_analysis.inp.
  (b) compute power, thrust and torque as a function of wind speed.

All the ouputs are save in folder output.

The example case consist of
(a) KDWT25: optimization of blades for a small wind turbine with
a rotor diameter of 0.25 m. The turbine uses SD7003 airfoil from the root to tip
of the blade.

(b) KDWT300kW: 300 kW wind turbine with rotor diameter of 33 m. 
