# Numerical Simulation and Analysis of Mutual Inductance Circuits

A numerical analysis and simulation study of a non-linear mutual inductance RLC circuit, implemented in MATLAB and independently verified with LTspice.

## Project Overview
This project models a coupled RLC circuit described by a system of three state variables: primary current (i_1), secondary current (i_2), and capacitor voltage (u_c). The study evaluates numerical methods across convergence, computational cost, and numerical stability under both linear and non-linear mutual inductance (M_n) models.

## Implemented Numerical Methods
* **ODE Solvers:** Euler (O(h)), Heun (O(h^2)), and Fourth-Order Runge-Kutta (RK4, O(h^4)).
* **Interpolation and Approximation:** Lagrange interpolation, Natural Cubic Spline, and Least Squares Polynomial Fitting (n=3, 5).
* **Numerical Integration:** Composite Rectangular Rule (O(h)) and Composite Trapezoidal Rule (O(h^2)) applied to compute total dissipated energy (P).
* **Root-Finding Algorithms:** Bisection, Secant, and Quasi-Newton methods applied to determine resonant excitation frequencies.

## Verification and Validation
The algorithms were independently benchmarked against two standard baselines:
* **MATLAB Built-in `ode45`:** Achieved agreement with the hand-written RK4 solver with a maximum relative error below 0.0001%.
* **LTspice Simulation:** Energy calculations matched LTspice's trapezoidal solver within 0.86%, confirming the physical consistency of the differential equation system.

## Repository Structure
```text
├── LTspice/                                    # LTspice schematics and simulation files
├── matlab codes/                               # Source code for ODE solvers, interpolation, and root-finding
├── Project – mutual inductance circuit report.pdf # Full technical documentation and results
└── README.md
