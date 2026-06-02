# Two-node explicit Petrov–Galerkin time integrator

Developed by Jae Young Lee, Department of Mechanical Engineering, Ajou University.

## Integrators

  - `integrator_PG2e_setup.m` — builds the reusable operator once (prefactorization).
  - `integrator_PG2e.m` — advances the solution one step, using the precomputed operators.

## Usage

The code is a two-call **setup → step** pattern (see the files for full docs):

```matlab
op             = integrator_PG2e_setup(M, Z, DT);	% once, before the loop
[U_Tp1, V_Tp1] = integrator_PG2e(op, F_int, U_T, V_T, R_T, R_Tp1);	% each step
```

`integrator_PG2e_setup` prefactorizes the damping operators once and returns an operator struct `op`.
`integrator_PG2e` then advances one step per call. No matrix is formed or factorized inside the loop. 
The internal force is supplied as a handle `F_int(U,V)`.

An example script for a linear 2-DOF system in free vibration:

```matlab
addpath('integrators/PG2e');

% --- Structural system:  M a + Z v + K u = R(t)  (lumped/diagonal M recommended) ---
M = [2 0; 0 1];  K = [6 -2; -2 4];  Z = [0.2 0; 0 0.1];

DT = 1e-3;

% --- Build reusable operators once (prefactorization) ---
op    = integrator_PG2e_setup(M, Z, DT);

% --- Internal force handle:  F_int(U,V) = K*U  for a linear system ---
F_int = @(U,V) K*U;

% --- Initial conditions ---
U = [0; 0];  V = [1; 0];

% --- External load R(t)  (here: free vibration, R = 0) ---
Rfun = @(t) [0; 0];

% --- Time integration loop (the integrator advances one step per call) ---
T = 5;  nT = round(T/DT);  t = 0;
for n = 1:nT
    [U, V] = integrator_PG2e(op, F_int, U, V, Rfun(t), Rfun(t+DT));
    t = t + DT;
end

disp(U);   % displacement at final time
disp(V);   % velocity at final time
```

## Requirements

- MATLAB R2017b or later. The setup routine uses the `decomposition` object to prefactorize the damping operators.
- No additional toolboxes are required. 
- The code uses only base MATLAB (sparse matrix operations, the backslash solver, and function handles).

## Citation

If you use this code, please cite:

> J. Y. Lee, An explicit Petrov–Galerkin time finite element integrator for structural dynamics with physical damping,
> *Computers & Structures*, submitted.

## License

MIT License — see [LICENSE](LICENSE).
