# structural-time-integrators

An explicit two-node Petrov–Galerkin time finite element integrator for structural dynamics with physical damping. 
The scheme advances the second-order equation of motion `M a + Z v + K u = R(t)` one step at a time using a single-step predictor–corrector form
that is explicit with respect to the stiffness matrix — the only linear solves involve the (typically small) damping operator. 
In the undamped limit it reduces to velocity-Verlet.

Developed by Jae Young Lee, SunLab, Ajou University.

## Integrators

- **PG2e** — two-node explicit Petrov–Galerkin time integrator:
  - `integrators/PG2e/integrator_PG2e_setup.m` — builds the reusable operator struct once (prefactorizes the damping operators).
  - `integrators/PG2e/integrator_PG2e.m` — advances the solution one step, using the precomputed operators.

## Usage

The API is a two-call **setup → step** pattern (see the files for full docs):

```matlab
op             = integrator_PG2e_setup(M, Z, DT)             % once, before the loop
[U_Tp1, V_Tp1] = integrator_PG2e(op, F_int, U_T, V_T, R_T, R_Tp1)   % each step
```

`integrator_PG2e_setup` prefactorizes the damping operators once and returns an operator
struct `op`; `integrator_PG2e` then advances one step per call, with no matrix formed or
factorized inside the loop. The internal force is supplied as a handle `F_int(U, V)`.

A minimal driver for a linear 2-DOF system in free vibration:

```matlab
addpath('integrators/PG2e');

% --- Structural system:  M a + Z v + K u = R(t)  (lumped/diagonal M recommended) ---
M = [2 0; 0 1];  K = [6 -2; -2 4];  Z = [0.2 0; 0 0.1];

DT = 1e-3;

% --- Build reusable operators once (prefactorization; requires R2017b+ for decomposition) ---
op    = integrator_PG2e_setup(M, Z, DT);

% --- Internal (elastic) force handle:  F_int(U,V) = K*U  for a linear system (K captured) ---
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

For an undamped system, pass `Z = []` (or an all-zeros matrix) to `integrator_PG2e_setup`;
the formulas reduce to velocity-Verlet with no special case. For nonlinear or
velocity-dependent internal forces, supply an `F_int` handle that uses both `U` and `V`; the
right-endpoint force is evaluated with the predicted velocity.

## Requirements

- MATLAB R2017b or later. The setup routine uses the `decomposition` object (introduced in
  R2017b) to prefactorize the damping operators. No additional toolboxes required — the code
  uses only base MATLAB (sparse matrix operations, the backslash solver, and function handles).

## Citation

If you use this code, please cite:

> J. Y. Lee, "An explicit Petrov–Galerkin time finite element integrator for structural dynamics with physical damping,"
> *Computers & Structures*, submitted.

```bibtex
% TODO: update with final publication details (volume, pages, year, doi) once available.
@article{Lee2026,
  author  = {Lee, Jae Young},
  title   = {An explicit Petrov--Galerkin time finite element integrator for structural dynamics with physical damping},
  journal = {Computers \& Structures},
  year    = {},
  volume  = {},
  pages   = {},
  doi     = {}
}
```

## License

MIT License — see [LICENSE](LICENSE).
