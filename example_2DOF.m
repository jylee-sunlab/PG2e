%example_2DOF  An example script for a linear 2-DOF system in free vibration

% --- Structural system:  M a + Z v + K u = R(t) ---
M = [2 0; 0 1];  K = [6 -2; -2 4];  Z = [0.2 0; 0 0.1];

DT = 1e-3;

% --- Build reusable operators once (pre-factorization) ---
op    = integrator_PG2e_setup(M, Z, DT);

% --- Internal force handle:  F_int(U,V) = K*U  for a linear system ---
F_int = @(U,V) K*U;

% --- Initial conditions ---
U = [0; 0];  V = [1; 0];

% --- External load R(t) (free vibration, R = 0) ---
Rfun = @(t) [0; 0];

% --- Time integration loop (the integrator advances one step per call) ---
T = 5;  nT = round(T/DT);  t = 0;
for n = 1:nT
    [U, V] = integrator_PG2e(op, F_int, U, V, Rfun(t), Rfun(t+DT));
    t = t + DT;
end

disp(U);   % displacement at final time
disp(V);   % velocity at final time

disp(U);   % displacement at final time
disp(V);   % velocity at final time

% Expected (approx, DT = 1e-3):
%   U = [0.1451; 0.4877]
%   V = [0.4076; 0.2495]

