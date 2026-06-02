function [U_Tp1, V_Tp1] = integrator_PG2e(op, F_int, U_T, V_T, R_T, R_Tp1)
%INTEGRATOR_PG2E  an explicit two-node Petrov-Galerkin time integrator.
%
%   Inputs
%     OP     operator struct returned by INTEGRATOR_PG2E_SETUP(M, Z, dt)
%     F_INT  function handle  F_int(U, V) -> N x 1 internal (elastic) force.
%              Linear:    F_int = @(U,V) K*U;
%              Nonlinear: F_int = @(U,V) myInternalForce(U,V);
%            (K and any other data are captured inside the handle.)
%     U_T    N x 1 displacement at time t
%     V_T    N x 1 velocity     at time t
%     R_T    N x 1 external force at time t
%     R_Tp1  N x 1 external force at time t+dt
%
%   Outputs
%     U_Tp1  N x 1 displacement at time t+dt
%     V_Tp1  N x 1 velocity     at time t+dt
%
%   Update (mass-normalized;  r = M^{-1}R,  f = M^{-1}F_int,  zeta = M^{-1}Z):
%     predictor      p * Vp = (3I - 2*dt*zeta)*V + 3*dt*(r - f)
%     displacement   U1     = U + (dt/2)*(V + Vp)
%     corrector      q * V1 = (2I -   dt*zeta)*V +   dt*((r - f) + (r1 - f1))

dt = op.DT;

RmF_T   = op.Msolve(R_T   - F_int(U_T,   V_T));           % r - f at t
zV      = op.zeta * V_T;                                  % zeta * V_T (mat-vec)

V_Tpp   = op.dP \ (3*V_T - 2*dt*zV + 3*dt*RmF_T);         % predictor (reuses factor)
U_Tp1   = U_T + 0.5*dt*(V_T + V_Tpp);                     % displacement

RmF_Tp1 = op.Msolve(R_Tp1 - F_int(U_Tp1, V_Tpp));         % r - f at t+dt (predicted V)
V_Tp1   = op.dQ \ (2*V_T - dt*zV + dt*(RmF_T + RmF_Tp1)); % corrector (reuses factor)

end
