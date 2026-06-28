function [U_Tp1, V_Tp1] = integrator_PG2e(op, F_int, U_T, V_T, R_T, R_Tp1)
%INTEGRATOR_PG2E  One step of the explicit two-node Petrov-Galerkin time integrator.
%
%   [U1, V1] = INTEGRATOR_PG2E(OP, F_INT, U, V, R, R1) advances 
%   the second-order system   M*A + Z*V + F_int(U,V) = R   by one time step.
%     predictor      P * Vp = 3*M*V - 2*dt*Z*V + 3*dt*(R - F)
%     displacement   U1     = U + (dt/2)*(V + Vp)
%     corrector      Q * V1 = 2*M*V -   dt*Z*V +   dt*((R - F) + (R1 - F1))
%   with P = 3M + dt*Z and Q = 2M + dt*Z prefactorized once by INTEGRATOR_PG2E_SETUP.
%   The right-endpoint internal force is evaluated with the predicted velocity Vp,
%   which matters for velocity-dependent internal forces.
%   The undamped limit Z -> 0 reduces to velocity Verlet.
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
%   See also INTEGRATOR_PG2E_SETUP.

dt = op.DT;

RmF_T   = R_T - F_int(U_T, V_T);                          % R - F at t (non-normalized)
MV      = op.Mmul(V_T);                                   % M * V_T
ZV      = op.Z * V_T;                                     % Z * V_T (mat-vec)

V_Tpp   = op.dP \ (3*MV - 2*dt*ZV + 3*dt*RmF_T);          % predictor (reuses factor)
U_Tp1   = U_T + 0.5*dt*(V_T + V_Tpp);                     % displacement

RmF_Tp1 = R_Tp1 - F_int(U_Tp1, V_Tpp);                    % R - F at t+dt (predicted V)
V_Tp1   = op.dQ \ (2*MV - dt*ZV + dt*(RmF_T + RmF_Tp1));  % corrector (reuses factor)

end
