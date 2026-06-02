function op = integrator_PG2e_setup(M, Z, DT)
%INTEGRATOR_PG2E_SETUP  Build operators for the PG2e time integrator.
%
%   Inputs
%     M    mass matrix:  N x N (full or sparse), or an N x 1 vector of lumped masses
%     Z    N x N damping matrix  (use [] or all-zeros for the undamped case)
%     DT   scalar time step
%
%   Output
%     OP   struct with fields:
%            DT      the time step
%            N       number of degrees of freedom
%            Msolve  function handle applying M^{-1}
%            zeta    mass-normalized damping  M^{-1} Z
%            dP, dQ  factorizations of p = 3I + DT*zeta and q = 2I + DT*zeta

if nargin < 3 || isempty(DT)
    error('integrator_PG2e_setup:DT', 'A scalar time step DT is required.');
end

if isvector(M)
    N = numel(M);
else
    N = size(M, 1);
end
if isempty(Z)
    Z = sparse(N, N);
end

op.DT = DT;
op.N  = N;

% --- efficient M^{-1} application and zeta = M^{-1} Z ---
if isvector(M)                          % lumped mass supplied as a vector
    md = M(:);
    op.Msolve = @(x) x ./ md;
    zeta = (1 ./ md) .* Z;              % row scaling (implicit expansion), stays sparse
elseif isdiag(M)                        % diagonal mass matrix
    md = full(diag(M));
    op.Msolve = @(x) x ./ md;
    zeta = (1 ./ md) .* Z;
else                                    % consistent (non-diagonal) mass matrix
    dM = decomposition(M, 'chol');      % SPD mass matrix: factorize once
    op.Msolve = @(x) dM \ x;
    zeta = dM \ Z;                      % dense if M and Z are non-diagonal
    warning('integrator_PG2e_setup:consistentMass', ...
        ['Consistent mass matrix detected: M^{-1}Z may be dense. ', ...
         'Lumped (diagonal) mass is recommended for large models.']);
end
op.zeta = zeta;

% --- damping operators p, q: factorize once, reuse every step ---
I = speye(N);
op.dP = decomposition(3*I + DT*zeta);   % p = 3I + DT*zeta
op.dQ = decomposition(2*I + DT*zeta);   % q = 2I + DT*zeta

end
