function op = integrator_PG2e_setup(M, Z, DT)
%INTEGRATOR_PG2E_SETUP  Build reusable operators for the PG2e time integrator.
%
%   OP = INTEGRATOR_PG2E_SETUP(M, Z, DT) precomputes once, every step-independent
%   quantity needed by INTEGRATOR_PG2E.
%
%   The method is presented in mass-normalized form with
%       p = 3I + DT*zeta,  q = 2I + DT*zeta,  zeta = M^{-1} Z.
%   For implementation this routine uses the EQUIVALENT non-normalized form
%       P = 3M + DT*Z,     Q = 2M + DT*Z,
%   obtained by multiplying the mass-normalized predictor and corrector equations
%   through by M. The two forms are algebraically identical. The non-normalized
%   form keeps P and Q symmetric whenever M and Z are symmetric, and keeps them
%   banded/sparse for banded/sparse models, so a structure-aware factorization
%   can be used.
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
%            M       the mass matrix as supplied (vector or matrix), for forming M*x
%            Mmul    function handle applying M (i.e. x -> M*x)
%            Z       the damping matrix (sparse), for forming Z*x
%            dP, dQ  factorizations of P = 3M + DT*Z and Q = 2M + DT*Z
%            sym     true if P and Q are symmetric
%          Pass OP to INTEGRATOR_PG2E.
%
%   See also INTEGRATOR_PG2E.

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
op.Z  = Z;

% --- mass application x -> M*x, and a sparse mass matrix for forming P, Q ---
if isvector(M)                          % lumped mass supplied as a vector
    md      = M(:);
    op.M    = md;
    op.Mmul = @(x) md .* x;             % M*x component-wise
    Mmat    = spdiags(md, 0, N, N);     % diagonal mass as a sparse matrix
elseif isdiag(M)                        % diagonal mass matrix
    md      = full(diag(M));
    op.M    = md;
    op.Mmul = @(x) md .* x;
    Mmat    = spdiags(md, 0, N, N);
else                                    % consistent (non-diagonal) mass matrix
    op.M    = M;
    op.Mmul = @(x) M * x;
    Mmat    = M;
end

% --- non-normalized operators P, Q: factorize ONCE, reuse every step ---
% P = 3M + DT*Z,  Q = 2M + DT*Z. Symmetric whenever M and Z are symmetric.
P = 3*Mmat + DT*Z;
Q = 2*Mmat + DT*Z;

% Symmetrize numerically tiny asymmetry from assembly so the automatic
% decomposition can take a symmetric (banded / Cholesky / LDL) path.
isSym  = @(A) (nnz(A) == 0) || (normest(A - A.', 1) <= 1e-12 * max(1, normest(A, 1)));
op.sym = isSym(P) && isSym(Q);
if op.sym
    P = (P + P.')/2;
    Q = (Q + Q.')/2;
end

% Automatic type selection (banded / Cholesky / LDL / LU as appropriate).
op.dP = decomposition(P);
op.dQ = decomposition(Q);

end
