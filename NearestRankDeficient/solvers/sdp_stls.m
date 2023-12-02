% SDP relaxation for structured total least squares (STLS)
%
% Let k,m,n integers (m<=n)
% Let SS: R^k -> R^{m x n} be an affine map
% Let u1 be a vector in R^k
% Consider the optimization problem
%
% min_u     |u - u1|^2
% s.t.      SS(u) is rank deficient
%
% More generally, the cost can be a quadratic function
%           (u - u1)' * W * (u - u1)
% for some symmetric weight matrix W of size k x k.
% If W is diagonal and the entry W_ii is equal to zero then
% the respective u1(i) is ignored (as if it was missing).
%
% The function also accepts problems with complex data.
% They are internally converted to real valued problems.
%
% Input:
% S - matrix of size (k+1)m x n describing the affine map SS
% u1 - a vector of length k, possibly containing 'nan' entries
%
% Optional inputs:
% W - weight matrix of size k x k (default W_ii=1 except for 'nan' entries)
% solver - SDP solver (default 'sdpt3')
% quiet - set to false to display information (default true)
%
% Output:
% opt - optimal value of SDP relaxation
% u - the minimizer of the problem
% U - the matrix SS(u)
% z - a vector in the left kernel of SS(u)
% X - the PSD matrix (relaxation is exact if rank(X)=1)
%
% Usage:
% [opt,u,U,z,X] = sdp_stls(S,u1,W,solver,quiet)
% The optional arguments can be omitted or set to empty values.
% For instance: sdp_stls(S,u1,[],solver)

function [opt,u,U,z,X] = sdp_stls(S,u1,W,solver,quiet)

narginchk(2,5);
if nargin<3||isempty(W); W = diag(~isnan(u1)); end
if nargin<4||isempty(solver); solver = 'sdpt3'; end
if nargin<5||isempty(quiet); quiet = false; end

k = length(u1);
n = size(S,2);
m = size(S,1)/(k+1);
if floor(m)~=m; error('mismatch in dimensions'); end
u1(isnan(u1)) = 0;
u1 = reshape(u1,[k,1]);

iscomplex = ~(isreal(S) && isreal(u1) && isreal(W));
if iscomplex; [k,m,n,S,u1,W] = data2real(k,m,n,S,u1,W); end
W = .5*(W+W');

% fprintf('PSD matrix size: %d\n',(k+1)*m)

g = [eye(k) -u1];
G0 = g'*W*g;
E0 = sparse(k+1,k+1,1);

Im = eye(m);
G = kron(G0,Im);
E = kron(E0,Im);

[opt, X] = primal_cvx(k,m,n,S,G,E,solver,quiet);

[~,E] = eig(full(X));
r=nnz(diag(E)>1e-4);
if (~iscomplex && r>1) || (iscomplex && r>2)
    warning('sdp might not be exact');
end

J0 = k*m+1:(k+1)*m;
z = recoverSol(X(J0,J0));
u = zeros(k,1);
for i = 1:k
    Ji = (i-1)*m+1:i*m;
    u(i) = trace(X(J0,Ji));
end
U = applyAffineMap(S,u);

if iscomplex
    [u,U,z] = data2complex(k,m,n,u,U,z);
end

function [k,m,n,S,u1,W] = data2real(k,m,n,S,u1,W)
toReal = @(a) [real(a) -imag(a); imag(a) real(a)];
u1 = [real(u1); imag(u1)];
if norm(W-W')>1e-4; warning('W forced to be hermitian'); end
W = toReal(W);
SS = mat2cell(S,m*ones(k+1,1),n);
B = SS(1:k); A = SS(k+1);
iB = cellfun(@(x) 1i*x, B,'Unif',0);
SS = cellfun(toReal, [B; iB; A],'Unif',0);
S = cell2mat(SS);
k = 2*k; m = 2*m; n = 2*n;

function [u,U,z] = data2complex(k,m,n,u,U,z)
k = k/2; m = m/2; n = n/2;
u = u(1:k) + 1i * u(k+1:2*k);
U = U(1:m,1:n) + 1i * U(m+1:2*m,1:n);
z = z(1:m) + 1i * z(m+1:2*m);

% dual sdp
function [opt, X] = dual_cvx(k,m,n,S,G,E,solver,quiet)
N = (k+1)*m;
SC = cell(n,1);

if quiet
cvx_begin sdp quiet
else
cvx_begin sdp
end
    cvx_solver(solver)
    variable t(1,1);
    variable Cvec(N,n)
    dual variable X
    for i=1:n
        Ci = Cvec(:,i);
        si = S(:,i);
        SC{i} = si*Ci';
    end
    sSC = sum(cat(3,SC{:}),3);
    sSC = blksym(k+1,m,sSC);
    maximize(t);
    Q = G - t*E + sSC;
    X: Q >= 0;
cvx_end

opt = cvx_optval;

% primal sdp
function [opt, X] = primal_cvx(k,m,n,S,G,E,solver,quiet)
N = (k+1)*m;

e = speye(N);
SC = zeros(N,N,n*N);
for i=1:n
    for j=1:N
        SC(:,:,(i-1)*N+j) = S(:,i)*e(j,:);
    end
end
SC = blksym(k+1,m,SC);
A = zeros(n*N,N*(N+1)/2);
for l=1:n*N
    A(l,:) = smat2vec(SC(:,:,l));
end
A = sparse(A);

if quiet
cvx_begin sdp quiet
else
cvx_begin sdp
end
    cvx_solver(solver)
    variable X(N,N) symmetric
    dual variable Q
    y = smat2vec(X);
    minimize(smat2vec(G)'*y);
    smat2vec(E)'*y == 1;
    A*y == 0;
    Q: X >= 0;
cvx_end

opt = cvx_optval;

% recover minimizer from moment matrix
function [x,e] = recoverSol(X)
N = size(X,1);
if any(isnan(X))
    x = nan(N,1);
    e = inf;
else
    [V,E] = eig(full(X));
    e=diag(E);
    x = sqrt(e(N))*V(:,N);
    e = e(N-1);
end

% vector to symmetric matrix
function M = vec2smat(v)

N = length(v);
k = (-1+sqrt(1+8*N))/2;

M = repmat(0*v(1:k),[1,k]);
I = triu(true(k,k),0);
I2 = triu(true(k,k),1);
M(I) = v/2;
M(I2) = M(I2)*sqrt(2);
M = M + M.';

% symmetric matrix to vector
function v = smat2vec(M)

k = size(M,1);
I = triu(true(k,k),0);
I2 = triu(true(k,k),1);
M(I2) = M(I2)*sqrt(2);
v = M(I);

% Block symmetrization of a matrix
% Input: mk x mk matrix A
% Output: symmetric matrix As such that
%         each m x m block is also symmetric (there are n^2 blocks)
function As = blksym(n,m,A)

As = mysym(A);
if isempty(A); return; end

for i=1:n
    for j=1:n
        I = m*(i-1) + (1:m);
        J = m*(j-1) + (1:m);
        As(I,J,:) = mysym(As(I,J,:));
    end
end

% symmetrize matrix
function As = mysym(A)
At = permute(A,[[2,1],3:ndims(A)]);
As = .5*(A+At);
