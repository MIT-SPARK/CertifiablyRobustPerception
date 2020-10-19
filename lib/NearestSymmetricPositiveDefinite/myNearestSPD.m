function [Ahat,H] = myNearestSPD(A,nrEigToFix)
% nearestSPD - the nearest (in Frobenius norm) Symmetric Positive Definite matrix to A
% usage: Ahat = nearestSPD(A)
%
% From Higham: "The nearest symmetric positive semidefinite matrix in the
% Frobenius norm to an arbitrary real matrix A is shown to be (B + H)/2,
% where H is the symmetric polar factor of B=(A + A')/2."
%
% http://www.sciencedirect.com/science/article/pii/0024379588902236
%
% arguments: (input)
%  A - square matrix, which will be converted to the nearest Symmetric
%    Positive Definite Matrix.
%
% Arguments: (output)
%  Ahat - The matrix chosen as the nearest SPD matrix to A.

% test for a square matrix A
[r,c] = size(A);
if r ~= c
  error('A must be a square matrix.')
elseif (r == 1) && (A <= 0)
  % A was scalar and non-positive, so just return eps
  Ahat = eps;
  return
end

% symmetrize A into B
B = (A + A')/2;

[Ve,De] = eig(B);
indNegEig = find(diag(De)<1e-10);

if nargin ==2
    de = diag(De);
    [sortedDe,indSorted] = sort(de);
    Ve = Ve(:,indSorted);
    De = De(indSorted,indSorted);
    indNegEig = [1:nrEigToFix];% fix nrEigToFix smallest eig
end

H = Ve(:,indNegEig) * abs(De(indNegEig,indNegEig)) *Ve(:,indNegEig)'; 
Ahat = B + H;

% ensure symmetry
Ahat = (Ahat + Ahat')/2;

% [sort(eig(B)) sort(eig(Ahat))]

% test that Ahat is in fact PD. if it is not so, then tweak it just a bit.
p = 1;
k = 0;
while p ~= 0
  [R,p] = chol(Ahat);
  k = k + 1;
  if p ~= 0
    % Ahat failed the chol test. It must have been just a hair off,
    % due to floating point trash, so it is simplest now just to
    % tweak by adding a tiny multiple of an identity matrix.
    mineig = min(eig(Ahat));
    Ahat = Ahat + (-mineig*k.^2 + eps(mineig))*eye(size(A));
  end
end
