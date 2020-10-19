function Ahat = nearestSPD(A)
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

if nargin ~= 1
  error('Exactly one argument must be provided.')
end

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

% Compute the symmetric polar factor of B. Call it H.
% Clearly H is itself SPD.
[U,Sigma,V] = svd(B);
H = V*Sigma*V';

% get Ahat in the above formula
Ahat = (B+H)/2;

% ensure symmetry
Ahat = (Ahat + Ahat')/2;

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

%% %% test by Luca:
% %% 1- find negative eigenvalues and corresponding U/V columns
% [Ve,De] = eig(B);
% dDe= diag(De);
% dSigma = diag(Sigma);
% indNegSigmas = [];
% for i=1:r
%    if dDe(i)<0
%        for j=1:r
%           if abs(dSigma(j) + dDe(i))<1e-5
%             indNegSigmas = [indNegSigmas j];  
%           end
%        end
%    end
% end
% 
% %% 2- shift negative eigenvalues to make them 0
% % norm(U*Sigma*V' - B) = 0
% AhatLC = B - 0.5 * U(:,indNegSigmas)*Sigma(indNegSigmas,indNegSigmas)*V(:,indNegSigmas)'...
%     +0.5 * V(:,indNegSigmas)*Sigma(indNegSigmas,indNegSigmas)*V(:,indNegSigmas)';
% if norm(AhatLC - Ahat)>1e-10
%    error('AhatLC ~= Ahat') 
% end
% 
% %% try even simpler version:
% % projection is just taking left singular vectors and replacing them 
% % with right singular vectors
% AhatLC2 = B + 0.5*( V(:,indNegSigmas)-U(:,indNegSigmas) )...
%     * Sigma(indNegSigmas,indNegSigmas)*V(:,indNegSigmas)' ;
% if norm(AhatLC2 - Ahat)>1e-10
%    error('AhatLC2 ~= Ahat') 
% end
% 
% if norm(AhatLC2 - Ahat)>1e-10
%    error('AhatLC2 ~= Ahat') 
% end
% 
% %% last test for fun:
% indPosSigmas = setdiff([1:r],indNegSigmas);
% if norm(V(:,indPosSigmas)-U(:,indPosSigmas)) > 1e-10
%    disp('V(:,indPosSigmas) ~= U(:,indPosSigmas)') 
%    % sometimes the sign of a single colum is flipped, interesting
% end






