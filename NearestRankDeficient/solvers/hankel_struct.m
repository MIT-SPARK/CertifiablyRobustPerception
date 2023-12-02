% Constructs the affine map of Hankel structure
% SS: R^k -> R^{m x n},
%    u -> hankel(u)
% S = Scell{1}; Scell{2}; ... Scell{k}; Scell{0}

function [S,k,Scell] = hankel_struct(m,n)

k = m+n-1;
Im = eye(m);

SS = zeros(m,k+1,n);
S = zeros(m*(k+1),n);
for i=1:n
    SS(:,i:i+m-1,i) = Im;
    S(:,i) = vect(SS(:,:,i));
end

S   = sparse(S);

if nargout > 2
    Scell = {};
    for i = 1:k+1
        Scell{end+1} = sparse(squeeze(SS(:,i,:)));
    end
end

end

function v = vect(M)
v = M(:);

end
