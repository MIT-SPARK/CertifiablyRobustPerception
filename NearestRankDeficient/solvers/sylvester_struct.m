% Constructs the affine map of Sylvester structure
% SS: R^D -> R^{m x n},
%    (u1,u2) -> Syl(u1,u2)

function [S,k,m,n] = sylvester_struct(D,d)
k = sum(D+1);
m = sum(D)-2*d+2;
n = sum(D)-d+1;

S = zeros(m*(k+1),n);
for i=1:sum(D)-d+1
    Pi = get_matrixGCD_i(D,d,i);
    Pi = [Pi zeros(m,1)];
    S(:,i) = Pi(:);
end

function Pi = get_matrixGCD_i(D,d,i)
k1 = D(2)-d+1;
k2 = D(1)-d+1;
P1 = (eye(D(1)+1, D(1)+1));
P1 = [zeros(D(2)-d,D(1)+1); P1; zeros(D(2)-d,D(1)+1)];
P2 = P1(end-i+2-k1:end-i+1,:);
Q1 = (eye(D(2)+1, D(2)+1));
Q1 = [zeros(D(1)-d,D(2)+1); Q1; zeros(D(1)-d,D(2)+1)];
Q2 = Q1(end-i+2-k2:end-i+1,:);
Pi = blkdiag(P2,Q2);
