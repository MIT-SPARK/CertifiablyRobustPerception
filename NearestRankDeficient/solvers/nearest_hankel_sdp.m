function SDP = nearest_hankel_sdp(S,u1,W)
%% The nearest rank deficient problem
%% min_u  || u - u1 ||
%% s.t.  S(u) is rank deficient
%% Implement the SDP relaxation in Diego Cifuentes
%% A CONVEX RELAXATION TO COMPUTE THE NEAREST STRUCTURED RANK DEFICIENT MATRIX
%% Possibly allow for weighted norm by W

k       = length(u1);
if nargin < 3
    W   = eye(k);
end
W       = sparse(0.5*(W+W')); % make sure W is symmetric
n       = size(S,2);
m       = size(S,1)/(k+1);
if floor(m)~=m; error('mismatch in dimensions'); end

g       = [speye(k) -u1];
G0      = g'*W*g;
E0      = sparse(k+1,k+1,1,k+1,k+1);

Im      = speye(m);
G       = kron(G0,Im);
E       = kron(E0,Im);

N       = (k+1)*m;
blk{1,1}= 's';
blk{1,2}= N;
SDP.n   = N; 
SDP.blk = blk;
SDP.C   = G;

%% construct At and b for SDP
Acell   = {E};
% Construct rank-deficient constraints
e       = speye(N);
fprintf('Build rank-deficient constraints, progress ...')
for i=1:n
    fprintf('%d/%d ',i,n);
    for j=1:N
        tmp          = S(:,i) * (e(:,j))';
        if nnz(tmp)  == 0
            fprintf('Warning: A matrix all zero, skip.\n')
        else
            Acell{end+1} = 0.5*(tmp + tmp'); % make symmetric
        end
    end
end
fprintf('Done.\n')

% construct blk symmetric constraints
fprintf('Build blk sym constraints, progress ...')
for i = 1:k
    fprintf('%d/%d ',i,k);
    for j = i+1:k+1
        % each m by m block is symmetric
        for ii = 1:m-1
            for jj = ii+1:m
                row_shift = m*(i-1);
                col_shift = m*(j-1);
                tmp = sparse([ii+row_shift,jj+row_shift],[jj+col_shift,ii+col_shift],[-1,1],N,N);
                tmp = 0.5 * (tmp + tmp'); % make it symmetric
                Acell{end+1} = tmp;
            end
        end
    end
end
fprintf('Done.\n');

b         = zeros(length(Acell),1);
b(1)      = 1;

% At        = svec(blk,Acell,1);
% At        = At{1};

% AAt       = At'*At + (1e-14 * speye(size(At,2)));
% R         = chol(AAt);
% dR        = diag(R);
% smtol     = 1e-6;
% idxkp     = (abs(dR) > smtol);
% 
% fprintf('Keep %d / %d constraints.\n',full(sum(idxkp)),length(idxkp));

% At        = At(:,idxkp);
% Acell     = Acell(idxkp);
% b         = b(idxkp);

if norm(b) < 1e-6
    error('b is almost zero.');
end

% SDP.Acell = Acell;
SDP.m     = length(Acell);

At        = sparsesvec(blk,Acell);
Avect     = sparsevec(blk,Acell);

SDP.b     = b;
SDP.At    = {At};
SDP.Avect = Avect;
SDP.C     = {SDP.C};

% SDP.A     = At';

end



% % Block symmetrization of a matrix
% % Input: mk x mk matrix A
% % Output: symmetric matrix As such that
% %         each m x m block is also symmetric (there are k^2 blocks)
% function As = blksym(n,m,A)
% 
% As = mysym(A);
% if isempty(A); return; end
% 
% for i=1:n
%     for j=1:n
%         I = m*(i-1) + (1:m);
%         J = m*(j-1) + (1:m);
%         As(I,J,:) = mysym(As(I,J,:));
%     end
% end
% 
% end
% 
% 
% % symmetrize matrix
% function As = mysym(A)
% At = permute(A,[[2,1],3:ndims(A)]);
% As = .5*(A+At);
% end

