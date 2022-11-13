function problem = QUASAR_Problem(v1,v2,barc2,varargin)
% return the cost function and constraint matrices of QUASAR

fprintf('Relax Wahba problem to standard linear SDP...')

params = inputParser;
params.CaseSensitive = false;
params.addParameter('computeConstraintMatrix', 'true', @(x) islogical(x));
params.addParameter('TraceAll','true',@(x) islogical(x));
params.parse(varargin{:});

computeConstraintMatrix = params.Results.computeConstraintMatrix;
TraceAll                = params.Results.TraceAll;

problem.computeConstraintMatrix = computeConstraintMatrix;

N = size(v1,2); % nr vector
% X = [q\tran, q_1\tran, q_2\tran, ..., q_N\tran] has size 1 by Npm, where Npm = a+b+4*N
Npm=4+4*N; 

% coefficient matrix that maps vec(qq\tran) to vec(R)
P=[1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1;
   0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0;
   0, 0, 1, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, -1, 0, 0;
   0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0;
   -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1;
   0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0;
   0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0;
   0, 0, 0, -1, 0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0, 0;
   -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
P=sparse(P);

% build the cost matrix Q_cost
Q_1=zeros(Npm,Npm);
for k=1:N
    idx = 4+blkIndices(k,4);
    P_k = reshape(P'*reshape(v2(:,k)*v1(:,k)',[9,1]),[4,4]);
    ck = 0.5 * ( v1(:,k)'*v1(:,k)+v2(:,k)'*v2(:,k) - barc2 );
    Q_1((1:4),idx) = Q_1((1:4),idx)-0.5*P_k+ck/2*eye(4);
    Q_1(idx,(1:4)) = Q_1(idx,(1:4))-0.5*P_k+ck/2*eye(4);
end

Q_2=zeros(Npm,Npm);
for k=1:N
%     idx = 4+blkIndices(k,4);
    idx = blkIndices(1,4);
    P_k = reshape(P'*reshape(v2(:,k)*v1(:,k)',[9,1]),[4,4]);
    ck = 0.5 * ( v1(:,k)'*v1(:,k)+v2(:,k)'*v2(:,k) + barc2 );
    Q_2(idx,idx) = Q_2(idx,idx) - P_k + ck*eye(4);
end
Q_cost=Q_1+Q_2;
Q_cost=sparse(Q_cost);

Q_cost=Q_cost/barc2;

problem.P = P;
problem.C = {Q_cost};
problem.n = size(Q_cost,1);
n = problem.n;

if computeConstraintMatrix
    % compute all the constraint matrices
    nr_constraints = 3*N*(N+1) + 10*N + 1;
    % nr_licq = 4*N+1;
    nr_licq = 10*N+1;
    nr_redundant = nr_constraints - nr_licq;
    b = sparse([1],[1],[1],nr_constraints,1);
    A = {};
    % trace([Z]_qq) = 1
    if TraceAll
        A_trace  = speye(n,n);
        b(1)     = N + 1;
    else
        A_trace = sparse([1,2,3,4],[1,2,3,4],[1,1,1,1],n,n);
        b(1) = 1;
    end
    A{end+1} = A_trace;
    
    % [Z]_qiqi - [Z]_qq = 0
    for k = 1:N
        for i = 1:4
            for j = i:4
                row_shift = 4*k;
                col_shift = 4*k;
                tmp = sparse([i,i+row_shift],[j,j+col_shift],[-1,1],n,n);
                tmp = tmp + tmp'; % make it symmetric
                tmp = tmp / norm(tmp,'fro'); % normalize to frobenius norm=1
                A{end+1} = tmp;
            end
        end
    end
  
    % [Z]_qiqj symmetric
    fprintf('\nBuild symmetric constraints, progress ... ')
    for k1 = 1:N
        fprintf('%d/%d ',k1,N);
        for k2 = k1+1:N+1
            for i = 1:3
                for j = i+1:4
                    row_shift = 4*(k1-1);
                    col_shift = 4*(k2-1);
                    tmp = sparse([i+row_shift,j+row_shift],[j+col_shift,i+col_shift],[-1,1],n,n);
                    tmp = tmp + tmp'; % make it symmetric
                    tmp = tmp / 2; % normalize to frobenius norm=1
                    A{end+1} = tmp;
                end
            end
        end
    end
    fprintf('Done.\n')
    %{
%*************************************************************
%% added by Kim-Chuan Toh
%*************************************************************  
    cnt = 1;
    len = 3*N*(N+1); rr = 1/sqrt(2); 
    row = zeros(len,1); col=zeros(len,1); val=zeros(len,1);
    fprintf('Build symmetric constraints, progress ... ')
    for k1 = 1:N
        fprintf('%d/%d ',k1,N);
        for k2 = k1+1:N+1
            for i = 1:3
                for j = i+1:4
                    row_shift = 4*(k1-1);
                    col_shift = 4*(k2-1);     
                    ii = i+row_shift; ii2 = j+row_shift;
                    jj = j+col_shift; jj2 = i+col_shift; 
                    idx = 2*cnt-1:2*cnt; 
                    row(idx) = [ii+(jj-1)*jj/2; ii2+(jj2-1)*jj2/2];
                    col(idx) = [cnt; cnt];
                    val(idx) = [-rr; rr];  
                    cnt = cnt+1;
                end
            end
        end
    end
    problem.Bt = spconvert([row,col,val;n*(n+1)/2,len,0]); 
    fprintf('Done.\n')
%************************************************************* 
%************************************************************* 
    %}
    assert(length(A) == nr_constraints, 'Incorrect number of constraint matrices in A.');
    problem.Acell = A;
    problem.b = b;
    problem.m_localize = nr_licq;
    problem.m_moment = nr_redundant;
    problem.m = length(A);
    problem.n = size(A{1},1);
    ndelta    = triangle_number(n);
    l         = ndelta - problem.m;
    blk{1,1} = 's';
    blk{1,2} = n;
    problem.blk = blk;
    
    At       = sparsesvec(blk,problem.Acell);
    problem.At = {At};
    problem.M  = N + 1;
    
    fprintf('Done.\n')
    fprintf('Linear SDP: n = %d, m = %d, m_localize = %d, m_moment = %d, l = %d.\n',...
        problem.n,problem.m,problem.m_localize,problem.m_moment,l);
    
%     problem      = rmfield(problem,{'P','computeConstraintMatrix','Acell','m_localize','m_moment'});
end








