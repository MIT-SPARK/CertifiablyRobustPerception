function SDP = chordal_relax_absolute_pose_estimation(problem)
%% Apply a chordal sparse second-order relaxation to absolute pose estimation
%% Full perspective camera, R and t relaxed together
%% Depending on multivariate polynomial package in SPOT
%% Translation bounded modeled as an extra PSD constraint
%% t(3) >= depthBound modelled as an extra PSD constraint
%% tan(FOV/2)^2 * t(3)^2 >= t(1)^2 + t(2)^2 in FOV cone, extra PSD constraint (July 01)
%% Generate multiple smaller PSD blocks, can be solved efficiently by IPMs
%% But typically not tight compared to the one-block SDP relaxation
%% Heng Yang, July 01, 2021

fprintf('\n===================================================================')
fprintf('\nApplying chordal SDP relaxation to absolute pose estimation')
fprintf('\n===================================================================\n')
t0              = tic;

N               = problem.N;
meas3D          = problem.X;
meas2D          = problem.x;
noiseBoundSq    = problem.noiseBoundSq;
tBound          = problem.translationBound;
tBoundSq        = tBound^2; % t'*t <= tBoundSq
dBound          = problem.depthBound;
barc2           = 1.0;
FOV             = deg2rad(problem.FOV);

%% define POP variables
nrPrimalVars    = 9+3+N;
p               = msspoly('p',nrPrimalVars);
r               = p(1:9);
R               = reshape(r,3,3); 
col1            = r(1:3);
col2            = r(4:6);
col3            = r(7:9);
row1            = R(1,:)';
row2            = R(2,:)';
row3            = R(3,:)';
t               = p(10:12);
x               = [r;t];
theta           = p(13:nrPrimalVars);

%% compute the cost function
residuals = [];
for i = 1:N 
    bearingi    = [meas2D(:,i);1];
    bearingi    = bearingi / norm(bearingi);
    pointi      = R * meas3D(:,i) + t;
    res         = pointi' * (eye(3) - bearingi*bearingi') * pointi; % point to line distance
    residuals   = [residuals; res / noiseBoundSq];
end
f_cost = [];
for i = 1:N 
    f_cost = [f_cost; (1+theta(i))/2 * residuals(i) + (1-theta(i))/2 * barc2];
end

%% define equality constraints
% h_x  = [1.0-col1'*col1;...
%        1.0-col2'*col2;...
%        1.0-col3'*col3;... % columns unit length
%        col1'*col2;...
%        col2'*col3;...
%        col3'*col1;... % columns orthogonal
%        cross(col1,col2) - col3;...
%        cross(col2,col3) - col1;...
%        cross(col3,col1) - col2]; % columns righthandedness
   
h_x  = [1.0-col1'*col1;...
       1.0-col2'*col2;...
       1.0-col3'*col3;... % columns unit length
       1.0-row1'*row1;...
       1.0-row2'*row2;...
       1.0-row3'*row3;... % rows unit length
       col1'*col2;...
       col2'*col3;...
       col3'*col1;... % columns orthogonal
       row1'*row2;...
       row2'*row3;...
       row3'*row1;... % rows orthogonal
       cross(row1,row2) - row3;...
       cross(row2,row3) - row1;...
       cross(row3,row1) - row2;... % rows righthandedness
       cross(col1,col2) - col3;...
       cross(col2,col3) - col1;...
       cross(col3,col1) - col2]; % columns righthandedness

h_theta = [];
for i = 1:N 
    h_theta =[h_theta; 1-theta(i)^2];
end

g_x = [tBoundSq - t'*t; ... % bounded translation
       t(3)-dBound; ... % in front of camera
       (tan(FOV/2))^2 * t(3)^2 - (t(1)^2 + t(2)^2)]; % inside FOV cone
       

%% Formulate the chordal sparse second-order relaxation
%% the 0-th block [1;x] * [1;x]'
basis0          = [1;x];
n0              = length(basis0);
basis_x0        = 1;

pop0            = [mykron(basis_x0,h_x);...
                   mykron(basis0,basis0)];
[~,degmat,coef_all] = decomp(pop0);
coef_all            = coef_all';
dim_loc0        = length(basis_x0) * length(h_x);
n0delta         = triangle_number(n0);
nterms          = size(degmat,1);   
m_mom0          = n0delta - nterms;

assert(m_mom0==0,'The zero-th blk should have 0 moment constraints.')

coef_mom    = coef_all(:,dim_loc0+1:end);
coef_mom    = coef_mom';
B           = {};
B_normalize = {};

for i = 1:nterms
    [row,~,~]   = find(coef_mom(:,i));
    SDP_coli    = floor((row-1)./n0) + 1;
    SDP_rowi    = mod(row-1,n0) + 1;
    nnz         = length(SDP_rowi);
    
    Bi          = sparse(SDP_rowi,SDP_coli,ones(nnz,1),n0,n0);
    B{end+1}    = Bi;
    B_normalize{end+1} = Bi/nnz;
end

coef_loc0       = coef_all(:,1:dim_loc0);
A0_local        = {};

for i = 1:dim_loc0
    [rowi,~,vi] = find(coef_loc0(:,i));
    Ai      = sparse(n0,n0);
    for j   = 1:length(rowi)
        Ai  = Ai + vi(j) * B_normalize{rowi(j)};
    end
    A0_local = [A0_local;{Ai}];
end
A0_0    = sparse([1],[1],[1],n0,n0);
% The first block satisfies A0(X0) = b0;
A0      = [{A0_0};A0_local];
b0      = sparse(1,1,1,length(A0),1);

%% the 1-N blocks [1;x;theta(i);theta(i)*x] * [1;x;theta(i);theta(i)*x]'
%% Since there are three inequality constraints, it will generate 4*N blocks
Aall            = {};
A1all           = {};
A2all           = {};
A3all           = {};
ball            = [];
Call            = {};
A0append        = {};
for blkidx = 1:N
    basis       = [1;x;theta(blkidx);theta(blkidx)*x];
    n           = length(basis);
    basis_x     = monomials(theta(blkidx),1:2); % perhaps only need 1 since theta^2 = 1
    basis_theta = monomials(x,0:2);
    basis_g     = [1;theta(blkidx)];
    
    out         = gen_chordal_subblk_ape(...
        basis,basis_x,h_x,basis_theta,h_theta(blkidx),g_x,basis_g,f_cost(blkidx));
    
    Acell       = out.A;
    A1          = out.A1;
    A2          = out.A2;
    A3          = out.A3;
    b           = out.b;
    C           = out.C;
    
    % add constraint that the top-left [1;x]*[1;x]' block is the same as
    % the 0-th block
    A0blk = {};
    for i = 1:n0
        for j = i:n0
            if i == j
                A0i  = sparse(i,j,-1,n0,n0);
                Ai   = sparse(i,j,1,n,n);
            else
                A0i  = sparse([i,j],[j,i],[-0.5,-0.5],n0,n0);
                Ai   = sparse([i,j],[j,i],[0.5,0.5],n,n);
            end
            A0blk    = [A0blk;{A0i}];
            Acell    = [Acell;{Ai}];
        end
    end
    
    ball             = [ball;b;sparse(n0delta,1)];
    A0append{end+1}  = A0blk;
    Aall{end+1}      = Acell;
    Call             = [Call;C];
    A1all{end+1}     = A1;
    A2all{end+1}     = A2;
    A3all{end+1}     = A3;
end

%% Convert to standard SDPT3 format
b           = [b0;ball];
blk         = cell((length(g_x)+1*N)+1,2);
blk{1,1}    = 's'; blk{1,2} = n0;
n1          = out.blk{2,2};
n1delta     = triangle_number(n1);
for i = 1:N
    blk{4*i-2,1} = 's';
    blk{4*i-2,2} = n;
    
    blk{4*i-1,1} = 's';
    blk{4*i-1,2} = n1;
    
    blk{4*i,1} = 's';
    blk{4*i,2} = n1;
    
    blk{4*i+1,1} = 's';
    blk{4*i+1,2} = n1;
end

A0t     = sparsesvec(blk(1,:),A0);
for i = 1:N
    A0t = [A0t,...
           sparse(n0delta,out.m),...
           sparsesvec(blk(1,:),A0append{i})];
end
ndelta  = triangle_number(n);
At    = {A0t};

for i = 1:N
    Ait = [sparse(ndelta,length(b0)),...
           sparse(ndelta,(i-1)*length(Aall{i})),...
           sparsesvec(blk(4*i-2,:),Aall{i}),... % the moment matrix block
           sparse(ndelta,(N-i)*length(Aall{i}))];
       
    A1it = [sparse(n1delta,length(b0)),...
            sparse(n1delta,(i-1)*length(Aall{i})),...
            sparse(n1delta,out.m_mom+out.m_loc),... % the moment constraint
            sparsesvec(blk(4*i-1,:),A1all{i}),... % first sub psd blk
            sparse(n1delta,n1delta),... % second sub psd blk
            sparse(n1delta,n1delta),... % third sub psd blk
            sparse(n1delta,n0delta),...
            sparse(n1delta,(N-i)*length(Aall{i}))];
    
    A2it = [sparse(n1delta,length(b0)),...
            sparse(n1delta,(i-1)*length(Aall{i})),...
            sparse(n1delta,out.m_mom+out.m_loc),... % the moment constraint
            sparse(n1delta,n1delta),... % first sub psd blk
            sparsesvec(blk(4*i,:),A2all{i}),... % second sub psd blk
            sparse(n1delta,n1delta),... % third sub psd blk
            sparse(n1delta,n0delta),...
            sparse(n1delta,(N-i)*length(Aall{i}))];
    
    A3it = [sparse(n1delta,length(b0)),...
            sparse(n1delta,(i-1)*length(Aall{i})),...
            sparse(n1delta,out.m_mom+out.m_loc),... % the moment constraint
            sparse(n1delta,n1delta),... % first sub psd blk
            sparse(n1delta,n1delta),... % second sub psd blk
            sparsesvec(blk(4*i+1,:),A3all{i}),... % third sub psd blk
            sparse(n1delta,n0delta),...
            sparse(n1delta,(N-i)*length(Aall{i}))];
        
    At  = [At;{Ait};{A1it};{A2it};{A3it}];
end


SDP.blk = blk;
SDP.At  = At;
SDP.n   = n;
SDP.m   = length(b);
SDP.C   = [{sparse(n0,n0)};Call];
SDP.b   = b;

tf    = toc(t0);
fprintf('\nDone in %g seconds.\n',tf);
fprintf('===================================================================\n')

%% Convert to sedumi format
fprintf('\nConvert to sedumi format.\n');
t0    = tic;

sK.s  = [n0];
for i = 1:N
    sK.s        = [sK.s,n,n1,n1,n1];
end

A0t     = sparsevec(blk(1,:),A0);
n0sq    = n0^2;
for i = 1:N
    A0t = [A0t,...
           sparse(n0sq,out.m),...
           sparsevec(blk(1,:),A0append{i})];
end

nsq     = n^2;
n1sq    = n1^2;
n1delta = triangle_number(n1);

At    = {A0t};

for i = 1:N
    Ait = [sparse(nsq,length(b0)),...
           sparse(nsq,(i-1)*length(Aall{i})),...
           sparsevec(blk(4*i-2,:),Aall{i}),... % the moment matrix block
           sparse(nsq,(N-i)*length(Aall{i}))];
       
    A1it = [sparse(n1sq,length(b0)),...
            sparse(n1sq,(i-1)*length(Aall{i})),...
            sparse(n1sq,out.m_mom+out.m_loc),... % the moment constraint
            sparsevec(blk(4*i-1,:),A1all{i}),... % first sub psd blk
            sparse(n1sq,n1delta),... % second sub psd blk
            sparse(n1sq,n1delta),... % third sub psd blk
            sparse(n1sq,n0delta),...
            sparse(n1sq,(N-i)*length(Aall{i}))];
    
    A2it = [sparse(n1sq,length(b0)),...
            sparse(n1sq,(i-1)*length(Aall{i})),...
            sparse(n1sq,out.m_mom+out.m_loc),... % the moment constraint
            sparse(n1sq,n1delta),... % first sub psd blk
            sparsevec(blk(4*i,:),A2all{i}),... % second sub psd blk
            sparse(n1sq,n1delta),... % third sub psd blk
            sparse(n1sq,n0delta),...
            sparse(n1sq,(N-i)*length(Aall{i}))];
    
    A3it = [sparse(n1sq,length(b0)),...
            sparse(n1sq,(i-1)*length(Aall{i})),...
            sparse(n1sq,out.m_mom+out.m_loc),... % the moment constraint
            sparse(n1sq,n1delta),... % first sub psd blk
            sparse(n1sq,n1delta),... % second sub psd blk
            sparsevec(blk(4*i+1,:),A3all{i}),... % third sub psd blk
            sparse(n1sq,n0delta),...
            sparse(n1sq,(N-i)*length(Aall{i}))];
        
    At  = [At;{Ait};{A1it};{A2it};{A3it}];
end

sdata.K     = sK;
sdata.At    = cat(1,At{:});
sdata.b     = b;

sc          = [];
for i = 1:length(SDP.C)
    sc      = [sc;sparsevec(blk(i,:),SDP.C(i))];
end
sdata.c     = sc;

SDP.sedumi  = sdata;


tf    = toc(t0);
fprintf('\nDone in %g seconds.\n',tf);
fprintf('===================================================================\n')
end