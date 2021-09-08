function [R_est,t_est,c_est,out] = outlier_free_category_registration(problem,path,varargin)
%% Solve outlier-free category registration problem as weighted least squares
%% Using a small SDP relaxation of fixed size
params = inputParser;
params.CaseSensitive = false;
params.addParameter('lambda',0.1, @(x) isscalar(x));
params.parse(varargin{:});
lambda         = params.Results.lambda;

N              = problem.N;
K              = problem.K;
scene          = problem.scene;
shapes         = problem.shapes;
noiseBoundSq   = problem.noiseBoundSq;
if ~isfield(problem,'weights'); problem.weights = ones(N,1); end
weights        = problem.weights;
% !!! This is important, the effective weights should be scaled by the
% noiseBoundSq !!!
weights        = weights / noiseBoundSq;

%% compute centers and demeaned versions of shapes and scenes
[yw,ybar]      = weighted_center(scene,weights); % ybar 3 by N
bbar           = zeros(3,N,K); 
bw             = zeros(3,K);
for k = 1:K 
    [center, bared] = weighted_center(shapes(:,:,k),weights);
    bw(:,k)         = center;
    bbar(:,:,k)     = bared;
end

% Ybar           = ybar(:); % Ybar 3N by 1
Bbar           = reshape(bbar,3*N,K); % Bbar 3N by K

%% c can be solved in closed-form from R
BtB            = Bbar' * Bbar + lambda * eye(K);
M              = (BtB) \ (Bbar');
H              = [Bbar * M - eye(3*N);...
                  sqrt(lambda)*M];
HtH            = H'*H;
P              = vecRt_r;
br             = kron(ybar',eye(3))*P; % b = br * r
Q              = br'*HtH*br; % 9 by 9
C              = [1,zeros(1,9);
                  zeros(9,1),Q];

%% Standard SDP data in SDPT3 format
addpath(genpath(path.stridepath))
C              = {sparse(C)};
Acell          = get_Amap;
blk{1,1}       = 's';
blk{1,2}       = 10;
At             = svec(blk,Acell);
b              = sparse(1,1,1,length(Acell),1);
rmpath(genpath(path.stridepath))
%% Solve using MOSEK
addpath(genpath(path.stridepath))
[At,b,c,K] = SDPT3data_SEDUMIdata(blk,At,C,b);
rmpath(genpath(path.stridepath))
prob       = convert_sedumi2mosek(At, b, c, K);
addpath(genpath(path.mosekpath))
[~,res]    = mosekopt('minimize info echo(0)',prob);
rmpath(genpath(path.mosekpath))
[Xopt,~,~,obj] = recover_mosek_sol_blk(res,blk);

Xopt           = Xopt{1,1};
[V,~]          = sorteig(Xopt);
v              = V(:,1); v = v/v(1);
R_est          = project2SO3(reshape(v(2:10),3,3));
r_est          = R_est(:);
f_est          = [1;r_est]' * C{1} * [1;r_est];
c_est          = M * (br * r_est);
t_est          = yw - R_est * (bw * c_est);
gap            = abs(f_est - obj(1))/(1+abs(f_est)+abs(obj(1)));

out.f_est      = f_est;
out.f_sdp      = obj(1);
out.gap        = gap;
end


function [center,A_decenter] = weighted_center(A,weights)
    % A size d by N, weights size N by 1
    center      = A * weights / sum(weights);
    A_decenter  = (A - center) .* sqrt(weights');
end

function Acell = get_Amap()
A0 = sparse(1,1,1,10,10);
A1 = sparse([1,2,3,4],...
            [1,2,3,4],...
            [1,-1,-1,-1],...
            10,10);
A2 = sparse([1,5,6,7],...
            [1,5,6,7],...
            [1,-1,-1,-1],...
            10,10);
A3 = sparse([1,8,9,10],...
            [1,8,9,10],...
            [1,-1,-1,-1],...
            10,10);
A4 = sparse([2,3,4],...
            [5,6,7],...
            [1,1,1],...
            10,10);
A4 = (A4 + A4')/2;
A5 = sparse([2,3,4],...
            [8,9,10],...
            [1,1,1],...
            10,10);
A5 = (A5+A5')/2;
A6 = sparse([5,6,7],...
            [8,9,10],...
            [1,1,1],...
            10,10);
A6 = (A6+A6')/2;
A7 = sparse([3,4,1],...
            [7,6,8],...
            [1,-1,-1],...
            10,10);
A7 = (A7 + A7')/2;
A8 = sparse([4,2,1],...
            [5,7,9],...
            [1,-1,-1],...
            10,10);
A8 = (A8+A8')/2;
A9 = sparse([2,1,3],...
            [6,10,5],...
            [1,-1,-1],...
            10,10);
A9 = (A9 + A9')/2;
A10= sparse([6,1,7],...
            [10,2,9],...
            [1,-1,-1],...
            10,10);
A10= (A10+A10')/2;
A11= sparse([7,5,1],...
            [8,10,3],...
            [1,-1,-1],...
            10,10);
A11= (A11 + A11')/2;
A12= sparse([5,1,6],...
            [9,4,8],...
            [1,-1,-1],...
            10,10);
A12= (A12 + A12')/2;
A13= sparse([4,3,1],...
            [9,10,5],...
            [1,-1,-1],...
            10,10);
A13= (A13 + A13')/2;
A14= sparse([2,1,4],...
            [10,6,8],...
            [1,-1,-1],...
            10,10);
A14= (A14 + A14')/2;
A15= sparse([3,2,1],...
            [8,9,7],...
            [1,-1,-1],...
            10,10);
A15= (A15+A15')/2;

Acell = {A0;A1;A2;A3;A4;A5;A6;A7;A8;A9;A10;A11;A12;A13;A14;A15};
end