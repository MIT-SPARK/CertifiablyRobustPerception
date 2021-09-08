function solution = outlier_free_absolute_pose_estimation(problem,varargin)
%% Solve absolute pose estimation using a weighted least squares formualation
%% The residual function is a point to line distance (point to bearing vector)
%% The weights of correspondences are given
%% The problem is solved using a small SDP relaxation, modelled by YALMIP and solved using MOSEK
%% Used as a non-minimal solver for graduated non-convexity
%% Heng Yang, June 28, 2021

params = inputParser;
params.CaseSensitive = false;
params.addParameter('relaxOrder',1, @(x) isscalar(x));
params.parse(varargin{:});

N               = problem.N;
meas3D          = problem.X;
meas2D          = problem.x;
tBound          = problem.translationBound;
dBound          = problem.depthBound;
noiseBoundSq    = problem.noiseBoundSq;
relaxOrder      = params.Results.relaxOrder;

if ~isfield(problem,'weights'); problem.weights = ones(N,1); end
weights         = problem.weights;

%% define decision variables
yalmip('clear')
time0       = tic;
r           = sdpvar(9,1);
R           = reshape(r,[3,3]);
r1          = R(:,1); 
r2          = R(:,2); 
r3          = R(:,3);
t           = sdpvar(3,1);

%% calculate cost function
f_cost      = 0;
for i = 1:N
    bearingi    = [meas2D(:,i);1];
    bearingi    = bearingi / norm(bearingi);
    pointi      = R * meas3D(:,i) + t;
    f_cost      = f_cost + weights(i) * ...
                  (pointi' * (eye(3) - bearingi*bearingi') * pointi)/noiseBoundSq;
end

%% define equality and inequality constraints
calr = [r1'*r1 - 1 == 0;...
        r2'*r2 - 1 == 0;...
        r3'*r3 - 1 == 0;...
        r1'*r2 == 0;...
        r2'*r3 == 0;...
        r3'*r1 == 0;...
        cross(r1,r2) - r3 == 0;...
        cross(r2,r3) - r1 == 0;...
        cross(r3,r1) - r2 == 0];
calt = [t'*t <= tBound^2;
        t(3) >= dBound];

calX = [calr;calt];

%% Generate dense first-order relaxation
[Fnew,f_cost_new,momentMat] = momentmodel(calX,f_cost,relaxOrder);
time_relax = toc(time0); % this is the time used to generate the relaxation

%% solve the resulting SDP
options = sdpsettings('solver', 'mosek', 'verbose', 0, 'debug',1);
optimize(Fnew,f_cost_new,options);
time_sdp = toc(time0) - time_relax; % this is the time used to solve the SDP

%% extract solutions and do certification
f_sdp = value(f_cost_new);
momentMat_val={};
for i = 1:length(momentMat)
    momentMat_val{i} = value(momentMat{i});
end

M1              = momentMat_val{2};
rank_moment(1)  = rank(M1,1e-3);

[V,D]       = eig(M1);
[~,idxsort] = sort(diag(D),'descend');
V           = V(:,idxsort);
v_monomials = V(:,1);
v_monomials = v_monomials/v_monomials(1);
r_est       = v_monomials(2:10);
R_est       = project2SO3( reshape(r_est,[3,3]) );
r_est       = R_est(:);
t_est       = v_monomials(11:13);


f_est           = replace(f_cost,[r;t],[r_est;t_est]);
absDualityGap   = abs(f_est - f_sdp);
relDualityGap   = absDualityGap / (1+abs(f_est)+abs(f_sdp));

residuals_est = zeros(N,1);
for i = 1:N
    bearingi    = [meas2D(:,i);1];
    bearingi    = bearingi / norm(bearingi);
    pointi      = R_est * meas3D(:,i) + t_est;
    residuals_est(i) = (pointi' * (eye(3) - bearingi*bearingi') * pointi)/noiseBoundSq;
end


solution.type               = 'LS';
solution.f_est              = f_est;
solution.f_sdp              = f_sdp;
solution.relaxOrder         = relaxOrder;
solution.absDualityGap      = absDualityGap;
solution.relDualityGap      = relDualityGap;
solution.R_est              = R_est;
solution.t_est              = t_est;
solution.Moments            = momentMat_val;
solution.rank_moment        = rank_moment;
solution.time_relax         = time_relax;
solution.time_sdp           = time_sdp;
solution.residuals          = residuals_est;
end