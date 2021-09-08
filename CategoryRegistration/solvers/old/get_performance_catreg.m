function info = get_performance_catreg(Xopt,yopt,Sopt,SDP,problem,path)
%% Given (X,y,S) as output of the SDP solver, compute performance metrics
%% (1) Standard SDP relative KKT residuals 
%% (2) Rotation and translation errors
%% (3) Relative suboptimality of the rounded rank-one solution for the original nonconvex problem
%% Heng Yang, July 05, 2021
if iscell(path); for i = 1:length(path); addpath(genpath(path{i})); end; else; addpath(genpath(path)); end

tBound      = problem.translationBound;
cBound      = problem.cBound;
R_gt        = problem.R_gt;
t_gt        = problem.t_gt;
c_gt        = problem.c_gt;

blk         = SDP.blk;
At          = SDP.At;
b           = SDP.b;
C           = svec(blk,SDP.C);
Xoptmat     = Xopt;
Xopt        = svec(blk,Xopt);
Sopt        = svec(blk,Sopt);
Amap  = @(X) AXfun(blk,At,X);
ATmap = @(y) Atyfun(blk,At,y); 

%% standard SDP KKT residuals
Rp          = Fnorm(b - Amap(Xopt))/(1+Fnorm(b));  
Rd          = Fnorm(ops(ops(C,'-',ATmap(yopt)),'-',Sopt))/(1+Fnorm(C));
pobj        = blktrace(blk,C,Xopt);
dobj        = b'*yopt;
gap         = abs(pobj - dobj)/(1+abs(pobj) + abs(dobj));

%% round a rank-one solution
[R_est,t_est,c_est,theta_est] = round_catreg(Xoptmat,problem.N,problem.K,1);
v_est       = lift_catreg(R_est(:),t_est,c_est,theta_est,cBound,tBound);
f_est       = v_est{1}' * SDP.C{1} * v_est{1};

%% errors
R_err       = getAngularError(R_gt,R_est);
t_err       = getTranslationError(t_gt,t_est);
c_err       = getTranslationError(c_gt,c_est);

S_mineig    = min(eig(SDP.C{1} - smat(blk(1,:),At{1}*yopt)));
M           = SDP.M;
f_lb        = SDP.b'*yopt + M*S_mineig;
eta         = abs(f_est - f_lb)/(1+abs(f_est)+abs(f_lb));


info.Rp     = Rp;
info.Rd     = Rd;
info.Rg     = gap;
info.Rs     = eta;
info.R_err  = R_err;
info.t_err  = t_err;
info.c_err  = c_err;
info.pobj   = pobj;
info.dobj   = dobj;
info.f_est      = f_est;
info.S_mineig   = S_mineig;
info.f_lb       = f_lb;
info.theta_est  = theta_est;
info.R_est      = R_est;
info.t_est      = t_est;
info.c_est      = c_est;

fprintf('\n============ Performance of Category Registration ===========\n')
fprintf('Rp: %3.2e, Rd: %3.2e, Rg: %3.2e, Rs: %3.2e.\n',Rp,Rd,gap,eta);
fprintf('f_est: %3.4e, f_lb: %3.4e, pobj: %3.4e, dobj: %3.4e.\n',f_est,f_lb,pobj,dobj);
fprintf('R_err: %3.4e, t_err: %3.4e, c_err: %3.4e.\n',R_err,t_err,c_err);
fprintf('==============================================================\n')

if iscell(path); for i = 1:length(path); rmpath(genpath(path{i})); end; else; rmpath(genpath(path)); end

end