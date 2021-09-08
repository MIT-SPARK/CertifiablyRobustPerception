function info = get_performance_pcr(Xmat,y,Smat,SDP,problem,path)
%% Given (X, y, S) as the solution of the SDP relaxation
%% (1) Compute the standard relative KKT residuals of the SDP
%% (2) Round an approximate solution for point cloud registration
%% (3) Compute relative suboptimality gap from the rounded solution
%% (4) Compute rotation and translation estimation errors
%% Heng Yang
%% June 25, 2021

if iscell(path)
    for i = 1:length(path)
        addpath(genpath(path{i}))
    end
else
    addpath(genpath(path))
end

tBound      = problem.translationBound;
R_gt        = problem.R_gt;
t_gt        = problem.t_gt;

blk         = SDP.blk;
At          = SDP.At;
b           = SDP.b;
C           = svec(blk,SDP.C);

X           = svec(blk,Xmat);
S           = svec(blk,Smat);

Amap  = @(X) AXfun(blk,At,X);
ATmap = @(y) Atyfun(blk,At,y); 

Rp          = Fnorm(b - Amap(X))/(1+Fnorm(b));  
Rd          = Fnorm(ops(ops(C,'-',ATmap(y)),'-',S))/(1+Fnorm(C));

pobj        = blktrace(blk,C,X);
dobj        = b'*y;
gap         = abs(pobj - dobj)/(1+abs(pobj) + abs(dobj));

[R_est,t_est,theta_est] = round_pcr_v4(Xmat,tBound);
v           = lift_pcr_v3(R_est(:),t_est,theta_est);
f_est       = v'*SDP.C{1}*v;

R_err       = getAngularError(R_gt,R_est);
t_err       = getTranslationError(t_gt,t_est);

% S_mineig    = min(eig(SDP.C{1} - smat(blk(1,:),At{1}*y)));
S_mineig    = mineig( smat(blk, ops(C,'-',ATmap(y)) ) );
S_mineig_1  = min(S_mineig,0);
M           = SDP.M;
f_lb        = SDP.b'*y + M'*S_mineig_1;

eta         = abs(f_est - f_lb)/(1+abs(f_est)+abs(f_lb));


info.Rp     = Rp;
info.Rd     = Rd;
info.Rg     = gap;
info.Rs     = eta;
info.R_err  = R_err;
info.t_err  = t_err;
info.pobj   = pobj;
info.dobj   = dobj;
info.f_est  = f_est;
info.S_mineig = S_mineig;
info.f_lb   = f_lb;
info.theta_est = theta_est;
info.R_est  = R_est;
info.t_est  = t_est;

fprintf('\n=========== Point Cloud Registration Performance ===============\n')
fprintf('Rp: %3.2e, Rd: %3.2e, Rg: %3.2e, Rs: %3.2e.\n',Rp,Rd,gap,eta);
fprintf('f_est: %3.4e, f_lb: %3.4e, pobj: %3.4e, dobj: %3.4e.\n',f_est,f_lb,pobj,dobj);
fprintf('R_err: %3.4e, t_err: %3.4e.\n',R_err,t_err);
fprintf('================================================================\n')

if iscell(path)
    for i = 1:length(path)
        rmpath(genpath(path{i}))
    end
else
    rmpath(genpath(path))
end
end