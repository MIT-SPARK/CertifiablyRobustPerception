function info = get_performance_sra(Xopt,yopt,Sopt,SDP,problem,pgdpath)
%% given (X,y,S) as output from the SDP solver
%% compute standard relative KKT residuals of the SDP
%% compute relative suboptimality of the original nonconvex promblem
%% compute rotation estimation error
%% Heng Yang, June 29, 2021

addpath(genpath(pgdpath))
R_gt        = problem.R_gt;

Xopt        = Xopt{1};
Sopt        = Sopt{1};
At          = SDP.At{1};

blk         = SDP.blk;
Rp          = norm(At'*svec(blk,Xopt) - SDP.b)/(1+norm(SDP.b));
Rd          = Fnorm(SDP.C{1} - smat(blk,At*yopt) - Sopt)/(1+Fnorm(SDP.C));
pobj        = trace(SDP.C{1}*Xopt);
dobj        = SDP.b'*yopt;
gap         = abs(pobj - dobj)/(1+abs(pobj) + abs(dobj));

[R_est,theta_est] = round_sra({Xopt});
v_est       = lift_sra(R_est(:),theta_est);
f_est       = v_est{1}' * SDP.C{1} * v_est{1};
R_err       = getAngularError(R_gt,R_est);
S_mineig    = min(eig(SDP.C{1} - smat(blk,At*yopt))); % single block
M           = SDP.M;
f_lb        = SDP.b'*yopt + M * min(S_mineig,0);
eta         = abs(f_est - f_lb)/(1+abs(f_est)+abs(f_lb));


info.Rp     = Rp;
info.Rd     = Rd;
info.Rg     = gap;
info.Rs     = eta;
info.R_err  = R_err;
info.pobj   = pobj;
info.dobj   = dobj;
info.theta_est = theta_est;
info.f_est  = f_est;
info.f_lb   = f_lb;


fprintf('\n========== Performance of Single Rotation Averaging ==========\n')
fprintf('Rp: %3.2e, Rd: %3.2e, Rg: %3.2e, Rs: %3.2e.\n',Rp,Rd,gap,eta);
fprintf('f_est: %3.4e, f_lb: %3.4e, pobj: %3.4e, dobj: %3.4e.\n',f_est,f_lb,pobj,dobj);
fprintf('R_err: %3.4e.\n',R_err);
fprintf('===============================================================\n')

rmpath(genpath(pgdpath))
end