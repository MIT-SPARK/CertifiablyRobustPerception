function info = get_performance_quasar(Xopt,yopt,Sopt,SDP,R_gt)

%% Compute standard KKT residuals
blk         = SDP.blk;
At          = SDP.At;
b           = SDP.b;
C           = svec(blk,SDP.C);

X           = svec(blk,Xopt);
S           = svec(blk,Sopt);

Amap        = @(X) AXfun(blk,At,X);
ATmap       = @(y) Atyfun(blk,At,y); 

Rp          = Fnorm(b - Amap(X))/(1+Fnorm(b));  
Rd          = Fnorm(ops(ops(C,'-',ATmap(yopt)),'-',S))/(1+Fnorm(C));

pobj        = blktrace(blk,C,X);
dobj        = b'*yopt;
gap         = abs(pobj - dobj)/(1+abs(pobj) + abs(dobj));

%% Round a solution from Xopt
n     = size(Xopt{1},1);
N     = n/4 - 1; % number of measurements
[V,D] = eig(Xopt{1});
[~,I] = sort(diag(D),'descend');
V     = V(:,I);
x     = V(:,1);
q     = x(blkIndices(1,4));
q     = q / norm(q);
theta = zeros(N,1);
for i = 1:N
    inprod = q' * x(blkIndices(i+1,4));
    if inprod > 0
        theta(i) = 1;
    else
        theta(i) = -1;
    end
end

%% Compute suboptimality of the rounded solution
x_est       = kron([1;theta],q);
f_est       = x_est' * SDP.C{1} * x_est;
S_mineig    = mineig( smat(blk, ops(C,'-',ATmap(yopt)) ) );
M           = SDP.M;
f_lb        = dobj + M * min(S_mineig,0);
eta         = abs(f_est - f_lb)/(1+abs(f_est)+abs(f_lb));

R_err       = getAngularError(R_gt,quat2rot(q));

info.Rp     = Rp;
info.Rd     = Rd;
info.Rg     = gap;
info.Rs     = eta;
info.R_err  = R_err;
info.pobj   = pobj;
info.dobj   = dobj;
info.R      = quat2rot(q);
info.q      = q;
info.theta  = theta;

fprintf('\n=================== QUASAR Performance =======================\n')
fprintf('Rp: %3.2e, Rd: %3.2e, Rg: %3.2e, Rs: %3.2e.\n',Rp,Rd,gap,eta);
fprintf('f_est: %3.4e, f_lb: %3.4e, pobj: %3.4e, dobj: %3.4e.\n',f_est,f_lb,pobj,dobj);
fprintf('R_err: %3.4e.\n',R_err);
fprintf('================================================================\n')
end