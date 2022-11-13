function info = get_performance_q4s(Xopt,yopt,Sopt,SDP,POP,path)

if iscell(path)
    for i = 1:length(path)
        addpath(genpath(path{i}))
    end
else
    addpath(genpath(path))
end

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
Xmom        = Xoptmat{1};
[V,~]       = sorteig(Xmom);
v           = V(:,1)/V(1,1);
x_est       = v(2:1+POP.d);
x_est       = x_est / norm(x_est);
f_est       = POP.f.coefficient * ( prod(x_est.^(POP.f.degmat),1) )';

S_mineig    = mineig( smat(blk, ops(C,'-',ATmap(yopt)) ) );
S_mineig_1  = min(S_mineig,0);
M           = SDP.M;
f_lb        = SDP.b'*yopt + M'*S_mineig_1;
eta         = abs(f_est - f_lb)/(1+abs(f_est)+abs(f_lb));


info.Rp     = Rp;
info.Rd     = Rd;
info.Rg     = gap;
info.Rs     = eta;
info.pobj   = pobj;
info.dobj   = dobj;
info.f_est      = f_est;
info.S_mineig   = S_mineig;
info.f_lb       = f_lb;
info.x_est      = x_est;

fprintf('\n========== Performance of Quartic Opt Over Sphere ==========\n')
fprintf('Rp: %3.2e, Rd: %3.2e, Rg: %3.2e, Rs: %3.2e.\n',Rp,Rd,gap,eta);
fprintf('f_est: %3.4e, f_lb: %3.4e, pobj: %3.4e, dobj: %3.4e.\n',f_est,f_lb,pobj,dobj);
fprintf('==============================================================\n')

if iscell(path)
    for i = 1:length(path)
        rmpath(genpath(path{i}))
    end
else
    rmpath(genpath(path))
end
end