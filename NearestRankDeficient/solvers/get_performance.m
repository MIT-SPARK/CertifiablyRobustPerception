function info = get_performance(Xopt,yopt,Sopt,SDP,n1,k,pgdpath)
addpath(genpath(pgdpath))

pobj        = trace(SDP.C * Xopt);
dobj        = SDP.b' * yopt;

X           = svec(SDP.blk,Xopt);
S           = svec(SDP.blk,Sopt);
C           = svec(SDP.blk,SDP.C);
Rp          = norm(SDP.At' * X - SDP.b)/(1+norm(SDP.b));
Rd          = Fnorm(C - SDP.At*yopt - S) / (1+Fnorm(C));
Rg          = (pobj - dobj)/(1+abs(pobj)+abs(dobj));

mineig      = min(eig(smat(SDP.blk,C-SDP.At*yopt)));

f_lb        = dobj + SDP.M * mineig;

[~,~,f_est] = recover_solution(Xopt,SDP.C,n1,k);

Rs          = abs(f_est - f_lb) / (1+abs(f_est)+abs(f_lb));

info.pobj   = pobj;
info.dobj   = dobj;
info.Rp     = Rp;
info.Rd     = Rd;
info.Rg     = Rg;
info.Rs     = Rs;

rmpath(genpath(pgdpath))
end