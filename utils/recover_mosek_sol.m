function [Xopt,yopt,Sopt,obj] = recover_mosek_sol(res,n)
idxtril         = tril(true(n,n));
Xopt            = zeros(n,n);
Xopt(idxtril)   = res.sol.itr.barx;
Xopt            = Xopt + Xopt';
Xopt            = Xopt - 0.5*diag(diag(Xopt));
yopt            = res.sol.itr.y;
Sopt            = zeros(n,n);
Sopt(idxtril)   = res.sol.itr.bars;
Sopt            = Sopt + Sopt';
Sopt            = Sopt - 0.5*diag(diag(Sopt));
obj             = [res.sol.itr.pobjval;res.sol.itr.dobjval];
end