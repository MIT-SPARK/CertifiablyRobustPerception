function [Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,blk)
nblks           = size(blk,1);
if nblks == 1
    n               = blk{1,2};
    idxtril         = tril(true(n,n));
    Xopt            = zeros(n,n);
    Xopt(idxtril)   = res.sol.itr.barx;
    Xopt            = Xopt + Xopt';
    Xopt            = Xopt - 0.5*diag(diag(Xopt));
    Sopt            = zeros(n,n);
    Sopt(idxtril)   = res.sol.itr.bars;
    Sopt            = Sopt + Sopt';
    Sopt            = Sopt - 0.5*diag(diag(Sopt));
    
    Xopt            = {Xopt};
    Sopt            = {Sopt};
else
    Xopt  = cell(size(blk,1),1);
    Sopt  = cell(size(blk,1),1);
    cid   = 0;
    for blkid = 1:nblks
        n               = blk{blkid,2};
        ndelta          = triangle_number(n);
        idxtril         = tril(true(n,n));
        Xopti           = zeros(n,n);
        Xopti(idxtril)  = res.sol.itr.barx(cid+1:cid+ndelta);
        Xopti           = Xopti + Xopti';
        Xopti           = Xopti - 0.5*diag(diag(Xopti));
        Sopti           = zeros(n,n);
        Sopti(idxtril)  = res.sol.itr.bars(cid+1:cid+ndelta);
        Sopti           = Sopti + Sopti';
        Sopti           = Sopti - 0.5*diag(diag(Sopti));
        
        cid             = cid + ndelta;
        
        Xopt{blkid}     = Xopti;
        Sopt{blkid}     = Sopti;
    end
end
yopt            = res.sol.itr.y;
obj             = [res.sol.itr.pobjval;res.sol.itr.dobjval];
end