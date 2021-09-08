function [Xmat,pobjround,info] = local_search_catreg_v2(Zmat,C,rrPar,opt,roundonly)
%% local search for category registration
%% used as a subroutine for STRIDE
%% Heng Yang
%% July 05, 2021

if nargin < 5
    roundonly = false;
end
if nargin < 4
    % default round the first two eigenvectors
    opt = [1,2]; 
end

tBound        = rrPar.translationBound;
cBound        = rrPar.cBound;
N             = rrPar.N;
K             = rrPar.K;
blk           = rrPar.blk;

if roundonly
    [R,t,c,theta] = round_catreg_v2(Zmat,N,K,[1]);
    xtld          = lift_catreg_v2(R(:),t,c,theta,cBound,tBound);
    Xmat          = rank_one_lift(xtld);
    pobjround     = blktrace(blk,Xmat,C);
else
    [R,t,c,theta] = round_catreg_v2(Zmat,N,K,opt);
    Zmatround     = {};
    pobjround     = zeros(length(opt),1);
    for i = 1:length(opt)
        Ri              = squeeze(R(:,:,i));
        ti              = squeeze(t(:,i));
        ci              = squeeze(c(:,i));
        thetai          = squeeze(theta(:,i));
        [Ztmp,pobjtmp]  = nlp_catreg_v2(C,rrPar,Ri,ti,ci,thetai);
        Zmatround{end+1}   = Ztmp;
        pobjround(i)       = pobjtmp;
    end
    pobjs            = pobjround;
    [pobjround,idx]  = min(pobjround);
    if pobjround == inf
        fprintf('        NLP fails to find a good solution, return a rounded solution only.\n');
        Ropt         = squeeze(R(:,:,1));
        ropt         = Ropt(:);
        topt         = t(:,1);
        copt         = c(:,1);
        thetaopt     = theta(:,1);
        xopt         = lift_catreg_v2(ropt,topt,copt,thetaopt,cBound,tBound);
        Xmat         = rank_one_lift(xopt);
        pobjround    = blktrace(blk,Xmat,C);
        nlpsuccess   = false;
    else
        Xmat = Zmatround{idx};
        nlpsuccess   = true;
    end
end

if nargout > 2
    info.minidx     = idx;
    info.nlpsuccess = nlpsuccess;
    info.pobjs      = pobjs;
    info.diffpobj   = pobjs(1) - pobjround;
end


end
    
    
    