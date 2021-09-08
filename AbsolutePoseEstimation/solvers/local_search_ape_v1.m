function [Xmat,pobjround,info] = local_search_ape_v1(Zmat,C,rrPar,opt,roundonly)
%% local search for absolute pose estimation
%% used as a subroutine for STRIDE
%% Heng Yang, June 28, 2021

if nargin < 5
    roundonly = false;
end
if nargin < 4
    % default round the first two eigenvectors
    opt = [1,2]; 
end

tBound        = rrPar.translationBound;
dBound        = rrPar.depthBound;
blk           = rrPar.blk;
FOV           = rrPar.FOV;

if roundonly
    [R,t,theta] = round_ape_v1(Zmat,tBound,dBound,[1]);
    xtld        = lift_ape_v1(R(:),t,theta,tBound,dBound,FOV);
    Xmat        = {xtld{1}*xtld{1}';xtld{2}*xtld{2}';xtld{3}*xtld{3}';xtld{4}*xtld{4}'};
    pobjround   = blktrace(blk,Xmat,C);
else
    [R,t,theta] = round_ape_v1(Zmat,tBound,dBound,opt);
    Zmatround   = {};
    pobjround   = zeros(length(opt),1);
    for i = 1:length(opt)
        Ri              = squeeze(R(:,:,i));
        ti              = squeeze(t(:,i));
        thetai          = squeeze(theta(:,i));
        [out,pobjtmp]   = nlp_ape(C,rrPar,Ri,ti,thetai);
        Zmatround{end+1}   = out.X;
        pobjround(i)       = pobjtmp;
    end
    pobjs            = pobjround;
    [pobjround,idx]  = min(pobjround);
    if pobjround == inf
        fprintf('        NLP fails to find a good solution, return a rounded solution only.\n');
        Ropt         = squeeze(R(:,:,1));
        ropt         = Ropt(:);
        topt         = t(:,1);
        thetaopt     = theta(:,1);
        xopt         = lift_ape_v1(ropt,topt,thetaopt,tBound,dBound,FOV);
        Xmat         = {xopt{1} * xopt{1}';xopt{2} * xopt{2}';xopt{3}*xopt{3}';xopt{4}*xopt{4}'};
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
        
            