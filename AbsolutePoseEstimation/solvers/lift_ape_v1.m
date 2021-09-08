function vcell = lift_ape_v1(r,t,theta,tBound,depthBound,FOV)
    r       = r(:);
    t       = t(:);
    x       = [r;t];
    theta   = theta(:);
    v       = [1;x;theta;kron(theta,x)];
    if tBound^2 < t'*t
        v1  = [1;theta] * 0;
    else
        v1  = [1;theta] * sqrt(tBound^2 - t'*t);
    end
    
    if t(end) < depthBound
        v2  = [1;theta] * 0; 
    else
        v2  = [1;theta] * sqrt(t(end)-depthBound);
    end
    
    scale = tan(deg2rad(FOV)/2);
    if scale^2 * t(3)^2 < t(1)^2 + t(2)^2
        v3  = [1;theta] * 0; 
    else
        v3  = [1;theta] * sqrt(scale^2 * t(3)^2 - t(1)^2 - t(2)^2);
    end
    
    vcell   = {v;v1;v2;v3};
end