function vcell = lift_pcr_v4(r,t,theta,tBound)
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
    vcell   = {v;v1};
end