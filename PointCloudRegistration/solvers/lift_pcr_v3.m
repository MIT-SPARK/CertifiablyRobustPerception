function v = lift_pcr_v3(r,t,theta)
    r       = r(:);
    t       = t(:);
    x       = [r;t];
    theta   = theta(:);
    v       = [1;x;theta;kron(theta,x)];
end