function v = lift_catreg_v2(r,t,c,theta,cBound,tBound)
%% Lifting for category registration
%% Heng Yang, July 05, 2021
r           = r(:);
t           = t(:);
c           = c(:);
K           = length(c);
theta       = theta(:);
cBoundSq    = cBound^2;
tBoundSq    = tBound^2;

x           = [r;t;c];
v1          = [1;x;theta;kron(theta,x)];

if tBoundSq < t'*t
    v2      = [1;theta] * 0;
else
    v2      = [1;theta] * sqrt(tBoundSq - t'*t);
end

if cBoundSq < c'*c 
    tmp     = [1;theta] * 0;
else
    tmp     = [1;theta] * sqrt(cBoundSq - c'*c);
end

v           = {v1;v2;tmp};

for k = 1:K 
    if c(k) < 0
        tmp = [1;theta] * 0;
    else
        tmp = [1;theta] * sqrt(c(k));
    end
    v       = [v;{tmp}];
end

end