function v = lift_catreg(r,t,c,theta,cBound,tBound)
%% Lifting for category registration
%% Heng Yang, July 05, 2021
r           = r(:);
t           = t(:);
c           = c(:);
K           = length(c);
theta       = theta(:);
cBoundSq    = cBound^2;
tBoundSq    = tBound^2;

cr          = kron(c,r);
x           = [r;t];
v1          = [1;x;c;cr;theta;kron(theta,x);kron(theta,cr)];

if tBoundSq < t'*t
    v2      = [1;theta] * 0;
else
    v2      = [1;theta] * sqrt(tBoundSq - t'*t);
end

if cBoundSq < c'*c 
    tmp     = [1;r] * 0;
else
    tmp     = [1;r] * sqrt(cBoundSq - c'*c);
end

v           = {v1;v2;tmp};

for k = 1:K 
    if c(k) < 0
        tmp = [1;r] * 0;
    else
        tmp = [1;r] * sqrt(c(k));
    end
    v       = [v;{tmp}];
end

end