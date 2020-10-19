function x = hatmap(u)
%HATMAP Summary of this function goes here
%   return skew symmetric matrix from a vector
x = [0, -u(3), u(2); 
    u(3), 0, -u(1); 
    -u(2), u(1), 0];
end

