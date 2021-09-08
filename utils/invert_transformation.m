function [R1,t1] = invert_transformation(R,t)
R1      = R';
t1      = - R'*t;
end