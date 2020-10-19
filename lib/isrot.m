function isR = isrot(a,display)

if nargin < 2
    display = 0;
end

[r c] = size(a);
isR = false;
if r ~= 3
    if(display) disp('Matrix has not 3 rows'); end
elseif c ~= 3
    if(display) disp('Matrix has not 3 columns'); end
elseif norm ( a * a' - eye(3) ) > 1E-10
    if(display) disp('Matrix is not orthonormal, i.e. ||(R''R-I)|| > 1E-10'); end
elseif det (a) < 0
    if(display) disp('Matrix determinant is -1'); end
else isR = true;
end