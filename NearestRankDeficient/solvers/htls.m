% interpmethod is either 'linear' or 'spline'
function [z,Hu] = htls(k,y,interpmethod)
% "Algorithm for time-domain NMR data fitting based on total least squares"
% Vanhuffel, S., Chen, H., Decanniere, C. and Vanhecke, P., (1994)
if nargin < 3; interpmethod = 'linear'; end

y = complete_signal(y,interpmethod);
n = length(y);
y = reshape(y,[n,1]);
l = ceil(n/2);
Hy = hankel(y(1:l),y(l:n));
[U,~,~] = svd(Hy,'econ');
Ur = U(:,1:k);
Ur1 = Ur(2:end,:);
Ur2 = Ur(1:end-1,:);
[~,~,V] = svd([Ur2 Ur1]);
I = 1:k;
V12 = V(I,k+I);
V22 = V(k+I,k+I);
Q = -V12/V22;
e = eig(Q);
A = (e.^(0:n-1)).';
c = A \ y;
u = A*c;
Hu = hankel(u(1:k+1),u(k+1:n));
[U,~,~] = svd(Hu,'econ');
z = U(:,end);

function y = complete_signal(y,interpmethod)
I = isnan(y);
if ~any(I); return; end
t = 1:length(y);
y(I) = interp1(t(~I),y(~I),t(I),interpmethod);