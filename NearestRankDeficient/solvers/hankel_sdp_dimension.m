function [n,m] = hankel_sdp_dimension(N1,N2)
N       = N1 + N2 - 1;
n       = (N+1)*N1;
m       = 1 + n*N2 + (N1-1)*N1/2 * (N+1)*N/2;
end