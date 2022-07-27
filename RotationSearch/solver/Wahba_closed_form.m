function R = Wahba_closed_form(v1, v2)
% closed form solution of the Wahba problem
if size(v1,1) ~= 3 || size(v2,1) ~= 3
    error('Wahab_closed_form works on N 3D vectors (3xN matrices)')
end

B = v2 * v1';
[U,~,V] = svd(B);
M = diag([1,1,det(U)*det(V)]);
R = U*M*V';