function mat_K = vec_K_to_mat_K(vec_K)
mat_K           = eye(3);
mat_K(1,1)      = vec_K(1);
mat_K(2,2)      = vec_K(2);
mat_K(1,3)      = vec_K(3);
mat_K(2,3)      = vec_K(4);
end