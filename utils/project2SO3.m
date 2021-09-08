function R = project2SO3(M)

if size(M,1) ~= 3 || size(M,2) ~= 3
    error('project2SO3 requires a 3x3 matrix as input')
end

[U,S,V] = svd(M);
R = U*V';
if(det(R)<0)
  R = U * diag([1 1 -1]) * V';  
end

end