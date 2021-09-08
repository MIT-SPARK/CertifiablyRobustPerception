function shape = combine_shapes(shapes,c)
assert(size(shapes,3) == size(c,1), 'combine shapes dimension wrong!')
K       = size(c,1);
shape   = zeros(size(shapes,1),size(shapes,2));
for k = 1:K
    shape = shape + squeeze( shapes(:,:,k) * c(k) );
end
end