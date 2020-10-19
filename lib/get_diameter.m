function max_pairwise_d = get_diameter(B)
    nrFeatures = size(B,2);
    D = zeros(nrFeatures,nrFeatures);
    for i = 1:nrFeatures
        for j = i:nrFeatures
            Bi = B(:,i);
            Bj = B(:,j);
            D(i,j) = norm(Bi-Bj);
        end
    end
    max_pairwise_d = max(D(:));
end