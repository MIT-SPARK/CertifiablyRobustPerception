function lams = mineig(C)
if ~iscell(C)
    lams = min(eig(C));
else
    lams = zeros(length(C),1);
    for i = 1:length(C)
        lams(i) = min(eig(C{i})); 
    end
end
end