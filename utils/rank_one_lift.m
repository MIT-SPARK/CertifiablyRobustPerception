function X = rank_one_lift(v)

assert(iscell(v),'Rank one lift inputs a cell of vectors.')

X = cell(length(v),1);
for i = 1:length(v)
    X{i} = v{i}*v{i}';
end

end