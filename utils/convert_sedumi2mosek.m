function prob = convert_sedumi2mosek(A, b, c, K)
%% Credit to: https://github.com/unc-optimization/SieveSDP

m = length(b);  % Number of constraints
% Make sure that sizes of A, b and c are consistent
if size(A, 2) ~= m
    A = A';
end
if size(b, 2) ~= 1
    b = b';
end
if size(c, 2) ~= 1
    c = c';
end

% Get the number of free variables
n_fre = 0;
if isfield(K, 'f') && ~isempty(K.f) && (K.f > 0)
    n_fre = K.f;
else
    K.f = 0;
end

% Get the number of positive variables
n_pos = 0;
if isfield(K, 'l') && ~isempty(K.l) && (K.l > 0)
    n_pos = K.l;
else
    K.l = 0;
end

% Get the orders of SDP variables
if isfield(K, 's') && ~isempty(K.s) && (K.s(1) > 0)
    nj = K.s;
else
    nj = [];
end
n_sdp = length(nj);

% prepare for blocks
njs = nj.^2;
Ainds = zeros(n_sdp + 1, 1);
for j = 1:n_sdp
   Ainds(j + 1) = Ainds(j) + njs(j);
end
prob.bardim = nj;

% Convert c
prob.c = c(1:(n_fre + n_pos));
c(1:(n_fre + n_pos)) = [];
I = find(c);
count = length(I);
subj = zeros(1, count);
subk = zeros(1, count);
subl = zeros(1, count);
val = c(I)';
if ~isempty(I)
I(count + 1) = Ainds(end) + 1;
k = 1;
while I(1) > Ainds(k + 1)
    k = k + 1;
end
start = 1;
stop = 1;
for vv = 2:count
    j = k;
    while I(vv) > Ainds(j + 1)
        j = j + 1;
    end
    if j > k
        [row, col] = ind2sub([nj(k), nj(k)], I(start:stop) - Ainds(k));
        subj(start:stop) = k*ones(1, stop - start + 1);
        subk(start:stop) = row;
        subl(start:stop) = col;
        start = vv;
        k = j;
    end
    stop = vv;
end
[row, col] = ind2sub([nj(k), nj(k)], I(start:stop) - Ainds(k));
subj(start:stop) = k*ones(1, stop - start + 1);
subk(start:stop) = row;
subl(start:stop) = col;
I = find(subk < subl);
subj(I) = [];
subk(I) = [];
subl(I) = [];
val(I) = [];
end
prob.barc.subj = subj;
prob.barc.subk = subk;
prob.barc.subl = subl;
prob.barc.val = full(val);

% Convert b
prob.blc = b;
prob.buc = prob.blc;

% Convert A
prob.a = A(1:(n_fre + n_pos), :)';
A(1:(n_fre + n_pos), :) = [];
subi = cell(1, m);
subj = cell(1, m);
subk = cell(1, m);
subl = cell(1, m);
val = cell(1, m);
for i = 1:m
    ai = A(:, i);
    I = find(ai);
    count = length(I);
    subi{i} = i*ones(1, count);
    subj{i} = zeros(1, count);
    subk{i} = zeros(1, count);
    subl{i} = zeros(1, count);
    val{i} = ai(I)';
    if isempty(I)
        continue;
    end
    I(count + 1) = Ainds(end) + 1;
    k = 1;
    while I(1) > Ainds(k + 1)
        k = k + 1;
    end
    start = 1;
    stop = 1;
    for vv = 2:count
        j = k;
        while I(vv) > Ainds(j + 1)
            j = j + 1;
        end
        if j > k
            [row, col] = ind2sub([nj(k), nj(k)], I(start:stop) - Ainds(k));
            subj{i}(start:stop) = k*ones(1, stop - start + 1);
            subk{i}(start:stop) = row;
            subl{i}(start:stop) = col;
            start = vv;
            k = j;
        end
        stop = vv;
    end
    [row, col] = ind2sub([nj(k), nj(k)], I(start:stop) - Ainds(k));
    subj{i}(start:stop) = k*ones(1, stop - start + 1);
    subk{i}(start:stop) = row;
    subl{i}(start:stop) = col;
    I = find(subk{i} < subl{i});
    subi{i}(I) = [];
    subj{i}(I) = [];
    subk{i}(I) = [];
    subl{i}(I) = [];
    val{i}(I) = [];
end
prob.bara.subi = horzcat(subi{1:m});
prob.bara.subj = horzcat(subj{1:m});
prob.bara.subk = horzcat(subk{1:m});
prob.bara.subl = horzcat(subl{1:m});
prob.bara.val = full(horzcat(val{1:m}));

% Convert x
prob.blx = [-inf(1, n_fre), zeros(1, n_pos)];
prob.bux = [];
end