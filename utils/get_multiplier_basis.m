function basis = get_multiplier_basis(x,v,h,verbose)
%% x is a vector of polynomial variables
%% v is the sparse monomial basis of the moment matrix: M = v*v';
%% h is an equality constraint
%% Find a set of basis monomials such that h*basis \in v
%% Heng Yang, July 05, 2021
if nargin < 4
    verbose = 1;
end

degv    = deg(v);
degh    = deg(h);
degv2   = 2*degv;

if degh > degv2
    error('The equality constraint h cannot be described by the basis v.');
end

if degh == degv2
    basis = [1];
else

Mflat        = mykron(v,v);
[~,degmat,~] = decomp(Mflat);

dense               = monomials(x,0:degv2-degh);
[~,degmat1,coef]    = decomp(h*dense);

lia          = ismember(degmat1,degmat,'rows');

valid        = ones(length(dense),1);

for i = 1:length(lia)
    if ~lia(i)
        [rows,~,~]  = find(coef(:,i));
        valid(rows) = 0;
    end
end

basis        = dense(logical(valid));
end

if verbose > 0
    fprintf('Found a basis of %d elements.\n',length(basis));
end
end