function SDP = relax_single_rotation_averaging(problem)
%% Apply a sparse second-order relaxation to single rotation averaging
%% Depending on multivariate polynomial package SPOT
%% Heng Yang, June 29, 2021

fprintf('\n===================================================================')
fprintf('\nApplying SDP relaxation to single rotation averaging problem')
fprintf('\n===================================================================\n')
t0              = tic;

%% define POP variables
N               = problem.N;
noiseBoundSq    = problem.noiseBoundSq;
R_measurements  = problem.R_measurements;

nrPrimalVars    = 9 + N;
p               = msspoly('p',nrPrimalVars);
r               = p(1:9); 
col1            = r(1:3);
col2            = r(4:6);
col3            = r(7:9);
theta           = p(10:nrPrimalVars);

%% compute the cost function
residuals = {};
for i = 1:N 
    ri               = reshape(R_measurements(:,:,i),[9,1]);
    residuals{end+1} = (6-2*ri'*r) / noiseBoundSq;
end
f_cost = 0;
for i = 1:N 
    f_cost = f_cost + (1+theta(i))/2 * residuals{i} + (1-theta(i))/2 * 1.0;
end

%% Define the equality constraints
h_r = [1.0-col1'*col1;...
       1.0-col2'*col2;...
       1.0-col3'*col3;... % columns unit length
       col1'*col2;...
       col2'*col3;...
       col3'*col1;... % colums orthogonal
       cross(col1,col2) - col3;...
       cross(col2,col3) - col1;...
       cross(col3,col1) - col2]; % columns righthandedness
   
h_theta = [];
for i = 1:N 
    h_theta =[h_theta; 1-theta(i)^2];
end

%% Formulate the sparse second-order relaxation
basis_p         = [1;r;theta;mykron(theta,r)];
n               = length(basis_p);
basis_r         = monomials(theta,0:2);
basis_theta     = monomials(r,0:2);

fprintf('Computing localizing and moment polynomials ...')
time_start  = tic;
pop = [mykron(basis_r,h_r);...
       mykron(basis_theta,h_theta);...
       mykron(basis_p,basis_p);
       f_cost];
[~,degmat,coef_all] = decomp(pop);
coef_all            = coef_all';
time_prep   = toc(time_start);
fprintf(' Done in %g seconds.\n',time_prep);

dim_loc = length(basis_r)*length(h_r) + length(basis_theta)*length(h_theta);

%% generate standard SDP data from degmat and coefficients
[Acell,b,C] = poly2SDP(n,degmat,coef_all,dim_loc);

%% svec in sdpt3 format and output standard data
fprintf('Converting Acell to svec...')
blk{1,1}        = 's';
blk{1,2}        = n;
At              = sparsesvec(blk,Acell);

SDP.blk = blk;
SDP.At  = {At};
SDP.n   = n;
SDP.m   = length(b);
SDP.C   = {C};
SDP.b   = b;
SDP.M   = 4*N + 4;

tf    = toc(t0);
fprintf('\nDone in %g seconds.\n',tf);
fprintf('===================================================================\n')

end