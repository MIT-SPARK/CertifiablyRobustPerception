function R = rand_rotation(varargin)
% generate random rotation matrix

params = inputParser;
params.CaseSensitive = false;

params.addParameter('RotationBound',2*pi,...
    @(x) 0.0<=x && x<=2*pi);

params.parse(varargin{:});

RotationBound = params.Results.RotationBound;

angle = RotationBound*rand - RotationBound/2;
axis  = randn(3,1);
axis  = axis / norm(axis);
R     = axang2rotm([axis' angle]);
end