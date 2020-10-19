function [R]=quat2rot(q,varargin)
% convert quaternion to rotation matrix
% q has the form [q1, q2, q3, q4] and does not need to be normalized, q4 is the real part

params=inputParser;
params.CaseSensitive=false;

params.addParameter('method','matrixMulti',...
	@(x) ischar(x) && any(strcmpi({'matlab','matrixMulti','monomial'},x)));

params.parse(varargin{:});

method=params.Results.method;

q=q(:)/norm(q);

switch method
	case 'matlab'
		q=quaternion([q(4),q(1:3)']);
		R=rotmat(q,'point');
	case 'matrixMulti'
		P=[1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1;
		   0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0;
		   0, 0, 1, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, -1, 0, 0;
		   0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0;
		   -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1;
		   0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0;
		   0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0;
		   0, 0, 0, -1, 0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0, 0;
		   -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
		P=sparse(P);
		R=reshape(P*reshape(q*q',[16,1]), [3,3]);
    case 'monomial'
        A = zeros(9,10);
        A(1,1)=1; A(1,2)=-1; A(1,3)=-1; A(1,4)=1;
        A(2,5)=2; A(2,10)=2;
        A(3,6)=2; A(3,9)=-2;
        A(4,5)=2; A(4,10)=-2;
        A(5,1)=-1; A(5,2)=1; A(5,3)=-1; A(5,4)=1;
        A(6,8)=2; A(6,7)=2;
        A(7,6)=2; A(7,9)=2;
        A(8,8)=2; A(8,7)=-2;
        A(9,1)=-1; A(9,2)=-1; A(9,3)=1; A(9,4)=1;
        
        alpha = [q(1)^2, q(2)^2, q(3)^2, q(4)^2, q(1)*q(2), q(1)*q(3), q(1)*q(4), q(2)*q(3), q(2)*q(4), q(3)*q(4)]';
        
        R = reshape(A*alpha, [3,3]);
        
	otherwise
		error('Unsupported conversion method.');	
end



