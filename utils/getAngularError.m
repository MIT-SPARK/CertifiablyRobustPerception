function rotError = getAngularError(R_gt,R_est)

rotError = abs(acos( (trace(R_gt' * R_est)-1) / 2 ));
rotError = rad2deg( rotError );
end