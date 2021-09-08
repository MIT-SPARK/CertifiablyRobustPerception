function tranError = getTranslationError(t_gt,t_est)
t_gt = t_gt(:);
t_est = t_est(:);
tranError = norm(t_gt-t_est);
end