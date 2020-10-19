function tranErr = getTranslationError(t_est, t_gt)

tranErr = norm(t_est - t_gt);