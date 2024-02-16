function R2 = coef_det(y,f)
% x - x values
% y - experimental values
% f - modeled values

y_m = mean(y);

SS_res = sum((y-f)  .^2);
SS_tot = sum((y-y_m).^2);

R2 = 1-SS_res/SS_tot;
end