function fval = bs(temp_p, temp_x, temp_y)
% calculate Brier score for single periods

% temp_p: prob
% temp_x: range of hist
% temp_y: realization


temp_I = zeros(size(temp_x,1),1);
temp_I((temp_y>temp_x(:,1)) & (temp_y<=temp_x(:,2))) = 1;
fval = mean((temp_p/100-temp_I).^2);

