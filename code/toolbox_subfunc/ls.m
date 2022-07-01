function fval = ls(temp_p, temp_x, temp_y)
% calculate log score for single period


temp_p = temp_p/100;

temp_y_bin_pos = temp_x(:,1) < temp_y & temp_x(:,2) >= temp_y;

fval = -log(temp_p(temp_y_bin_pos));





