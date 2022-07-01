function fval = rps(temp_p, temp_x, temp_y)
% Rankded prob score
% calculate Ranked prob score for multiple periods

% higher the better
% temp_p: prob 
% temp_x: range of hist
% temp_y: realization

fval = sum( (cumsum(temp_p)/100 - (temp_y<=temp_x(:,2))).^2 );