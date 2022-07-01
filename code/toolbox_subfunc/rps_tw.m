function fval = rps_tw(X,w)
% calculate Ranked prob score for multiple periods

% X is structure
% w is mixture weight (k by 1)
% y is realization

nt = size(X,1);
fval = zeros(nt,1);
for t=1:1:nt
    temp_XX = X(t);
    
    temp_p = temp_XX.hist_fixed_nozero'*w;
    temp_x = temp_XX.histx_fixed;
    temp_y = temp_XX.actual;
    
    fval(t) = sum( (cumsum(temp_p)/100 - (temp_y<=temp_x(:,2))).^2 );
end





