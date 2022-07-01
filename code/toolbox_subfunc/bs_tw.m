function fval = bs_tw(X,w)
% calculate Brier score for multiple periods

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
    
    
    temp_I = zeros(size(temp_x,1),1);
    temp_I((temp_y>temp_x(:,1)) & (temp_y<=temp_x(:,2))) = 1;
    
    fval(t) = mean((temp_p/100-temp_I).^2);
    
end





