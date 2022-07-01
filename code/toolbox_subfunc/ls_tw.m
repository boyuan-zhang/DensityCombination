function fval = ls_tw(X,w)
% calculate log score for multiple periods

% X is structure
% w is mixture weight (k by 1)
% y is realization

% X = xxx(t0:t1);
% w = ones(nf,1)/nf;

nt = size(X,1);
fval = zeros(nt,1);

for t=1:1:nt
    temp_XX = X(t);
    
    temp_p = temp_XX.hist_fixed_nozero/100;
    temp_x = temp_XX.histx_fixed;
    temp_y = temp_XX.actual;
    
    temp_y_bin_pos = temp_x(:,1) < temp_y & temp_x(:,2) >= temp_y;
    
    temp_p2 = temp_p(:,temp_y_bin_pos)' * w;
    fval(t) = -log(temp_p2);
end





