function [score_final,n_final,b_final] = SimplexRidgeExPostLambda_opt(xxx, w, ns0, ns1,score_type, lam)

nf = size(xxx(1).hist_fixed_nozero, 1);
ns = ns1;

mat_score_b = nan(nf,ns);
mat_score_n = nan(ns,1);
mat_score_f = nan(ns,1);

for t = ns0:ns1
    
    % period selection
    t0 = max(1,t-w);
    t1 = t-1;
    
    % estimation, we are solving a min problem
    if score_type == 'R'
        fun = @(b) sum(rps_tw(xxx(t0:t1),b))  + lam*sum( (b-1/nf).^2 );
    elseif score_type == 'B'
        fun = @(b) sum(bs_tw(xxx(t0:t1),b))  + lam*sum( (b-1/nf).^2 );
    elseif score_type == 'L'
        fun = @(b) sum(ls_tw(xxx(t0:t1),b))  + lam*sum( (b-1/nf).^2 );
    end
    
    bini = (1/nf)*ones(nf,1);
    
    options = optimoptions('fmincon','Display','off');
    
    b_reg = fmincon(fun,bini,[],[],ones(1,nf),1,zeros(nf,1),ones(nf,1),[],options);
    mat_score_b(:,t) = b_reg';
    mat_score_n(t,1) = sum(b_reg > (1/nf)/100);
    
    % evaluation
    temp_p = xxx(t).hist_fixed_nozero' * b_reg;
    temp_y = xxx(t).actual;
    temp_x = xxx(t).histx_fixed;
    temp_score = 0;
    
    if score_type == 'R'
        temp_i = temp_y<=temp_x(:,2);
        temp_score = sum((cumsum(temp_p)/100 - temp_i).^2,1);
        
    elseif score_type == 'B'
        temp_i = zeros(size(temp_x,1),1);
        temp_i((temp_y>temp_x(:,1)) & (temp_y<=temp_x(:,2))) = 1;
        temp_score = mean((temp_p/100 - temp_i).^2);
        
    elseif score_type == 'L'
        temp_i = zeros(size(temp_x,1),1);
        temp_i((temp_y>temp_x(:,1)) & (temp_y<=temp_x(:,2))) = 1;
        temp_score = -log( temp_p(temp_i == 1)/100 );
    end
    
    mat_score_f(t,1) = temp_score;
    
end

% score_final = mean(mat_score_f(ns0:ns1));
% n_final= mean(mat_score_n(ns0:ns1));

score_final = mat_score_f;
n_final= mat_score_n;
b_final = mat_score_b';

end