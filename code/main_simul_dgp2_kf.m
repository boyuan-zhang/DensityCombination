% Simulation for Simplex with L1, Entropy regularization

% This is the DGP2 in the paper
% Estimate using Kalman Filter

clc; clear all; close all;

workpath = pwd;
savepath = [workpath, filesep, 'graphics'];
latexpath = [workpath, filesep, 'latex'];
resultpath = [workpath, filesep, 'simul_results'];

addpath(genpath('toolbox_plot'));
addpath(genpath('toolbox_subfunc'));

chk_dir(savepath)
chk_dir(latexpath)
chk_dir(resultpath)

rng(1001)

%% Setting

nsimul = 10000;

K = 20; % # of forecasters
T = 20; % number of time periods
T2 = 1; % OOS evaluation (this experiment only makes sense when

lam_s = ones(K,1); %forecaster's conditional mean follows a factor structure
sig_s = 1; %forecaster's conditional variance (in mean)

rho_x = 0.9;
sig_x = 1;

sig_y = 0.5; %dgp's conditional variance

%% Shrinkage setting
ngrid = 10;
grid_lam_2 = [linspace(1e-15, 0.2, ngrid/2)';linspace(0.3, 20, ngrid/2)'];
grid_lam_3 = [linspace(1e-15, 10, ngrid/2)';linspace(15, 10000, ngrid/2)'];

%% average subset setting
Nmax = round(K);

%% Loop starts here
mat_score_ind = zeros(nsimul,K);
mat_score_avg = zeros(nsimul,1);
mat_score_reg = zeros(nsimul,2);
mat_n_reg = zeros(nsimul,2);
mat_n_reg2 = zeros(nsimul,ngrid);
mat_n_reg3 = zeros(nsimul,ngrid);

mat_score_reg2 = zeros(nsimul,ngrid); %KL

% mixture of the best N
mat_score_mixBestN = zeros(nsimul,Nmax); 
mat_n_mixBestN = zeros(nsimul,Nmax);

% mixture of the best N <= Nmax
mat_score_mixBestLessN = zeros(nsimul,Nmax);
mat_n_mixBestLessN = zeros(nsimul,Nmax);

% best N-misture
mat_score_BestNmix = zeros(nsimul,Nmax); 
mat_n_BestNmix = zeros(nsimul,Nmax);

% best less than N<= Nmax-misture
mat_score_BestLessNmix = zeros(nsimul,Nmax); 
mat_n_BestLEssNmix = zeros(nsimul,Nmax);

% individual for comparison
mat_score_individual = zeros(nsimul, 5);

parfor simulind = 1:nsimul
    
    simulind
    rng(simulind)

    %% Data - Conditional Mean
    x0 = sqrt( sig_x^2 / (1-rho_x^2) ) * randn(1,1);
    
    yt = zeros(T+T2,1);
    xt = zeros(T+T2,1); %conditional mean of dgp
    st = zeros(T+T2,K); %conditional mean of forecasts
    for t=1:(T+T2)
        x0 = rho_x*x0 + sig_x*randn(1,1);
        xt(t,1) = x0;
        %         st(t,:) = lam_s'*x0 + sig_s*randn(1,K);
        st(t,:) = lam_s'*x0 + [sig_s*ones(1,K/2), 5*sig_s*ones(1,K/2)].*randn(1,K);
        
        yt(t,1) = x0 + sig_y*randn(1,1);
    end
    
	% ----------------------------------------------------
    % run KF to obtain E[xt|zt,zt-1,...]
    x0 = 0;
	P0 = sig_x^2 / (1-rho_x^2);
    % ----------------------------------------------------
    %  y_{t} = A + B s_{t} + u
    %  s_{t} = Phi s_{t-1} + R*ep
    %  var(u) ~ H
    %  var(ep) ~ S2
    % ----------------------------------------------------
    xt_up = zeros(T+T2,1);
    for k = 1:K
        if k <= K/2
            [s_up, P_up, loglik] = kalman_filter_ini(x0,P0,0,1,rho_x,1,sig_s^2,sig_x^2, st(:,k));
        else
            [s_up, P_up, loglik] = kalman_filter_ini(x0,P0,0,1,rho_x,1,(5*sig_s)^2,sig_x^2, st(:,k));
        end
        xt_up(:,k) = s_up;
    end

    % height of the density
    %ft = normpdf(yt, st, sig_y); % original st = zt
    ft = normpdf(yt, xt_up, sig_y); % revised xt_up = E[xt|zt,zt-1,zt-2...]
    
    %% Weight computation
    data_est = ft(1:T,:);
    
    options = optimoptions('fmincon','Display','off');
    
    % Optimal linear opinion
    fun = @(b) sum(-log(data_est*b));
    bini = (1/K)*ones(K,1);
%     options.Display = 'off';
    [b_reg,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(fun,bini,[],[],ones(1,K),1,zeros(K,1),ones(K,1),[],options);
    
    % Dirichlet shrinkage
    temp_b_reg2 = [];
    for gg=1:ngrid
        lam = grid_lam_2(gg);
        fun = @(b) sum(-log(data_est*b)) - lam*sum(log(b));
        bini = (1/K)*ones(K,1);
%         options.Display = 'off';
        [b_reg2,FVAL2,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(fun,bini,[],[],ones(1,K),1,zeros(K,1),ones(K,1),[],options);
        temp_b_reg2 = [temp_b_reg2, b_reg2];
    end
    
    % L2 shrinkage
    temp_b_reg3 = [];
    for gg=1:ngrid
        lam = grid_lam_3(gg);
        fun = @(b) sum(-log(data_est*b))  + lam*sum( (b-1/K).^2 );
        bini = (1/K)*ones(K,1);
%         options.Display = 'off';
        [b_reg3,FVAL3,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(fun,bini,[],[],ones(1,K),1,zeros(K,1),ones(K,1),[],options);
        temp_b_reg3 = [temp_b_reg3, b_reg3];
    end
    
    % Equal weight
    b_reg4 = 1/K *ones(K,1);
    
    %     % Elastic-Net shrinkage (zero)
    %     lam = 70;
    %     fun = @(b) sum(-log(data_est*b))  + lam*sum( (b-0).^2 );
    %     bini = (1/K)*ones(K,1);
    %     options.Display = 'off';
    %     [b_reg5,FVAL5,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(fun,bini,[],[],ones(1,K),1,zeros(K,1),ones(K,1),[],options);
    
    %% Averaging methods - mixture of best N
    % select the best N forecastor based on rps, then take the averge of their
    % density forecasts.
    
    mat_ls = - log(data_est); % data_est is a T x K matrix containing log score for each forecastor
    
    % ls (sorted), smaller is better
    temp_ls_ind = sum(mat_ls,1);
    [~,b] = sort(temp_ls_ind, 'ascend');
    
    % mixture of best N
    ls_mixbestN = zeros(Nmax, 1); % record ls for each bestN
    id_mixbestN = nan(K, Nmax); % record the selected forecasts. Row: forecastors; Column: bestN
    
    for Nind = 1:1:Nmax
        
        best_density = mean(data_est(:,b(1:Nind)),2);
        ls_mixbestN(Nind,1) = sum(-log(best_density));
        
        id_mixbestN(b(1:Nind), Nind) = 1;
    end
    
    %% Averaging methods - Best N-Mixture
    % consider all possible combination of N forecastor, evaluate the score for
    % each combo and select the best.
    
    mat_bestNmix_ls = nan*ones(Nmax,1);
    mat_bestNmix_set = cell(Nmax,1);
    
    for Nind = 1:1:Nmax
        
        mixset = nchoosek( (1:1:K)', Nind ); % all possible combinations
        mixset_n = size(mixset,1);
        
        % historical performance of all N-mixtures
        temp_p = zeros(T, mixset_n); % T x # of combinations
        for t = 1:1:T
            temp_p(t,:) = mean(reshape(data_est(t,mixset),[],Nind),2); % calculate the avg pdf for each bin and each combination
        end
        temp_ls = sum( - log(temp_p), 1)';
        
        % ls (sorted), smaller is better
        [a,b] = sort(temp_ls, 'ascend');
        
        % record N-best-mix
        mat_bestNmix_ls(Nind, 1) = a(1);
        mat_bestNmix_set{Nind, 1} = mixset(b(1),:)';
        
    end
    
    %% OOS evaluation
    n_cut = (1/K)*1e-3;
    
    % evaluation sample
    data_eval = ft(T+1:end,:);
    
    % individual score
    score_ind = -log(data_eval);
    
    % regression score
    score_reg = [-log(data_eval*b_reg),-log(data_eval*b_reg4)];
    
    score_reg2 = -log(data_eval*temp_b_reg2);
    score_reg3 = -log(data_eval*temp_b_reg3);
    
    % number of selected forecasters
    n_reg = [sum(b_reg>n_cut), sum(b_reg4>n_cut)];
    n_reg2 = sum(temp_b_reg2>n_cut,1);
    n_reg3 = sum(temp_b_reg3>n_cut);

    
    % Averaging methods - mixture of best N
    score_mixBestN = zeros(Nmax, 1);
    n_mixBestN = 1:Nmax;
    
    for Nind = 1:1:Nmax
        temp_id = find(id_mixbestN(:,Nind) == 1);
        score_mixBestN(Nind) = -log( mean(data_eval(:,temp_id)) );
    end
    
    % Averaging methods - mixture of best N <=  Nmax
    score_mixBestLessN = zeros(Nmax, 1);
    n_mixBestLessN     = zeros(Nmax, 1);
    
    for Nind = 1:1:Nmax
        temp_ls = ls_mixbestN(1:Nind);
        [~,b]   = sort(temp_ls,'ascend');
        temp_id = find(id_mixbestN(:,b(1)) == 1);
        
        score_mixBestLessN(Nind) = -log( mean(data_eval(:,temp_id)) );
        n_mixBestLessN(Nind) = numel(temp_id);
    end
    
    % Averaging methods - Best N-Mixture
    score_BestmixN = mat_bestNmix_ls;
    n_BestmixN = 1:Nmax;
    
    for Nind = 1:1:Nmax
        temp_set = mat_bestNmix_set{Nind};
        score_BestmixN(Nind) = -log( mean(data_eval(:,temp_set)) );
    end
    
    % Averaging methods - Best N-Mixture <= Nmax
    score_BestmixLessN = nan(Nmax, 1);
    n_BestmixLessN = nan(Nmax, 1);
    
    for Nind = 1:1:Nmax
        temp_ls = mat_bestNmix_ls(1:Nind);
        [a,b] = sort(temp_ls,'ascend');
        temp_set = mat_bestNmix_set{b(1)};
        
        score_BestmixLessN(Nind,1) = -log( mean(data_eval(:,temp_set)) );
        n_BestmixLessN(Nind,1) = numel(temp_set);
    end
    
    % individual for comparison
    temp = sort(score_ind, 'descend');
    score_comparison = temp(flip([1, 5, K/2, K-4, K]));

    %% Store
    mat_score_ind(simulind,:)  = score_ind;
    mat_score_avg(simulind,:)  = mean(score_ind);
    mat_score_reg(simulind,:)  = score_reg;
    mat_score_reg2(simulind,:) = score_reg2;
    mat_score_reg3(simulind,:) = score_reg3;
    mat_n_reg(simulind,:)  = n_reg;
    mat_n_reg2(simulind,:) = n_reg2;
    mat_n_reg3(simulind,:) = n_reg3;
    
    mat_score_individual(simulind,:) = score_comparison;
    
    % mixture of best N
    mat_score_mixBestN(simulind,:) = score_mixBestN;
    mat_n_mixBestN(simulind,:) = n_mixBestN;
    
    % mixture of the best N <= Nmax
    mat_score_mixBestLessN(simulind,:) = score_mixBestLessN;
    mat_n_mixBestLessN(simulind,:) = n_mixBestLessN;
    
    % best N-misture
    mat_score_BestNmix(simulind,:) = score_BestmixN; 
    mat_n_BestNmix(simulind,:) = n_BestmixN;
    
    % best less than N<= Nmax-misture
    mat_score_BestLessNmix(simulind,:) = score_BestmixLessN; 
    mat_n_BestLessNmix(simulind,:) = n_BestmixLessN;
    
end

%% Summary
disp('results ... ');
res_m1 = mean(mat_score_reg,1)
mean(mat_n_reg,1)

res_m2 = mean(mat_score_reg2,1)
mean(mat_n_reg2,1)

res_m3 = mean(mat_score_reg3,1)
mean(mat_n_reg3,1)

mat_score_avg(mat_score_avg == -Inf) = min( mat_score_avg(mat_score_avg ~= -Inf) );
res_avg = mean(mat_score_avg,1)

res_mixBestN = mean(mat_score_mixBestN,1)
mean(mat_n_mixBestN,1)
[res_mixBestN_max,res_mixBestN_n] = min( res_mixBestN );

res_mixBestLessN = mean(mat_score_mixBestLessN,1)
[res_mixBestLessN_max,b] = min( res_mixBestLessN );
temp = mean(mat_n_mixBestLessN,1)
res_mixBestLessN_n = temp(b)

res_BestNmix = mean(mat_score_BestNmix,1)
mean(mat_n_BestNmix,1)
[res_BestNmix_max,res_BestNmix_n] = min( res_BestNmix );

res_BestLessNmix = mean(mat_score_BestLessNmix,1)
[res_BestLessNmix_max,b] = min( res_BestLessNmix );
temp = mean(mat_n_BestLessNmix,1)
res_BestLessNmix_n = temp(b)

mat_score_individual(mat_score_individual == Inf) = max(mat_score_individual(mat_score_individual ~= Inf));
res_comparison = mean(mat_score_individual,1)

%% Save
cd(resultpath)
savefilename = ['dgp2_K',num2str(K), 'T',num2str(T), '_KF.mat'];
save(savefilename);
cd(workpath)

%% Figure
% part 1
fig = figure(1);
position = [2,2,12,7];
set(fig, 'color', 'w');
set(fig, 'units', 'inches');
set(fig, 'position', position);

plot(res_m1(2)*ones(ngrid,1), 'linewidth', 3, 'color','black') % equal weight
hold on
plot(res_m1(1)*ones(ngrid,1), 'linewidth', 3, 'color','black') % Simplex
plot(res_m2, 'linewidth', 3, 'color','black') % Simplex + Entropy
plot(res_m3, 'linewidth', 3, 'color','black') % Simplex + Ridge
plot(min(res_BestNmix)*ones(ngrid,1), 'linewidth', 3, 'color','black')
hold off
ylim([1.05, 1.4]); yticks(1.0:0.1:1.4);
xlim([0.5,ngrid+0.5]); xticks(1:1:10);
xlabel('Penalty Strength');
ylabel('Average Log Score');

text(2,1.36,'Simple Average', 'fontsize', 20)
text(6.9,1.26,{'Simplex \newline+Entropy'}, 'fontsize', 20)
text(5.1,1.25,{'Simplex \newline+Ridge'}, 'fontsize', 20)
text(8.5,1.24,'Simplex', 'fontsize', 20)
text(8.5,1.14,'N-Average', 'fontsize', 20)

set(gca,'linewidth', 1, 'fontsize', 20);

box off
ax=gca;
axes('position',ax.Position,'box','on','ytick',[],'xtick',[],'color','none','LineWidth',1)

title('DGP 2', 'FontWeight','normal','FontSize',22);

cd(savepath);
saveas(fig,['fig_simul_dgp2_part1_K',num2str(K), '_T',num2str(T),'_KF.png']);
cd(workpath);

close all


% part 2
fig = figure(2);
setmyfig(fig, [2,2,8,6]);
plot(res_m1(2)*ones(Nmax,1), 'linewidth', 3)
hold on
plot(res_m1(1)*ones(Nmax,1), 'linewidth', 3)
plot(res_mixBestN, 'linewidth', 3)
plot(res_mixBestLessN, 'linewidth', 3)
plot(res_BestNmix, 'linewidth', 3)
plot(res_BestLessNmix, 'linewidth', 3)
hold off
xlim([1,Nmax+1]);
xlabel('Nmax');
l = legend('Simple Average', 'Simplex',...
    'Mixture of Best N','Mixture of Best N $\le$ Nmax',...
    'Best N-Mixture','Best N-Mixture, N $\le$ Nmax',...
    'Location','northeast','Interpreter','latex');
set(l, 'fontsize', 16);
set(gca,'linewidth',2,'fontsize', 20);
title(['DGP 2 (T=',num2str(T),')']);
grid on

cd(savepath);
saveas(fig,['fig_simul_dgp2_part2_K',num2str(K), '_T',num2str(T),'_KF.png']);
cd(workpath);

close all

%% save "data" in csv file
cd(resultpath)

% part 1
output_table = table((1:ngrid)', res_m1(2)*ones(ngrid,1), res_m1(1)*ones(ngrid,1), res_m2', res_m3',...
                'VariableNames',{'shrinkage','simple average', 'simplex', 'simplex+entropy', 'simplex+ridge'});
write(output_table, ['simul_result2plot_dgp2_part1_K',num2str(K), '_T',num2str(T),'_KF.csv'])



% part 2
output_table = table((1:Nmax)', res_m1(2)*ones(Nmax,1), res_m1(1)*ones(Nmax,1), res_mixBestN',...
                     res_mixBestLessN', res_BestNmix',res_BestLessNmix',...
                     'VariableNames',{'Nmax','Simple Average', 'Simplex',...
                     'Mixture of Best N','Mixture of Best N<=Nmax',...
                      'Best N-Mixture','Best N-Mixture, N<=Nmax'});
write(output_table, ['simul_result2plot_dgp2_part2_K',num2str(K), '_T',num2str(T),'_KF.csv'])

cd(workpath)

%% generate latex

% simplex and average
res_simplex_score = res_m1(1);
res_simplex_n     = mean(mat_n_reg(:,1),1);
res_average_score = res_m1(2);
res_average_n     = K;

% ridge
[a,b] = min(res_m3);
res_reg_L2_score  = a;
res_reg_L2_n      = mean(mat_n_reg3(:,b),1);
res_reg_L2_lambda = grid_lam_3(b);

% entropy
[a,b] = min(res_m2);
res_dirichlet_score  = a;
res_dirichlet_n      = mean(mat_n_reg2(:,b),1);
res_dirichlet_lambda = grid_lam_2(b);

% best mixture N
res_bestNmix_score = res_BestNmix;
res_bestNmix_n     = mean(mat_n_BestNmix,1);

% best mixture N<=Nmax
res_bestmixLessN_score = res_BestLessNmix;
res_bestmixLessN_n     = mean(mat_n_BestLessNmix,1);


%% --- Table 1: perfornmance of several methods


%--- Part 1: make matlab table

% sub-panel 1
panel1_name   = {'Simplex'; 'Simplex + Ridge'; 'Simplex + Entropy'};
panel1_score  = [res_simplex_score; res_reg_L2_score; res_dirichlet_score];
panel1_N      = [res_simplex_n; res_reg_L2_n; res_dirichlet_n];
panel1_lamdba = {'NA'; res_reg_L2_lambda; res_dirichlet_lambda};

panel1 = table(panel1_name, panel1_score, panel1_N, panel1_lamdba, 'VariableNames',{'title1','R','#','lambda'});

% sub-panel 2
N_list1 = [1:10 15 20];
N_ID_p1 = cell(numel(N_list1),1);
prefix = '\quad $N = ';
for n = 1:numel(N_list1)
    N_ID_p1{n,1} = [prefix num2str(N_list1(n)) '$'];
end

N_list2 = [2:3 5 10 15 20];
N_ID_p2 = cell(numel(N_list2),1);
prefix = 'Best $\leq ';
for n = 1:numel(N_list2)
    N_ID_p2{n,1} = [prefix num2str(N_list2(n)) '$-Average'];
end

panel2_name   = [N_ID_p1; N_ID_p2];
panel2_score  = [res_bestNmix_score(N_list1) res_bestmixLessN_score(N_list2)]';
panel2_N      = [res_bestNmix_n(N_list1) res_bestmixLessN_n(N_list2)]';
panel2_lambda = repmat('NA',numel(N_list1) + numel(N_list2),1);

panel2 = table(panel2_name, panel2_score, panel2_N, panel2_lambda, 'VariableNames',{'title1','R','#','lambda'});


% sub-panel 3
panel3_name   = {'Best'; '75\%'; 'Median'; '25\%'; 'Worst'};

panel3_score  = res_comparison';
panel3_N      = ones(numel(panel3_name),1);
panel3_lambda = repmat('NA',5,1);

panel3 = table(panel3_name, panel3_score, panel3_N, panel3_lambda, 'VariableNames',{'title1','R','#','lambda'});


%--- Part 2: generate latex table

fid = fopen([latexpath '/table_DGP2_ls_all_K',num2str(K), '_T',num2str(T),'_KF.tex'],'wt');

% header and others
fprintf(fid, '%s \n', '\begin{table}[htp]');
fprintf(fid, '%s \n', '\begin{center}');
fprintf(fid, '%s \n', '\caption{Ex Post Optimal Average Log Scores, DGP2}');
fprintf(fid, '%s \n', '\begin{tabular}{lrrr}');
fprintf(fid, '%s \n', '\toprule');fprintf(fid, '%s \n', '\toprule');
fprintf(fid, '%s \n', 'Regularization group & $ L $ & \# & $\lambda^*$ \\');
fprintf(fid, '%s \n', '\midrule');

% sub-panel 1: regularization
for n = 1:size(panel1,1)
    fprintf(fid, '%s ', char(panel1{n,1}));
    fprintf(fid, '%s ', '&');
    fprintf(fid, '%.2f ', panel1{n,2});
    fprintf(fid, '%s ', '&');
    fprintf(fid, '%.2f ', panel1{n,3});
    fprintf(fid, '%s ', '&');
    
    if n == 1
       fprintf(fid, '%s ', char(panel1{n,4}));
    else
       fprintf(fid, '%.2f ', cell2mat(panel1{n,4}));
    end
    
    fprintf(fid, '%s \n', '\\');
end
fprintf(fid, '%s \n', '\midrule'); fprintf(fid, '%s \n', '\midrule');

% sub-panel 2: mixture
fprintf(fid, '%s \n', 'Subset Averages & $ L $ & $\#$ & $\lambda^*$ \\');
fprintf(fid, '%s \n', '\midrule');
fprintf(fid, '%s \n', 'Best $N$-Average: &  &  & \\');

for n = 1:size(panel2,1)
    fprintf(fid, '%s ', char(panel2{n,1}));
    fprintf(fid, '%s ', '&');
    fprintf(fid, '%.2f ', panel2{n,2});
    fprintf(fid, '%s ', '&');
    fprintf(fid, '%.2f ', panel2{n,3});
    fprintf(fid, '%s ', '&');
    fprintf(fid, '%s ', panel2{n,4});
    fprintf(fid, '%s \n', '\\');
end
fprintf(fid, '%s \n', '\midrule'); fprintf(fid, '%s \n', '\midrule');

% sub-panel 3: indiviudal comparisons
fprintf(fid, '%s \n', 'Comparisons & $ L $ & $\#$ & $\lambda^*$ \\');
fprintf(fid, '%s \n', '\midrule');
for n = 1:size(panel3,1)
    fprintf(fid, '%s ', char(panel3{n,1}));
    fprintf(fid, '%s ', '&');
    fprintf(fid, '%.2f ', panel3{n,2});
    fprintf(fid, '%s ', '&');
    fprintf(fid, '%.f ', panel3{n,3});
    fprintf(fid, '%s ', '&');
    fprintf(fid, '%s ', char(panel3{n,4}));
    fprintf(fid, '%s \n', '\\');
end
fprintf(fid, '%s \n', '\midrule');

% fprintf(fid, '%s ', 'Mean');
% fprintf(fid, '%s ', '&');
% fprintf(fid, '%.2f ', res_avg);
% fprintf(fid, '%s ', '&');
% fprintf(fid, '%.f ', 1);
% fprintf(fid, '%s ', '&');
% fprintf(fid, '%s ', 'NA');
% fprintf(fid, '%s \n', '\\');
% fprintf(fid, '%s \n', '\midrule');

% sub-panel 4: Simple average
fprintf(fid, '%s ', 'Simple Average');
fprintf(fid, '%s ', '&');
fprintf(fid, '%.2f ', res_average_score);
fprintf(fid, '%s ', '&');
fprintf(fid, '%.f ', K);
fprintf(fid, '%s ', '&');
fprintf(fid, '%s ', 'NA');
fprintf(fid, '%s \n', '\\');
fprintf(fid, '%s \n', '\bottomrule'); fprintf(fid, '%s \n', '\bottomrule');

% footer
fprintf(fid, '%s \n', '\end{tabular}');
fprintf(fid, '%s \n', '\end{center}');
fprintf(fid, '%s \n', '\end{table}');


%% Latex table for Nmixture and mixture of N

% N_list = [1:10 15 20];
% 
% N_ID = strings;
% prefix = '$N = ';
% for n = 1:numel(N_list)
%     N_ID(n,:) = [prefix num2str(N_list(n)) '$'];
% end
% 
% Nmax_ID = strings;
% prefix = '$N_{max} = ';
% for n = 1:numel(N_list)
%     Nmax_ID(n,:) = [prefix num2str(N_list(n)) '$'];
% end
% 
% %%% (1)
% % average-best N
% output_table = table(N_ID, round(res_mixBestN(N_list), 2)', round(mean(mat_n_mixBestN(:,N_list),1), 2)',...
%     'VariableNames',{'title1','L','#'});
% 
% fid = fopen([latexpath '/table_DGP2_averageBestN_K',num2str(K), '_T',num2str(T),'.tex'],'wt');
% 
% for n = 1:numel(N_list)
%     fprintf(fid, '%s %s %.2f %s %.f %s %s %s \n', ...
%         [output_table{n,1} '&' output_table{n,2} '&' output_table{n,3} '&' 'NA' '\\'] );
% end
% 
% %%% (2)
% % average-best N <= Nmax
% output_table = table(Nmax_ID, round(res_mixBestLessN(N_list), 2)', round(mean(mat_n_mixBestLessN(:,N_list),1), 2)',...
%     'VariableNames',{'title1','L','#'});
% 
% fid = fopen([latexpath '/table_DGP2_averageBestLessN_K',num2str(K), '_T',num2str(T),'.tex'],'wt');
% 
% for n = 1:numel(N_list)
%     fprintf(fid, '%s %s %.2f %s %.2f %s %s %s \n', ...
%         [output_table{n,1} '&' output_table{n,2} '&' output_table{n,3} '&' 'NA' '\\'] );
% end
% 
% %%% (3)
% % best N-average
% output_table = table(N_ID, round(res_BestNmix(N_list), 2)', round(mean(mat_n_BestNmix(:,N_list),1), 2)',...
%     'VariableNames',{'title1','L','#'});
% 
% fid = fopen([latexpath '/table_DGP2_BestNAverage_K',num2str(K), '_T',num2str(T),'.tex'],'wt');
% 
% for n = 1:numel(N_list)
%     fprintf(fid, '%s %s %.2f %s %.0f %s %s %s \n', ...
%         [output_table{n,1} '&' output_table{n,2} '&' output_table{n,3} '&' 'NA' '\\'] );
% end
% 
% %%% (4)
% % best N-average N <= Nmax
% output_table = table(Nmax_ID, round(res_BestLessNmix(N_list), 2)', round(mean(mat_n_BestLessNmix(:,N_list),1), 2)',...
%     'VariableNames',{'title1','L','#'});
% 
% fid = fopen([latexpath '/table_DGP2_BestNAverageLessN_K',num2str(K), '_T',num2str(T),'.tex'],'wt');
% 
% for n = 1:numel(N_list)
%     fprintf(fid, '%s %s %.2f %s %.2f %s %s %s \n', ...
%         [output_table{n,1} '&' output_table{n,2} '&' output_table{n,3} '&' 'NA' '\\'] );
% end