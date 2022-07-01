% Empirical analysis of ECB data
% date: 2021/07/24
% Minchul Shin and Boyuan Zhang

clc; clear all; close all;
workpath = pwd;
datapath = '../data/';
addpath(datapath);

% savepath = [workpath, filesep, 'fig_slides'];
latexpath = [workpath, filesep, 'latex'];

addpath(genpath('toolbox_plot'));
addpath(genpath('toolbox_subfunc'));

chk_dir(latexpath)

%% load data
dname = 'infl_1y'; % infl_1y;
dnameL = ['ecbspf_',dname,'_bp_nozero']; % use the latest dataset
load([dnameL, '.mat']);
eval(['xxx = ', dnameL, ';']);

%% Set-up
ns0 = 9;  %2001Q1
ns1 = 83; %2019Q3
ns = ns1;

[nf,~] = size(xxx(1).hist);

% fixed window width
w = 20;

% score
score_type = 'L'; %'L','B','R'

%% Add one more hypothetical forecaster in the pool
version = 2; % 1 = simple average; 2 = uniform

for sind = 1:ns1
    
    temp_hist = xxx(sind).hist_fixed_nozero;
    
    if version == 1 % add simple average as the 19th guy
        temp_p = mean(xxx(sind).hist_fixed_nozero,1);
    else % add uniform as the 19th guy
        temp_p = ones(1,size(temp_hist,2));
        temp_p = temp_p/sum(temp_p) * 100;
    end
    new_hist = [temp_hist; temp_p];
    
    xxx(sind).hist_fixed_nozero = new_hist;
end

nf = nf + 1;
Nmax = nf;

%% Score measures
mat_rps = nan(ns,nf);
mat_bs  = nan(ns,nf);
mat_ls  = nan(ns,nf);

% average of 18 forecasters
mat_rps_avg1 = nan(ns,1);
mat_bs_avg1  = nan(ns,1);
mat_ls_avg1  = nan(ns,1);

% average of 19 forecasters
mat_rps_avg2 = nan(ns,1);
mat_bs_avg2  = nan(ns,1);
mat_ls_avg2  = nan(ns,1);

for sind = 1:1:ns1
    %    mat_rps(sind,:) = xxx(sind).rps';
    
    for i=1:nf
        mat_rps(sind,i) = rps(xxx(sind).hist_fixed_nozero(i,:)', xxx(sind).histx_fixed, xxx(sind).actual);
        mat_bs(sind,i) = bs(xxx(sind).hist_fixed_nozero(i,:)', xxx(sind).histx_fixed, xxx(sind).actual);
        mat_ls(sind,i) = ls(xxx(sind).hist_fixed_nozero(i,:)', xxx(sind).histx_fixed, xxx(sind).actual);
    end
    
    % average of the selected 18 (after extrapolation)
    temp_p1 = mean(xxx(sind).hist_fixed_nozero(1:end-1,:),1)';
    temp_x = xxx(sind).histx_fixed(:,:);
    temp_y = xxx(sind).actual;
    
    mat_rps_avg1(sind,:) = rps(temp_p1, temp_x, temp_y); %rps from Czado, Genieting ...
    mat_bs_avg1(sind,:) = bs(temp_p1, temp_x, temp_y);
    mat_ls_avg1(sind,:) = ls(temp_p1, temp_x, temp_y);
    
    % average of the selected 19 (with a fake forecaster)
    temp_p2 = mean(xxx(sind).hist_fixed_nozero,1)';
    mat_rps_avg2(sind,:) = rps(temp_p2, temp_x, temp_y); %rps from Czado, Genieting ...
    mat_bs_avg2(sind,:) = bs(temp_p2, temp_x, temp_y);
    mat_ls_avg2(sind,:) = ls(temp_p2, temp_x, temp_y);
end

%% Linear opinion rule

mat_simplex_score_b = nan(nf,ns);
mat_simplex_score_n = nan(ns,1);
mat_simplex_score_f = nan(ns,1);


for t = ns0:ns1
    
    % period selection
    t0 = max(1,t-w);
    t1 = t-1;
    
    % estimation
    if score_type == 'R'
        fun = @(b) sum(rps_tw(xxx(t0:t1),b)); % we are solving min problem
    elseif score_type == 'B'
        fun = @(b) sum(bs_tw(xxx(t0:t1),b)); % we are solving min problem
    elseif score_type == 'L'
        fun = @(b) sum(ls_tw(xxx(t0:t1),b)); % we are solving min problem
    end
    
    bini = (1/nf)*ones(nf,1);
    options = optimoptions('fmincon','Display','off');
    
    [b_reg,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(fun,bini,[],[],ones(1,nf),1,zeros(nf,1),ones(nf,1),[],options);
    mat_simplex_score_b(:,t) = b_reg';
    mat_simplex_score_n(t,1) = sum(b_reg > (1/nf)/100);
    
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
    
    mat_simplex_score_f(t,1) = temp_score;
    
end

% score
res_simplex_score = mean(mat_simplex_score_f(ns0:ns1));
% average number of forecastors
res_simplex_n = mean(mat_simplex_score_n(ns0:ns1));

mat_simplex_b = mat_simplex_score_b';

disp('Simplex is done')

%% Regularized Linear opinion rule - computation

fun = @(lam) SimplexRidgeExPostLambda(xxx, w, ns0, ns1, score_type, lam);
options = optimset('Display','iter');
[lam_opt,fval_reg_L2] = fminbnd(fun, 0.01, 1000, options);

% evaluate at optimal lambda
[mat_reg_L2_score, mat_reg_L2_n, mat_L2_b] = SimplexRidgeExPostLambda_opt(xxx, w, ns0, ns1, score_type, lam_opt);

% score & average number of forecastors
res_reg_L2_score = mean(mat_reg_L2_score(ns0:ns1));
res_reg_L2_n= mean(mat_reg_L2_n(ns0:ns1));

% lambda
res_reg_L2_lambda = lam_opt;

disp('Ridge is done')

%% Dirichlet shrinkage - computation

fun = @(lam) SimplexEntropyExPostLambda(xxx, w, ns0, ns1, score_type, lam);
options = optimset('Display','iter');
[lam_opt,fval_dirichlet] = fminbnd(fun, 0.01, 10, options);

% evaluate at optimal lambda
[mat_dirichlet_score, mat_dirichlet_n, mat_dir_b] = SimplexEntropyExPostLambda_opt(xxx, w, ns0, ns1, score_type, lam_opt);

% score & average number of forecastors
res_dirichlet_score = mean(mat_dirichlet_score(ns0:ns1));
res_dirichlet_n = mean(mat_dirichlet_n(ns0:ns1));

% lambda
res_dirichlet_lambda = lam_opt;

disp('Entropy is done')


%% Averaging methods - Best N-Mixture
% consider all possible combination of N forecastor, evaluate the score for
% each combo and select the best.

% matrix to store historical performance of mixtures
mat_bestNmix_score = nan*ones(ns,Nmax);
mat_bestNmix_set = cell(ns,Nmax);

parfor Nind = 1:Nmax
    
    disp(Nind)
    
    mixset = nchoosek( (1:1:nf)', Nind ); %all possible combinations
    mixset_n = size(mixset,1);
    
    for sind = ns0:ns1
        
        % set-up
        t0 = max(1,sind-w);
        
        % historical performance of all N-mixtures
        temp_score = zeros(mixset_n,(sind-1-t0+1));
        saveind = 1;
        
        for t = t0:1:sind-1
            
            temp_x = xxx(t).histx_fixed;
            temp_y = xxx(t).actual;
            temp_i = temp_y<=temp_x(:,2);
            temp_hist = xxx(t).hist_fixed_nozero;
            
            temp_p = zeros(size(temp_x,1), mixset_n); % # of bins x # of combinations
            for i=1:1:size(temp_x,1)
                temp_p(i,:) = mean(reshape(temp_hist(mixset',i),Nind,[]),1); % calculate the avg prob for each bin and each combination
            end
            
            if score_type == 'R'
                temp_i = temp_y<=temp_x(:,2);
                temp_score(:,saveind)  = sum((cumsum(temp_p,1)/100 - temp_i).^2,1);
                
            elseif score_type == 'B'
                temp_i = zeros(size(temp_x,1),1);
                temp_i((temp_y>temp_x(:,1)) & (temp_y<=temp_x(:,2))) = 1;
                temp_score(:,saveind)  = mean((temp_p/100 - temp_i).^2,1);
                
            elseif score_type == 'L'
                temp_i = zeros(size(temp_x,1),1);
                temp_i((temp_y>temp_x(:,1)) & (temp_y<=temp_x(:,2))) = 1;
                temp_score(:,saveind) = -log( temp_p(temp_i==1,:)/100 );
                
            end
            
            saveind = saveind + 1;
            
        end
        
        temp_score = mean(temp_score,2);
        
        % score (sorted) select the best N-mixture given Nind
        [a,b] = sort(temp_score,'ascend');
        
        % record N-best-mix
        mat_bestNmix_score(sind,Nind) = a(1);
        mat_bestNmix_set{sind,Nind} = mixset(b(1),:);
        
    end
end


% --- OOS (time-t) evaluation
mat_bestNmix_score_eva = nan(ns,Nmax);
mat_bestNmix_score_n   = nan(ns,Nmax);

% sind = ns0
% Nind = 2

for Nind = 1:1:Nmax
    for sind = ns0:1:ns1
        
        temp_hist  = xxx(sind).hist_fixed_nozero;
        temp_histx = xxx(sind).histx_fixed;
        temp_act   = xxx(sind).actual;
        temp_set   = mat_bestNmix_set{sind,Nind};
        
        % evaluate
        best_hist = mean(temp_hist(temp_set,:),1)';
        
        if score_type == 'R'
            temp_i = temp_act<=temp_histx(:,2);
            temp_score = sum((cumsum(best_hist)/100 - temp_i).^2);
            
        elseif score_type == 'B'
            temp_i = zeros(size(temp_histx,1),1);
            temp_i((temp_act>temp_histx(:,1)) & (temp_act<=temp_histx(:,2))) = 1;
            temp_score = mean((best_hist/100 - temp_i).^2);
            
        elseif score_type == 'L'
            temp_i = zeros(size(temp_histx,1),1);
            temp_i((temp_act>temp_histx(:,1)) & (temp_act<=temp_histx(:,2))) = 1;
            temp_score = -log( best_hist(temp_i==1)/100 );
        end
        
        % store
        mat_bestNmix_score_eva(sind,Nind) = temp_score;
        mat_bestNmix_score_n(sind,Nind) = length(temp_set);
    end
end


% average across surveys
res_bestNmix_score = mean(mat_bestNmix_score_eva(ns0:ns1,:),1);
res_bestNmix_n = mean(mat_bestNmix_score_n(ns0:ns1,:),1);


% Averaging methods - best mixture of N<=Nmax
mat_bestmixLessN_score_eva = nan(ns1,Nmax);
mat_bestmixLessN_score_n   = nan(ns1,Nmax);
mat_bestmixLessN_set = cell(ns,Nmax);

% time-t evaluation
for sind = ns0:1:ns1
    
    temp_hist  = xxx(sind).hist_fixed_nozero;
    temp_histx = xxx(sind).histx_fixed;
    temp_act   = xxx(sind).actual;
    
    temp_bestLessN_eva = zeros(Nmax,1);
    temp_bestLessN_n = zeros(Nmax,1);
    
    for i=1:1:Nmax
        
       temp_score = mat_bestNmix_score(sind, 1:i); % for N <= Nmax, fixed W
       [~,b]    = sort(temp_score,'ascend');
       temp_subset = mat_bestNmix_set{sind,b(1)};
       temp_bestLessN_n(i,1) = b(1); % record selected N
       mat_bestmixLessN_set{sind,i} = temp_subset;
        
        % calculate scores
        best_hist = mean(temp_hist(temp_subset,:),1)';
        
        if score_type == 'R'
            temp_i = temp_act<=temp_histx(:,2);
            temp_score = sum((cumsum(best_hist)/100 - temp_i).^2);
            
        elseif score_type == 'B'
            temp_i = zeros(size(temp_histx,1),1);
            temp_i((temp_act>temp_histx(:,1)) & (temp_act<=temp_histx(:,2))) = 1;
            temp_score = mean((best_hist/100 - temp_i).^2);
            
        elseif score_type == 'L'
            temp_i = zeros(size(temp_histx,1),1);
            temp_i((temp_act>temp_histx(:,1)) & (temp_act<=temp_histx(:,2))) = 1;
            temp_score = -log( best_hist(temp_i==1)/100 );
        end
        
        
        temp_bestLessN_eva(i,1) = temp_score;
        
    end
    
    mat_bestmixLessN_score_eva(sind,:) = temp_bestLessN_eva;
    mat_bestmixLessN_score_n(sind,:) = temp_bestLessN_n;
end

% average across survey
res_bestmixLessN_score = mean(mat_bestmixLessN_score_eva(ns0:ns1,:),1);
res_bestmixLessN_n     = mean(mat_bestmixLessN_score_n(ns0:ns1,:),1);

disp('Average of mix is done')


%% individual for comparison
if score_type == 'R'
    mat_score = mat_rps;
    mat_score_avg1 = mat_rps_avg1;
    mat_score_avg2 = mat_rps_avg2;
    
elseif score_type == 'B'
    mat_score = mat_bs;
    mat_score_avg1 = mat_bs_avg1;
    mat_score_avg2 = mat_bs_avg2;
    
elseif score_type == 'L'
    mat_score = mat_ls;
    mat_score_avg1 = mat_ls_avg1;
    mat_score_avg2 = mat_ls_avg2;
end

ind = mean(mat_score(ns0:ns1,:),1);
[a,b] = sort(ind, 'ascend');

res_comparison = quantile(ind, flip([1, 0.9, 0.7, 0.5, 0.3, 0.1, 0]));
disp('Individual is done')

res_simple_average_score1 = mean(mat_score_avg1(ns0:ns1));
res_simple_average_score2 = mean(mat_score_avg2(ns0:ns1));


%% save output
cd(datapath)
save(['current_empirics_', dname, '_', score_type, '_AugVerison_',num2str(version),'_fixedW.mat'])
cd(workpath)










%% --- Table 1: perfornmance of several methods

%--- Part 1: make matlab table

% sub-panel 1
panel1_name   = {'Simplex'; 'Simplex + Ridge'; 'Simplex + Entropy'};
panel1_score  = [res_simplex_score; res_reg_L2_score; res_dirichlet_score];
panel1_N      = [res_simplex_n; res_reg_L2_n; res_dirichlet_n];
panel1_lamdba = {'NA'; res_reg_L2_lambda; res_dirichlet_lambda};

panel1 = table(panel1_name, panel1_score, panel1_N, panel1_lamdba, 'VariableNames',{'title1','R','#','lambda'});

% sub-panel 2
N_list1 = [1:5 10 15 nf];
N_ID_p1 = cell(numel(N_list1),1);
prefix = '\quad $N = ';
for n = 1:numel(N_list1)
    N_ID_p1{n,1} = [prefix num2str(N_list1(n)) '$'];
end

N_list2 = [2:5 nf];
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
panel3_name   = {'Best'; '90\%'; '70\%'; 'Median'; '30\%'; '10\%'; 'Worst'};

panel3_score  = res_comparison';
panel3_N      = ones(numel(panel3_name),1);
panel3_lambda = repmat('NA',length(panel3_name),1);

panel3 = table(panel3_name, panel3_score, panel3_N, panel3_lambda, 'VariableNames',{'title1','R','#','lambda'});


%--- Part 2: generate latex table

fid = fopen([latexpath '/table_emp_' dname '_' score_type '_fixedW.tex'],'wt');

% header and others
fprintf(fid, '%s \n', '\begin{table}[htp]');
fprintf(fid, '%s \n', '\begin{center}');
fprintf(fid, '%s \n', '\caption{Log Scores for Eurozone Inflation Density Forecasts}');
fprintf(fid, '%s \n', '\begin{tabular}{lrr}');
fprintf(fid, '%s \n', '\toprule');fprintf(fid, '%s \n', '\toprule');
fprintf(fid, '%s \n', ['Regularized Mixtures & $ ', score_type ,' $ & \# \\']);
fprintf(fid, '%s \n', '\midrule');

% sub-panel 1: regularization
for n = 1:size(panel1,1)
    fprintf(fid, '%s ', char(panel1{n,1}));
    fprintf(fid, '%s ', '&');
    
    if score_type == 'B'
        fprintf(fid, '%.3f ', panel1{n,2}); % '%.3f ' for B
    else
        fprintf(fid, '%.2f ', panel1{n,2});
    end
    
    fprintf(fid, '%s ', '&');
    fprintf(fid, '%.2f ', panel1{n,3});
    fprintf(fid, '%s \n', '\\');
end
fprintf(fid, '%s \n', '\midrule'); fprintf(fid, '%s \n', '\midrule');

% sub-panel 2: mixture
fprintf(fid, '%s \n', ['Subset Averages & $ ', score_type, ' $ & $\#$ \\']);
fprintf(fid, '%s \n', '\midrule');
fprintf(fid, '%s \n', 'Best $N$-Average: &  & \\');

for n = 1:size(panel2,1)
    fprintf(fid, '%s ', char(panel2{n,1}));
    fprintf(fid, '%s ', '&');
    if score_type == 'B'
        fprintf(fid, '%.3f ', panel2{n,2}); % '%.3f ' for B
    else
        fprintf(fid, '%.2f ', panel2{n,2});
    end
    fprintf(fid, '%s ', '&');
    
    if n > numel(N_list1)
        fprintf(fid, '%.2f ', panel2{n,3});
    else
        fprintf(fid, '%.f ', panel2{n,3});
    end
    
    fprintf(fid, '%s \n', '\\');
end
fprintf(fid, '%s \n', '\midrule'); fprintf(fid, '%s \n', '\midrule');

% sub-panel 3: indiviudal comparisons
fprintf(fid, '%s \n', ['ECB/SPF Comparisons & $ ', score_type, ' $ & $\#$ \\']);
fprintf(fid, '%s \n', '\midrule');
for n = 1:size(panel3,1)
    fprintf(fid, '%s ', char(panel3{n,1}));
    fprintf(fid, '%s ', '&');
    
    if score_type == 'B'
        fprintf(fid, '%.3f ', panel3{n,2}); % '%.3f ' for B
    else
        fprintf(fid, '%.2f ', panel3{n,2});
    end
    
    fprintf(fid, '%s ', '&');
    fprintf(fid, '%.f ', panel3{n,3});
    fprintf(fid, '%s \n', '\\');
end
fprintf(fid, '%s \n', '\midrule');

% sub-panel 4: Simple average
fprintf(fid, '%s ', 'Simple Average');
fprintf(fid, '%s ', '&');
if score_type == 'B'
    fprintf(fid, '%.3f ', res_simple_average_score1); % '%.3f ' for B
else
    fprintf(fid, '%.2f ', res_simple_average_score1);
end
fprintf(fid, '%s ', '&');
fprintf(fid, '%.f ', nf-1);
fprintf(fid, '%s \n', '\\');

fprintf(fid, '%s \n', '\bottomrule'); fprintf(fid, '%s \n', '\bottomrule');

% footer
fprintf(fid, '%s \n', '\end{tabular}');
fprintf(fid, '%s \n', '\end{center}');
fprintf(fid, '%s \n', '\end{table}');
