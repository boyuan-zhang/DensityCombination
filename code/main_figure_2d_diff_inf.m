% Generate 2D heatmaps based on the Empirical analysis of ECB data
% Inflation

% 11/2/2020
% Minchul Shin and Boyuan Zhang

clc; clear all; close all;

% benchmark
% 1y
% log-score


%%  housekeeping
workpath00 = pwd;
savepath00 = [workpath00, filesep, 'graphics'];

datapath = '../data/';
addpath(datapath);

% -------------------------------------
% USER HAS TO CHANGE THIS
vname = 'infl_1y';
dname = 'ecbspf_infl_1y_bp_nozero';
% -------------------------------------

addpath(genpath('toolbox_subfunc'));
addpath(genpath('toolbox_plot'));

chk_dir(savepath00);

%% load data
load([dname, '.mat']);
eval(['xxx = ', dname, ';']);


%% Set-up
ns0 = 2; %2001Q1
ns1 = 83; %2019Q3
ns = ns1;

%% Score measures
[nf, ~] = size(xxx(1).hist);
mat_rps = nan(ns, nf);
mat_bs  = nan(ns, nf);
mat_ls  = nan(ns, nf);

mat_rps_avg = nan(ns, 1);
mat_bs_avg  = nan(ns, 1);
mat_ls_avg  = nan(ns, 1);

for sind = 1:1:ns1
    
    % individual rps, bs, ls
    for i=1:nf
        mat_rps(sind,i) = rps(xxx(sind).hist_fixed_nozero(i,:)', xxx(sind).histx_fixed, xxx(sind).actual);
        mat_bs(sind,i) = bs(xxx(sind).hist_fixed_nozero(i,:)', xxx(sind).histx_fixed, xxx(sind).actual);
        mat_ls(sind,i) = ls(xxx(sind).hist_fixed_nozero(i,:)', xxx(sind).histx_fixed, xxx(sind).actual);
    end
    
    % average of selected (after extrapolation)
    temp_p = mean(xxx(sind).hist_fixed_nozero, 1)';
    temp_x = xxx(sind).histx_fixed(:,:);
    temp_y = xxx(sind).actual;
    
    mat_rps_avg(sind,:) = rps(temp_p, temp_x, temp_y);
    mat_bs_avg(sind,:) = bs(temp_p, temp_x, temp_y);
    mat_ls_avg(sind,:) = ls(temp_p, temp_x, temp_y);
    
end

%% Load Regularization-based methods
load('current_empirics_infl_1y_L_AugVerison_2_fixedW.mat')

%% get b for the subset averaging
% -> we draw figure based on Best <=4-Average
bestNind = 4;

v_ones = ones(1,nf);

mat_subset_lessN_b = zeros(ns1,nf);

for sind = ns0:1:ns1
    temp_indx = mat_bestmixLessN_set{sind,bestNind};

    temp_b = zeros(1,nf);
    temp_b(temp_indx) = 1/numel(temp_indx);
    
    mat_subset_lessN_b(sind,:) = temp_b;
end


%% Heatmap

% ---------------------------
% Figure 3, 
% simple average
% (without hypothetical forecaster)

fig = figure(11);
setmyfig(fig, [2,2,12,6]);
mat_b = 1/(nf-1)*ones(size(mat_L2_b,1), nf-1);

temp_x = xxx(1).histx_fixed;
temp_x(1,1) = temp_x(1,2)-0.5;
temp_x(end,2) = temp_x(end,1)+0.5;
temp_mid = temp_x(:,1) + (temp_x(:,2)-temp_x(:,1))/2;

big_grid = temp_mid;
ngrid = numel(big_grid);
p_grid = zeros(ns,ngrid);
m_grid = zeros(ns,1);
s_grid = zeros(ns,1);

for sind = ns0:1:ns1
    
    % Simple average
    temp_x = xxx(sind).histx_fixed;
    temp_x(1,1) = temp_x(1,2)-0.5;
    temp_x(end,2) = temp_x(end,1)+0.5;
    temp_mid = temp_x(:,1) + (temp_x(:,2)-temp_x(:,1))/2;
    temp_hist0 = xxx(sind).hist_fixed_nozero;
    temp_hist = mean(temp_hist0,1);
    
    % probability
    p_grid(sind,:) = mat_b(sind,:) * xxx(sind).hist_fixed_nozero(1:end-1,:);
    
    % mean and sd
    m_grid(sind,1) = (temp_hist/100)*temp_mid;
    s_grid(sind,1) = sqrt( ((temp_hist/100)*temp_mid.^2)-m_grid(sind,1).^2);    
end

% actual figure
% mymap = cbrewer('seq', 'Greys', 150);
mymap = (cbrewer('seq', 'YlOrBr', 100));
% mymap = mymap(1:125,:);

imagesc(ns0:ns1, big_grid, p_grid(ns0:ns1,:)')

colormap(mymap)
cb = colorbar('eastoutside');
caxis([0,50]) ; %change color axis
% xlabel('Time');

cb.Label.String = '(%)';
cb.Ticks = (0:5:50);

xticks_id = 1:8:ns;
xticks_date = {};

for i = 1:numel(xticks_id)
   xticks_date{i} = xxx(xticks_id(i)).sdate;
end

set(gca,'YDir','normal')

xticks(xticks_id)
xticklabels(xticks_date)
yticks(0:4)
% set(gca,'TickLength',[0.1, 0.01], 'fontsize', 15)
set(gca,'fontsize', 18)

% view(-70,60) % fix the angles of the view

% xlabel('Survey date');
% ylabel('Inflation rate (%)')
ylabel('Forecast and Realized Inflation (%)')

% title('Simplex Minus Simple Average','FontWeight','Normal', 'fontsize', 20)

hold on

temp_act = [];
for i = 1:ns
   temp_act = [temp_act; xxx(i).actual];
end

plot(1:ns, temp_act, 'linewidth', 4, 'MarkerSize', 10, 'color', 'black'); % rgb('DimGray')
yline(0, 'black--', 'LineWidth', 1.5)
hold off


cd(savepath00)
saveas(fig, 'fig_heatmap_simple_avg.png');
cd(workpath00)

close all


% ---------------------------
% Figure 6
% subset (<=N) - simplex average
% (with hypothetical forecaster)

fig = figure(133);
setmyfig(fig, [2,2,12,5.75]);
mat_b1 = mat_subset_lessN_b;
mat_b2 = 1/(nf-1)*ones(size(mat_L2_b,1), nf-1);
script_heatmap_plot_diff_uniform;
% title('Best \leq 4-Average Minus Simple Average','FontWeight','Normal', 'fontsize', 20)

hold on

temp_act = [];
for i = 1:ns
   temp_act = [temp_act; xxx(i).actual];
end

plot(1:ns, temp_act, 'linewidth', 4, 'MarkerSize', 10, 'color', 'black'); % rgb('DimGray')
yline(0, 'black--', 'LineWidth', 1.5)
hold off

cd(savepath00)
saveas(fig, 'fig_heatmap_diff_subsetleN_avg_uni.png');
cd(workpath00)

close all

