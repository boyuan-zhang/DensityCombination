% Generate 3D plots based on the Empirical analysis of ECB data
% Real interest rate

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
load('data_nominal_rate.mat');

load([dname, '.mat']);
eval(['xxx = ', dname, ';']);

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


%% (1) Figure 8: Distribution over time, inflation, 3-D graph, average

temp_grid = [];
for sind=1:1:ns
    
    temp_x = xxx(sind).histx_fixed - mat_ni_data(sind);
    temp_x(1,1) = temp_x(1,2)-0.5;
    temp_x(end,2) = temp_x(end,1)+0.5;
    
    temp_x = -temp_x;
    
    temp_mid = temp_x(:,1) + (temp_x(:,2)-temp_x(:,1))/2;

    temp_grid = [temp_grid; temp_mid];
end


big_grid = unique(temp_grid);
ngrid = numel(big_grid);
p_grid = zeros(ns,ngrid);

for sind = 1:1:ns
    
    % grid x
    temp_x = xxx(sind).histx_fixed - mat_ni_data(sind);
    temp_x(1,1) = temp_x(1,2)-0.5;
    temp_x(end,2) = temp_x(end,1)+0.5;
    temp_x = -temp_x;
    temp_mid = temp_x(:,1) + (temp_x(:,2)-temp_x(:,1))/2;
    
    temp_hist = xxx(sind).hist_fixed_nozero;
    temp_hist = mean(temp_hist(1:end-1,:),1);
    
    % avereage probability
    for j=1:1:numel(temp_mid)
        loc_j = find(big_grid==temp_mid(j));
        p_grid(sind,loc_j) = temp_hist(j);
    end
 
end

% Adjust big_grid and p_grid
grid_max = round(max(big_grid)*4)/4;
grid_min = round(min(big_grid)*4)/4;

if max(big_grid) > grid_max || min(big_grid) < grid_min
   grid_max = grid_max + 0.25;
   grid_min = grid_min + 0.25;
end

% big_grid_new = grid_min:0.5:grid_max;
big_grid_new = grid_min:0.5:grid_max;
ngrid_new = numel(big_grid_new);
p_grid_new = zeros(ns, ngrid_new);

% assign each midpoint in big_grid to the nearest point in the new and small grid 
grid_convert_id = zeros(ngrid,1);
for i = 1:1:ngrid
   temp_dist = abs( big_grid(i) - big_grid_new );
   [a, b] = min(temp_dist);
   grid_convert_id(i) = b;
end

% aggregate prob
for sind = 1:1:ns
   for i = 1:1:ngrid_new
      temp_hist = p_grid(sind, grid_convert_id == i);
      p_grid_new(sind, i) = sum(temp_hist);
   end
end


% create map
fig = figure(1);
setmyfig(fig, [2,2,12,8]);

mymap = (cbrewer('seq', 'YlOrBr', 100));
surf(1:ns, big_grid_new, p_grid_new(:,:)','FaceAlpha',1);
colormap(mymap)

xticks_id = 11:10:ns;
xticks_date = {};

for i = 1:numel(xticks_id)
   xticks_date{i} = xxx(xticks_id(i)).sdate;
end

xticks(xticks_id)
xticklabels(xticks_date)
xlim([0 ns+1])
zlim([0 65])
set(gca,'TickLength',[0.1, 0.01], 'fontsize', 15, 'Ydir','reverse')

view(-70,60) % fix the angles of the view

xlabel('Survey date');
ylh = ylabel('Real interest rate (%)');
zlabel('Probability');

ylh.Position(2) = ylh.Position(2) + 0.6; 
ylh.Position(1) = ylh.Position(1) + 15; 

grid on
title('Simple Average','FontWeight','Normal', 'fontsize', 25, 'Position', [30, 2.5, 115])

cd(savepath00);
saveas(fig,['fig_real_',vname,'_Density_Average','.png']);
cd(workpath00);


%% (2) Figure 8: Distribution over time, real interest rate, 3-D graph, subset
mat_b = mat_subset_lessN_b;
ns0 = 9;

temp_grid = [];

for sind=ns0:1:ns
    
    temp_x = xxx(sind).histx_fixed - mat_ni_data(sind);
    temp_x(1,1) = temp_x(1,2)-0.5;
    temp_x(end,2) = temp_x(end,1)+0.5;
    
    temp_x = -temp_x;
    
    temp_mid = temp_x(:,1) + (temp_x(:,2)-temp_x(:,1))/2;

    temp_grid = [temp_grid; temp_mid];
end


big_grid = unique(temp_grid);
ngrid = numel(big_grid);
p_grid = zeros(ns,ngrid);

for sind = ns0:1:ns
    
    % grid x
    temp_x = xxx(sind).histx_fixed - mat_ni_data(sind);
    temp_x(1,1) = temp_x(1,2)-0.5;
    temp_x(end,2) = temp_x(end,1)+0.5;
    temp_x = -temp_x;
    temp_mid = temp_x(:,1) + (temp_x(:,2)-temp_x(:,1))/2;
    
    temp_hist = mat_b(sind,:) * xxx(sind).hist_fixed_nozero;
    
    % avereage probability
    for j=1:1:numel(temp_mid)
        loc_j = find(big_grid==temp_mid(j));
        p_grid(sind,loc_j) = temp_hist(j);
    end
 
end


% Adjust big_grid and p_grid
grid_max = round(max(big_grid)*4)/4;
grid_min = round(min(big_grid)*4)/4;

if max(big_grid) > grid_max || min(big_grid) < grid_min
   grid_max = grid_max + 0.25;
   grid_min = grid_min + 0.25;
end

% big_grid_new = grid_min:0.5:grid_max;
big_grid_new = grid_min:0.5:grid_max;
ngrid_new = numel(big_grid_new);
p_grid_new = zeros(ns, ngrid_new);

% assign each midpoint in big_grid to the nearest point in the new and small grid 
grid_convert_id = zeros(ngrid,1);
for i = 1:1:ngrid
   temp_dist = abs( big_grid(i) - big_grid_new );
   [a, b] = min(temp_dist);
   grid_convert_id(i) = b;
end

% aggregate prob
for sind = ns0:1:ns
   for i = 1:1:ngrid_new
      temp_hist = p_grid(sind, grid_convert_id == i);
      p_grid_new(sind, i) = sum(temp_hist);
   end
end


% create map
fig = figure(6);
setmyfig(fig, [2,2,12,8]);

mymap = (cbrewer('seq', 'YlOrBr', 100));
surf(ns0:ns, big_grid_new, p_grid_new(ns0:ns,:)','FaceAlpha',1);
colormap(mymap)

xticks_id = 11:10:ns;
xticks_date = {};

for i = 1:numel(xticks_id)
   xticks_date{i} = xxx(xticks_id(i)).sdate;
end

xticks(xticks_id)
xticklabels(xticks_date)
xlim([ns0-1 ns+1])
zlim([0 65])
set(gca,'TickLength',[0.1, 0.01], 'fontsize', 15, 'Ydir','reverse')

view(-70,60) % fix the angles of the view

xlabel('Survey date');
ylh = ylabel('Real interest rate (%)');
zlabel('Probability');

ylh.Position(2) = ylh.Position(2) + 0.5; 
ylh.Position(1) = ylh.Position(1) + 14; 

grid on
title('Best \leq4-Average','FontWeight','Normal', 'fontsize', 25, 'Position', [15, 3.5, 155])

cd(savepath00);
saveas(fig,['fig_real_',vname,'_Density_SubsetleN','.png']);
cd(workpath00);

close all
