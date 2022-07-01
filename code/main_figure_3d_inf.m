% Generate 3D plots based on the Empirical analysis of ECB data
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


%% (1a) Figure 5: Distribution over time, inflation, 3-D graph, average

temp_x = xxx(1).histx_fixed;
temp_x(1,1) = temp_x(1,2)-0.5;
temp_x(end,2) = temp_x(end,1)+0.5;
temp_mid = temp_x(:,1) + (temp_x(:,2)-temp_x(:,1))/2;

big_grid = temp_mid;
ngrid = numel(big_grid);
p_grid = zeros(ns,ngrid);

for sind = 1:1:ns
    
    temp_x = xxx(sind).histx_fixed;
    temp_x(1,1) = temp_x(1,2)-0.5;
    temp_x(end,2) = temp_x(end,1)+0.5;
    temp_mid = temp_x(:,1) + (temp_x(:,2)-temp_x(:,1))/2;
    
    temp_hist = xxx(sind).hist_fixed_nozero;
    temp_hist = mean(temp_hist(1:end-1,:),1);
    
    % probability
    for j=1:1:numel(temp_mid)
        loc_j = find(big_grid==temp_mid(j));
        p_grid(sind,loc_j) = temp_hist(j);
    end
    
end

%
fig = figure(19);
setmyfig(fig, [2,2,12,8]);

mymap = (cbrewer('seq', 'YlOrBr', 100));
s = surf(1:ns, big_grid, p_grid(:,:)','FaceAlpha',1);
colormap(mymap)

xticks_id = 11:10:ns;
xticks_date = {};

for i = 1:numel(xticks_id)
   xticks_date{i} = xxx(xticks_id(i)).sdate;
end

xticks(xticks_id)
xticklabels(xticks_date)
xlim([-1 ns+1])
zlim([0 60])
set(gca,'TickLength',[0.1, 0.01], 'fontsize', 15, 'Ydir','reverse')
 
view(-70,60) % fix the angles of the view
grid on

xlabel('Survey date');
ylabel('Inflation rate (%)')
zlabel('Probability');

% title('Simple Average','FontWeight','Normal', 'fontsize', 23, 'Position', [8, 1.65, 120])
title('Simple Average','FontWeight','Normal', 'fontsize', 25, 'Position', [8, 1.67, 125])

cd(savepath00);
saveas(fig,['fig_',vname,'_Density_Average','.png']);
cd(workpath00);

close all

%% (1b) Figure 3: Average distribution over time, inflation, 3-D graph, no text
%
fig = figure(19);
setmyfig(fig, [2,2,12,8]);

mymap = (cbrewer('seq', 'YlOrBr', 100));
s = surf(1:ns, big_grid, p_grid(:,:)','FaceAlpha',1);
colormap(mymap)

xticks_id = 11:10:ns;
xticks_date = {};

for i = 1:numel(xticks_id)
   xticks_date{i} = xxx(xticks_id(i)).sdate;
end

xticks(xticks_id)
xticklabels(xticks_date)
xlim([-1 ns+1])
zlim([0 60])
set(gca,'TickLength',[0.1, 0.01], 'fontsize', 19, 'Ydir','reverse')
 
view(-70,60) % fix the angles of the view
grid on

xlabel('Survey date');
ylabel('Inflation rate (%)')
zlabel('Probability');

cd(savepath00);
saveas(fig,['fig_',vname,'_Density_Average_less_text','.png']);
cd(workpath00);

close all

%% (2) Figure 4: Distribution over time, inflation, 3-D graph, simplex
mat_b = mat_simplex_b;
ns0 = 9;

temp_x = xxx(1).histx_fixed;
temp_x(1,1) = temp_x(1,2)-0.5;
temp_x(end,2) = temp_x(end,1)+0.5;
temp_mid = temp_x(:,1) + (temp_x(:,2)-temp_x(:,1))/2;

big_grid = temp_mid;
ngrid = numel(big_grid);
p_grid = zeros(ns,ngrid);

for sind = 1:1:ns
    
    temp_x = xxx(sind).histx_fixed;
    temp_x(1,1) = temp_x(1,2)-0.5;
    temp_x(end,2) = temp_x(end,1)+0.5;
    temp_mid = temp_x(:,1) + (temp_x(:,2)-temp_x(:,1))/2;
    
     temp_hist = mat_b(sind,:) * xxx(sind).hist_fixed_nozero;
    
    % probability
    for j=1:1:numel(temp_mid)
        loc_j = find(big_grid==temp_mid(j));
        p_grid(sind,loc_j) = temp_hist(j);
    end
    
end

% create map
fig = figure(15);
setmyfig(fig, [2,2,12,8]);

mymap = (cbrewer('seq', 'YlOrBr', 100));
surf(ns0:ns, big_grid, p_grid(ns0:ns,:)','FaceAlpha',1);
colormap(mymap)

xticks_id = 11:10:ns;
xticks_date = {};

for i = 1:numel(xticks_id)
   xticks_date{i} = xxx(xticks_id(i)).sdate;
end

xticks(xticks_id)
xticklabels(xticks_date)
xlim([ns0-1 ns+1])
zlim([0 60])
set(gca,'TickLength',[0.1, 0.01], 'fontsize', 15, 'Ydir','reverse')

view(-70,60) % fix the angles of the view

xlabel('Survey date');
ylabel('Inflation rate (%)')
zlabel('Probability');

grid on
title('Simplex','FontWeight','Normal', 'fontsize', 25, 'Position', [15, 2.4, 138])


cd(savepath00);
saveas(fig,['fig_',vname,'_Density_Simplex','.png']);
cd(workpath00);

close all

%% (3) Figure 4 & 5: Distribution over time, inflation, 3-D graph, SubsetleN
mat_b = mat_subset_lessN_b;
ns0 = 9;

temp_x = xxx(1).histx_fixed;
temp_x(1,1) = temp_x(1,2)-0.5;
temp_x(end,2) = temp_x(end,1)+0.5;
temp_mid = temp_x(:,1) + (temp_x(:,2)-temp_x(:,1))/2;

big_grid = temp_mid;
ngrid = numel(big_grid);
p_grid = zeros(ns,ngrid);

for sind = 1:1:ns
    
    temp_x = xxx(sind).histx_fixed;
    temp_x(1,1) = temp_x(1,2)-0.5;
    temp_x(end,2) = temp_x(end,1)+0.5;
    temp_mid = temp_x(:,1) + (temp_x(:,2)-temp_x(:,1))/2;
    
     temp_hist = mat_b(sind,:) * xxx(sind).hist_fixed_nozero;
    
    % probability
    for j=1:1:numel(temp_mid)
        loc_j = find(big_grid==temp_mid(j));
        p_grid(sind,loc_j) = temp_hist(j);
    end
    
end


% create map
fig = figure(17);
setmyfig(fig, [2,2,12,8]);

mymap = (cbrewer('seq', 'YlOrBr', 100));
surf(ns0:ns, big_grid, p_grid(ns0:ns,:)','FaceAlpha',1);
colormap(mymap)

xticks_id = 11:10:ns;
xticks_date = {};

for i = 1:numel(xticks_id)
   xticks_date{i} = xxx(xticks_id(i)).sdate;
end

xticks(xticks_id)
xticklabels(xticks_date)
xlim([ns0-1 ns+1])
zlim([0 60])
set(gca,'TickLength',[0.1, 0.01], 'fontsize', 15, 'Ydir','reverse')

view(-70,60) % fix the angles of the view

xlabel('Survey date');
ylabel('Inflation rate (%)')
zlabel('Probability');

grid on
title('Best \leq4-Average','FontWeight','Normal', 'fontsize', 25, 'Position', [15, 1.8, 133])

cd(savepath00);
saveas(fig,['fig_',vname,'_Density_SubsetleN','.png']);
cd(workpath00);

close all
