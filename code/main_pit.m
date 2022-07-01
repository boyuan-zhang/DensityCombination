% PIT analysis of ECB data

% 11/17/2020
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

%% Figure 7
fig = figure(1);

%% (1,1) PIT for simple average (18 forecasters)
subplot(2,2,1)

ns0 = 1;
ns1 = 36; %2007Q4
ns = ns1;

% randomized version of PIT 
nsim = 500;
mat_pit_random = nan(ns, nsim);

for sind = ns0:1:ns1
    
    % average of selected (after extrapolation)
    temp_p = mean(xxx(sind).hist_fixed_nozero(1:end-1,:), 1)';
    temp_x = xxx(sind).histx_fixed(:,:);
    temp_y = xxx(sind).actual;
    
    % Find bin position
    temp_y_bin_pos = temp_x(:,1) < temp_y & temp_x(:,2) >= temp_y;
    
    % compute randomized version of PIT (eqn(1) of Czado Gneiting, Held (2009) )
    temp_P1 = cumsum(temp_p);
    temp_P0 = [0; temp_P1(1:end-1)];
    pp0 = temp_P0(temp_y_bin_pos)/100;
    pp1 = temp_P1(temp_y_bin_pos)/100;
    temp_pit = pp0 + ( rand(1,nsim) * (pp1 - pp0) );
    
    % store
    mat_pit_random(sind,:) = temp_pit;
end

% PIT, uniform 
script_pit_uniform
text(0.05,0.45,'Simple Average','FontSize',24)

% title('PIT Histogram');
% t = title({'First','Subsample'},'FontSize',30, 'fontweight','normal');
t = title('2001Q1-2007Q4','FontSize',30, 'fontweight','normal');
set(t,'Position',get(t,'Position')+ [0 .01 0]);  % move up slightly
% title('Simple Average', 'Units', 'normalized', 'Position', [0.25, 0.85, 0], 'fontweight','normal')

%% (1,2) PIT for simple average (18 forecasters)
subplot(2,2,2)

ns0 = 37; %2008Q1
ns1 = 83; %2019Q3
ns = ns1;

% randomized version of PIT 
nsim = 500;
mat_pit_random = nan(ns, nsim);

for sind = ns0:1:ns1
    
    % average of selected (after extrapolation)
    temp_p = mean(xxx(sind).hist_fixed_nozero(1:end-1,:), 1)';
    temp_x = xxx(sind).histx_fixed(:,:);
    temp_y = xxx(sind).actual;
    
    % Find bin position
    temp_y_bin_pos = temp_x(:,1) < temp_y & temp_x(:,2) >= temp_y;
    
    % compute randomized version of PIT (eqn(1) of Czado Gneiting, Held (2009) )
    temp_P1 = cumsum(temp_p);
    temp_P0 = [0; temp_P1(1:end-1)];
    pp0 = temp_P0(temp_y_bin_pos)/100;
    pp1 = temp_P1(temp_y_bin_pos)/100;
    temp_pit = pp0 + ( rand(1,nsim) * (pp1 - pp0) );
    
    % store
    mat_pit_random(sind,:) = temp_pit;
end

% PIT, uniform 
script_pit_uniform
text(0.05,0.45,'Simple Average','FontSize',24)

% t = title({'Second','Subsample'},'FontSize',30, 'fontweight','normal');
t = title('2008Q1-2019Q3','FontSize',30, 'fontweight','normal');
set(t,'Position',get(t,'Position')+ [0 .01 0]);  % move up slightly

%% (2,1) PIT for subset <= 4
subplot(2,2,3)

ns0 = 9; %2001Q1
ns1 = 36; %2007Q4
ns = ns1;

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

% randomized version of PIT 
nsim = 500;
mat_pit_random = nan(ns, nsim);

for sind = ns0:1:ns1
    
    % average of selected (after extrapolation)
    temp_p = ( mat_subset_lessN_b(sind,:) * xxx(sind).hist_fixed_nozero )';
%     temp_p = xxx(sind).hist_fixed_nozero(12,:)'; %individual
    temp_x = xxx(sind).histx_fixed(:,:);
    temp_y = xxx(sind).actual;
    
    % Find bin position
    temp_y_bin_pos = temp_x(:,1) < temp_y & temp_x(:,2) >= temp_y;
    
    % compute randomized version of PIT (eqn(1) of Czado Gneiting, Held (2009) )
    temp_P1 = cumsum(temp_p);
    temp_P0 = [0; temp_P1(1:end-1)];
    pp0 = temp_P0(temp_y_bin_pos)/100;
    pp1 = temp_P1(temp_y_bin_pos)/100;
    temp_pit = pp0 + ( rand(1,nsim) * (pp1 - pp0) );
    
    % store
    mat_pit_random(sind,:) = temp_pit;
end

% PIT, uniform 
script_pit_uniform

% title('PIT Histogram');
% title('Best \leq 4-Average', 'Units', 'normalized', 'Position', [0.29, 0.85, 0], 'fontweight','normal')
text(0.05,0.45,'Best \leq4-Average','FontSize',24)

%% (2,2) PIT for subset <= 4
subplot(2,2,4)

ns0 = 37; %2008Q1
ns1 = 83; %2019Q3
ns = ns1;

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


% randomized version of PIT 
nsim = 500;
mat_pit_random = nan(ns, nsim);

for sind = ns0:1:ns1
    
    % average of selected (after extrapolation)
    temp_p = ( mat_subset_lessN_b(sind,:) * xxx(sind).hist_fixed_nozero )';
%     temp_p = xxx(sind).hist_fixed_nozero(12,:)'; %individual
    temp_x = xxx(sind).histx_fixed(:,:);
    temp_y = xxx(sind).actual;
    
    % Find bin position
    temp_y_bin_pos = temp_x(:,1) < temp_y & temp_x(:,2) >= temp_y;
    
    % compute randomized version of PIT (eqn(1) of Czado Gneiting, Held (2009) )
    temp_P1 = cumsum(temp_p);
    temp_P0 = [0; temp_P1(1:end-1)];
    pp0 = temp_P0(temp_y_bin_pos)/100;
    pp1 = temp_P1(temp_y_bin_pos)/100;
    temp_pit = pp0 + ( rand(1,nsim) * (pp1 - pp0) );
    
    % store
    mat_pit_random(sind,:) = temp_pit;
end

% PIT, uniform 
script_pit_uniform
text(0.05,0.45,'Best \leq4-Average','FontSize',24)


%%

% 
ha=get(gcf,'children');

set(ha(1),'position',[.55 .05 .4125 .375])
set(ha(2),'position',[.05 .05 .4125 .375])
set(ha(3),'position',[.55 .52 .4125 .375])
set(ha(4),'position',[.05 .52 .4125 .375])

setmyfig(fig, [1,1,16,12]);

cd(savepath00)
saveas(fig, 'fig_pit_uniform_4panel_v2.png');
cd(workpath00)
