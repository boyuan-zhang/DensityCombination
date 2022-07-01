% Descriptive analysis of ECB data 
% date: 2021/8/9
% Minchul Shin and Boyuan Zhang

% Individual and Average Density Forecasts at some point of time

% -------------------------------------
% _fixed: fixed bins with nozero bins
% _2yc: next year generated in every Q4
% -------------------------------------

clc; clear all; close all;

% -------------------------------------
% USER HAS TO CHANGE THIS
% vname = 'infl_2yc';
% dname = 'ecbspf_infl_2yc_bp_nozero';
vname = 'infl_1y';
dname = 'ecbspf_infl_1y_bp_nozero';
% -------------------------------------

workpath = pwd;
datapath = '../data/';
addpath(datapath);

savepath = [workpath, filesep, 'graphics'];

addpath(genpath('toolbox_subfunc'));
addpath(genpath('toolbox_plot'));

chk_dir(savepath)

%% load data

load('data_nominal_rate.mat');

load([dname, '.mat']);
eval(['xxx = ', dname, ';']);


%% Cut to have every 4Q (only for 2yc)

if strcmp(vname, 'infl_2yc')
   xxx_orig = xxx;
   xxx = xxx(4:4:84);
   
   cell_ni_data = cell_ni_data(4:4:84);
   mat_ni_data = mat_ni_data(4:4:84);
end


%% (1) inflation at 2018Q4

sind = 80;
temp_x = xxx(sind).histx_fixed;
temp_x(1,1) = temp_x(1,2)-0.5;
temp_x(end,2) = temp_x(end,1)+0.5;
temp_mid = temp_x(:,1) + (temp_x(:,2)-temp_x(:,1))/2;
temp_hist = xxx(sind).hist_fixed_nozero;
temp_hist_avg = mean(temp_hist,1);
temp_act = xxx(sind).actual;

fig = figure(1);
setmyfig(fig, [2,2,5,4]);

plot(temp_mid, temp_hist(:,:),'-*', 'color', rgb('grey'), 'linewidth', 1)
xlim([-1,4.5]);
ylim([0,80]);
title(xxx(sind).sdate,'FontWeight','normal')
set(gca,'linewidth', 1, 'fontsize', 15);
grid on
hold on
bar(temp_mid, temp_hist_avg, 'FaceAlpha', 0.5, 'FaceColor', rgb('orange'));
xline(temp_act,'LineWidth', 2.2,'color','black')

hold off
xlabel('Inflation rate (%)');
ylabel('Probability');

cd(savepath)
saveas(fig, ['fig_','avg_',vname,'_',xxx(sind).sdate, '.png']);
cd(workpath)

close all


%% (2) inflation at 2004Q4

sind = 24;
temp_x = xxx(sind).histx_fixed;
temp_x(1,1) = temp_x(1,2)-0.5;
temp_x(end,2) = temp_x(end,1)+0.5;
temp_mid = temp_x(:,1) + (temp_x(:,2)-temp_x(:,1))/2;
temp_hist = xxx(sind).hist_fixed_nozero;
temp_hist_avg = mean(temp_hist,1);
temp_act = xxx(sind).actual;


fig = figure(1);
setmyfig(fig, [2,2,5,4]);

plot(temp_mid, temp_hist(:,:),'-*', 'color', rgb('grey'), 'linewidth', 1)

xlim([-1,4.5]);
ylim([0,80]);

% title(['Survey from = ', xxx(sind).sdate, ', all versus average']);
title(xxx(sind).sdate,'FontWeight','normal')
set(gca,'linewidth', 1, 'fontsize', 15);
grid on
hold on
bar(temp_mid, temp_hist_avg, 'FaceAlpha', 0.5, 'FaceColor', rgb('orange'));
xline(temp_act,'LineWidth', 2.2,'color','black')
hold off
xlabel('Inflation rate (%)');
ylabel('Probability');


cd(savepath)
saveas(fig, ['fig_','avg_',vname,'_',xxx(sind).sdate, '.png']);
cd(workpath)


close all

