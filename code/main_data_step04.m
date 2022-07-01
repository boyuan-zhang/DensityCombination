% Check SPF data, bins
% ECB - inflation 

% Illustration of 
% - Time-varying bins and the correction
% - Zero probability issue

clc; clear all; close all;
workpath = pwd;
datapath = '../data/';

addpath(datapath);
addpath(genpath('toolbox_plot'))

data_type = 1; %1=1y, 2=2yc

if data_type == 1
    xxx = load('ecbspf_infl_1y_bp.mat');
    xxx = xxx.ecbspf_infl_1y_bp;
    varname = 'infl_1y';
elseif data_type == 2
    xxx = load('ecbspf_infl_2yc_bp.mat');
    xxx = xxx.ecbspf_infl_2yc_bp;
    varname = 'infl_2yc';
end

%% Histogram bins, inflation rate

fig = figure(1);
setmyfig(fig, [2,2,10,4]);

t = 1;
plot(t,0);
hold on

temp_xlab = {};
for t=1:size(xxx,1)
    temp_l = xxx(t).histx(:,1);
    temp_r = xxx(t).histx(:,2);
    plot(t, unique([temp_l; temp_r]), 'p', 'linewidth', 1, 'color','k')
    temp_xlab = [temp_xlab; xxx(t).sdate];
end
hold off
temp_xtic = (1:8:84)';
set(gca, 'Xtick', temp_xtic, 'XtickLabel', temp_xlab(temp_xtic));
xlim([-1, 86]);
ylim([-3,5]);
set(gca, 'fontsize', 15, 'linewidth', 2);
grid on


%% Merging histograms 

% fixed bins
fixed_histx = xxx(1).histx(:,:);
fixed_n = size(fixed_histx,1);

% loop to update
for t = 1:1:size(xxx,1)
    
    temp_prob = xxx(t).hist;
    temp_nf = size(temp_prob,1);
    index_low = xxx(t).histx(:,2) <= fixed_histx(1,2);
    index_up = xxx(t).histx(:,1) >= fixed_histx(end,1);
    
    % merge prob
    xxx(t).hist_fixed = nan(temp_nf, fixed_n);
    xxx(t).hist_fixed(:,1) = sum(temp_prob(:,index_low),2);
    xxx(t).hist_fixed(:,end) = sum(temp_prob(:,index_up),2);
    xxx(t).hist_fixed(:,2:end-1) = temp_prob(:,~index_low & ~index_up);
    
    % fixed bin
    xxx(t).histx_fixed = fixed_histx;
end


%% Histogram bins, inflation rate, corrected bins

fig = figure(2);
setmyfig(fig, [2,2,10,4]);

t = 1;
plot(t,0);
hold on

temp_xlab = {};
for t=1:size(xxx,1)
    temp_l = xxx(t).histx_fixed(:,1);
    temp_r = xxx(t).histx_fixed(:,2);
    plot(t, unique([temp_l; temp_r]), 'p', 'linewidth', 1, 'color','k')
     
    temp_xlab = [temp_xlab; xxx(t).sdate];
end
hold off
temp_xtic = (1:8:84)';
set(gca, 'Xtick', temp_xtic, 'XtickLabel', temp_xlab(temp_xtic));
xlim([-1, 86]);
ylim([-3,5]);
set(gca, 'fontsize', 15, 'linewidth', 2);
grid on


%% Zero probability

tab_prob = [];
temp_xlab = {};
for t = 1:1:size(xxx,1)-1
    temp_index_in = (xxx(t).histx_fixed(:,1) < xxx(t).actual) & (xxx(t).histx_fixed(:,2) >= xxx(t).actual);
    temp_prob = xxx(t).hist_fixed(:,temp_index_in);
    tab_prob = [tab_prob; temp_prob'];
    temp_xlab = [temp_xlab; xxx(t).sdate];
end


fig = figure(3);
setmyfig(fig, [2,2,10,4]);

nf = size(tab_prob,2);
plot(sum(tab_prob==0,2), '*-', 'linewidth', 2)
hold on
plot(-1:86,nf*ones(88,1), '--', 'linewidth', 3);
hold off

temp_xtic = (1:16:84)';
set(gca, 'Xtick', temp_xtic, 'XtickLabel', temp_xlab(temp_xtic));
xlim([-1, 86]);
ylim([-1,nf+5]);
set(gca, 'fontsize', 15, 'linewidth', 2);
grid on
ylabel('# of forecasters');


%% zero probability correction
nf = size(xxx(1).hist,1);
tab_prob = [];
temp_xlab = {};

for t = 1:1:size(xxx,1)-1
    temp_index_in = (xxx(t).histx_fixed(:,1) < xxx(t).actual) & (xxx(t).histx_fixed(:,2) >= xxx(t).actual);
    
%     temp_prob = xxx(t).hist_fixed(:,temp_index_in);
%     tab_prob = [tab_prob; temp_prob'];
%     temp_xlab = [temp_xlab; xxx(t).sdate];
    
    xxx(t).hist_fixed_nozero = xxx(t).hist_fixed;
    for i=1:nf
        
        temp_prob0 = xxx(t).hist_fixed(i,temp_index_in);
        
        if temp_prob0 == 0
    
           % distribute a share of prob to non-zero cell, so the zero cell
           % has 1%
           temp_prob = xxx(t).hist_fixed(i,:);
           temp_nonzero = temp_prob ~=0;
           temp_prob(temp_nonzero) = temp_prob(temp_nonzero) - 1/sum(temp_nonzero);
           temp_prob(temp_index_in) = 1;
           
           % update
           xxx(t).hist_fixed_nozero(i,:) = temp_prob;
        end
    end
end

tab_prob = [];
temp_xlab = {};

for t = 1:1:size(xxx,1)-1
    temp_index_in = (xxx(t).histx_fixed(:,1) < xxx(t).actual) & (xxx(t).histx_fixed(:,2) >= xxx(t).actual);
    temp_prob = xxx(t).hist_fixed_nozero(:,temp_index_in);
    tab_prob = [tab_prob; temp_prob'];
    temp_xlab = [temp_xlab; xxx(t).sdate];
end

%%
fig = figure(123);
setmyfig(fig, [2,2,10,4]);

nf = size(tab_prob,2);
plot(sum(tab_prob==0,2), '*-', 'linewidth', 2)
hold on
plot(-1:86,nf*ones(88,1), '--', 'linewidth', 3);
hold off

temp_xtic = (1:16:84)';
set(gca, 'Xtick', temp_xtic, 'XtickLabel', temp_xlab(temp_xtic));
xlim([-1, 86]);
ylim([-1,nf+5]);
set(gca, 'fontsize', 15, 'linewidth', 2);
grid on
ylabel('# of forecasters');

%% Save updated dataset
dname = ['ecbspf_', varname, '_bp_nozero'];
eval([dname, '=xxx;']);

cd(datapath)
save(['ecbspf_', varname, '_bp_nozero.mat'], dname);
cd(workpath)

close all