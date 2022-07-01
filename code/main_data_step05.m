% Data description 
% inflation and interest rate
% date: 2020/8/9
% minchul shin

% desc_v00
% figures for slides 
% scores are higher the better

clc; clear all; close all;
workpath = pwd;
datapath = '../data/';

addpath(datapath);
% savepath = [workpath, filesep, 'fig_paper_data'];

addpath(genpath('toolbox_fcst'));
addpath(genpath('toolbox_plot'))

%% Time series -- Nominal interest rate
% FRED Graph Observations
% Federal Reserve Economic Data
% Link: https://fred.stlouisfed.org
% Help: https://fred.stlouisfed.org/help-faq
% Economic Research Division
% Federal Reserve Bank of St. Louis
% 
% EUR12MD156N
% 12-Month London Interbank Offered Rate (LIBOR), based on Euro, Percent, Daily, Not Seasonally Adjusted

% There are occasional zeros, we take out them when we take the average


% interest rate dates
[~, ni_date, ~] = xlsread('EUR12MD156N.xls', 'A12:A5694');

% interest rate series
[num, txt, raw] = xlsread('EUR12MD156N.xls', 'B12:B5694');

% survey dates etc 
% survey date, sent out date, deadline date
[num2,txt2,raw2] = xlsread('SPF_rounds_dates.xlsx', 'A5:C91');

% ----
% calculate average interest rate 
ns = size(txt2,1);
cell_ni_data = cell(ns, 4); %survey date, average interest rate
mat_ni_data = zeros(ns,1); %average interest rate only

% loop over survey dates
temp_xlab = {};
for t = 1:ns
    temp_s0 = datenum(txt2{t, 2}); % sent out date
    temp_d0 = datenum(txt2{t, 3}); % deadline date
    temp_indx = (temp_s0 <= datenum(ni_date)) & (temp_d0 >= datenum(ni_date));
    
    % plug-in info
    cell_ni_data(t,1:3) = txt2(t,:);
    cell_ni_data{t,4}   = mean(num(temp_indx),'omitnan');
    mat_ni_data(t,1)    = mean(num(temp_indx),'omitnan');
    
    % for graph
    temp_xlab = [temp_xlab; txt2(t,1)];
end


%% Figure - Nominal interest rate 

fig = figure(2);
setmyfig(fig, [2,2,10,4]);

plot(zeros(size(mat_ni_data,1),1), '--', 'linewidth', 3, 'color', rgb('coral'));
hold on
plot(mat_ni_data, '*-', 'linewidth', 2, 'color', rgb('darkblue'));
hold off

temp_xtic = (1:16:84)';
set(gca, 'Xtick', temp_xtic, 'XtickLabel', temp_xlab(temp_xtic));
xlim([-1, ns+2]);
ylim([-1,6]);
set(gca, 'fontsize', 15, 'linewidth', 2);
grid on
ylabel('Nominal Interest Rate (%)');
xlabel('Survey date');

%% Save data file

cd(datapath)
save('data_nominal_rate.mat', 'mat_ni_data', 'cell_ni_data');
cd(workpath)

close all