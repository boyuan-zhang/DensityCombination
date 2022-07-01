% data management for ECP-SPF
% date: 2017-11-18 (revised: 2020-7-26)
% minchul shin
clc; clear all; close all;

workpath = pwd;
datapath = '../data/';

addpath(datapath);
addpath([datapath 'SPF_individual_forecasts']);
addpath(genpath('toolbox_subfunc'));

%% Preliminary -
% filename = survey date (in year-quarter format)
% start: 1999Q1.csv
% end  : 2017Q4.csv

surveydate = {};
filename = {};
counter = 1;
for yind = 1999:1:2019
    for qind = 1:1:4
        temp_name = [num2str(yind), 'Q', num2str(qind), '.csv'];
        
        filename = [filename; temp_name];
        surveydate = [surveydate; temp_name(1:6)];
        counter = counter + 1;
    end
end
nsurvey = size(filename,1); %number of total survey

%% Load
ecbspf_infl_1y = [];
ecbspf_infl_2y = [];
ecbspf_infl_1yc = [];
ecbspf_infl_2yc = [];
ecbspf_infl_5yc = [];

ecbspf_rgdp_1y = [];
ecbspf_rgdp_2y = [];
ecbspf_rgdp_1yc = [];
ecbspf_rgdp_2yc = [];
ecbspf_rgdp_5yc = [];

ecbspf_unemp_1y = [];
ecbspf_unemp_2y = [];

for sind = 1:1:nsurvey
    % sind = 70;
    temp_name = filename{sind,1};
    
    [a,b,c] = xlsread(filename{sind});
%     t = readtable(filename{sind});

    % location of inflation, gdp, unemployment
    [nr, nc] = size(c);
    indx_r = (1:1:nr)';
    indx_c = (1:1:nc)';
    
    infl_loc  = ~cellfun(@isempty, strfind(b, 'HICP'));
    cinfl_loc = ~cellfun(@isempty, strfind(b, 'CORE'));
    gdp_loc   = ~cellfun(@isempty, strfind(b, 'GROWTH EXPECTATIONS;'));
    unemp_loc = ~cellfun(@isempty, strfind(b, 'EXPECTED UNEMPLOYMENT RATE;'));
    
    % --------------------
    % inflation
    
    % location of inflation data
    ind_head = indx_r(infl_loc(:,1)) + 1;
    ind_0    = indx_r(infl_loc(:,1)) + 2;
    ind_1    = indx_r(cinfl_loc(:,1)) -2;
    script_data_manage;

    % ---
    % 1-year-ahead, rolling
    temp_data00 = temp_data(temp_indx1,:);
    temp_target00 = temp_target_r{1};
    ecbspf_infl_1y = fill_data_survey(ecbspf_infl_1y, sind, temp_data00, temp_target00, surveydate, c, ind_head, ind_hist, 1);
    
    % ---
    % 2-year-ahead, rolling
    temp_data00 = temp_data(temp_indx2,:);
    temp_target00 = temp_target_r{2};
    ecbspf_infl_2y = fill_data_survey(ecbspf_infl_2y, sind, temp_data00, temp_target00, surveydate, c, ind_head, ind_hist, 1);
    
    
    % ---
    % 1-year-ahead, calendar
    temp_data00 = temp_data(temp_indx3,:);
    temp_target00 = temp_target_c{1};
    ecbspf_infl_1yc = fill_data_survey(ecbspf_infl_1yc, sind, temp_data00, temp_target00, surveydate, c, ind_head, ind_hist, 1);
    
    % ---
    % 2-year-ahead, calendar
    temp_data00 = temp_data(temp_indx4,:);
    temp_target00 = temp_target_c{2};
    ecbspf_infl_2yc = fill_data_survey(ecbspf_infl_2yc, sind, temp_data00, temp_target00, surveydate, c, ind_head, ind_hist, 1);
    
    % ---
    % 5-year-ahead, calendar
    temp_data00 = temp_data(temp_indx5,:);
    temp_target00 = temp_target_c{end};
    ecbspf_infl_5yc = fill_data_survey(ecbspf_infl_5yc, sind, temp_data00, temp_target00, surveydate, c, ind_head, ind_hist, 1);
   
    
    % --------------------
    % RGDP
    
    % location of rgdp data
    ind_head = indx_r(gdp_loc(:,1)) + 1;
    ind_0    = indx_r(gdp_loc(:,1)) + 2;
    ind_1    = indx_r(unemp_loc(:,1)) -2;
    script_data_manage;
    
    % ---
    % 1-year-ahead, rolling
    temp_data00 = temp_data(temp_indx1,:);
    temp_target00 = temp_target_r{1};
    ecbspf_rgdp_1y = fill_data_survey(ecbspf_rgdp_1y, sind, temp_data00, temp_target00, surveydate, c, ind_head, ind_hist, 2);
    
    % ---
    % 2-year-ahead, rolling
    temp_data00 = temp_data(temp_indx2,:);
    temp_target00 = temp_target_r{2};
    ecbspf_rgdp_2y = fill_data_survey(ecbspf_rgdp_2y, sind, temp_data00, temp_target00, surveydate, c, ind_head, ind_hist, 2);
    
    % ---
    % 1-year-ahead, calendar
    temp_data00 = temp_data(temp_indx3,:);
    temp_target00 = temp_target_c{1};
    ecbspf_rgdp_1yc = fill_data_survey(ecbspf_rgdp_1yc, sind, temp_data00, temp_target00, surveydate, c, ind_head, ind_hist, 2);
    
    % ---
    % 2-year-ahead, calendar
    temp_data00 = temp_data(temp_indx4,:);
    temp_target00 = temp_target_c{2};
    ecbspf_rgdp_2yc = fill_data_survey(ecbspf_rgdp_2yc, sind, temp_data00, temp_target00, surveydate, c, ind_head, ind_hist, 2);
    
    % ---
    % 5-year-ahead, calendar
    temp_data00 = temp_data(temp_indx5,:);
    temp_target00 = temp_target_c{end};
    ecbspf_rgdp_5yc = fill_data_survey(ecbspf_rgdp_5yc, sind, temp_data00, temp_target00, surveydate, c, ind_head, ind_hist, 2);
    

end

%% Save
cd(datapath)
save('data_ecb_spf_2019Q4.mat', ...
    'ecbspf_infl_1y', 'ecbspf_infl_2y', 'ecbspf_infl_1yc', 'ecbspf_infl_2yc', 'ecbspf_infl_5yc', ...
    'ecbspf_rgdp_1y', 'ecbspf_rgdp_2y', 'ecbspf_rgdp_1yc', 'ecbspf_rgdp_2yc', 'ecbspf_rgdp_5yc');
cd(workpath)




