% data management for ECP-SPF
% date: 2017-11-18 (revised: 2020-07-27)
% minchul shin

% this is to fill-in "actual" value
% we treat actual value as "final" (most recent vintage, which is downloaded in 2017-11-18)

% ****************
% revised version: we treat actual value as "final" (most recent vintage,
% which is downloaded in 2020-07-27)
% ****************

% actual: https://sdw.ecb.europa.eu/home.do?chart=t1.3

clc; clear all; close all;
workpath = pwd;
datapath = '../data/';

addpath(datapath);

%% load data

% forecast
load('data_ecb_spf_2019Q4.mat');

% actual
[~,~,temp_rgdp] = xlsread('data_rgdp_20200727.csv');
[~,~,temp_infl] = xlsread('data_hicp_20200727.csv');

actual_rgdp = temp_rgdp(6:end,:);
actual_infl = temp_infl(6:end,:);

%% Fill actual into the survey data

% inflation
ecbspf_infl_1y = fill_data_actual(ecbspf_infl_1y,actual_infl,12);
ecbspf_infl_2y = fill_data_actual(ecbspf_infl_2y,actual_infl,12);
ecbspf_infl_1yc = fill_data_actual(ecbspf_infl_1yc,actual_infl,12);
ecbspf_infl_2yc = fill_data_actual(ecbspf_infl_2yc,actual_infl,12);
ecbspf_infl_5yc = fill_data_actual(ecbspf_infl_5yc,actual_infl,12);

% rgdp 
ecbspf_rgdp_1y = fill_data_actual(ecbspf_rgdp_1y,actual_rgdp,4);
ecbspf_rgdp_2y = fill_data_actual(ecbspf_rgdp_2y,actual_rgdp,4);
ecbspf_rgdp_1yc = fill_data_actual(ecbspf_rgdp_1yc,actual_rgdp,4);
ecbspf_rgdp_2yc = fill_data_actual(ecbspf_rgdp_2yc,actual_rgdp,4);
ecbspf_rgdp_5yc = fill_data_actual(ecbspf_rgdp_5yc,actual_rgdp,4);


%% Save
cd(datapath)
save('data_ecb_spf_2019Q4_v02.mat', ...
    'ecbspf_infl_1y', 'ecbspf_infl_2y', 'ecbspf_infl_1yc', 'ecbspf_infl_2yc', 'ecbspf_infl_5yc', ...
    'ecbspf_rgdp_1y', 'ecbspf_rgdp_2y', 'ecbspf_rgdp_1yc', 'ecbspf_rgdp_2yc', 'ecbspf_rgdp_5yc');
cd(workpath)
