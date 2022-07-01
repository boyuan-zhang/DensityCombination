% data management for ECP-SPF
% date: 2017-11-18 (revised: 2020-07-27)
% minchul shin

% preliminary analysis
clc; clear all; close all;
workpath = pwd;
datapath = '../data/';

addpath(datapath);
addpath(genpath('toolbox_plot'));

%% load data
% forecast
load('data_ecb_spf_2019Q4_v02.mat');

% pick a dataset
dname = 'ecbspf_infl_1y';
eval(['xxx = ', dname, ';']);

%% Checking bins
% how many realizations are outside of interval bins ...
ns0 = size(xxx,1); %

% computing actual number of surveys (surveys when actuals are available)
ns = 0;
for sind = 1:1:ns0
    if ~isnan(xxx(sind).actual)
        ns = ns+1;
    end
end

% reterive information
y0 = zeros(ns,1); %lowest bin
y1 = zeros(ns,1); %highest bin
ya = zeros(ns,1); %actual
nbins   = zeros(ns,1);
locbins = zeros(ns,1);
for sind = 1:1:ns;
    p  = xxx(sind).hist;
    px = xxx(sind).histx;
    y  = xxx(sind).actual;
    temp_ind = find((px(:,1)<y)&(y<=px(:,2)));
    
    nbins(sind,1) = size(px,1);
    locbins(sind,1) = temp_ind;
    
    
    y0(sind,1) = px(1,2);
    y1(sind,1) = px(end,1);
    ya(sind,1) = y;
    
end

figure
plot([ones(ns,1),nbins], 'r','linewidth', 2)
hold on
plot(locbins, 'b', 'linewidth',2)
hold off
title('red: bounded intervals, blue: location of realization');


fig = figure(2);
setmyfig(fig, [1.7, 1.2, 8, 4]);
plot([y0,y1], 'r', 'linewidth',2);
hold on
plot(ya, 'b', 'linewidth',2);
hold off
set(gca, 'linewidth',2, 'fontsize', 20);
xlabel('Time');
title('Range of intervals');


fig = figure(3);
setmyfig(fig, [1.7, 1.2, 8, 4]);
plot(nbins, 'linewidth',3)
ylim([7,23]);
set(gca, 'linewidth',2, 'fontsize', 20);
xlabel('Time');
title('number of bins');


%% Unique number of ids and their appearances
% how many realizations are outside of interval bins ...
ns0 = size(xxx,1); %

% computing actual number of surveys (surveys when actuals are available)
ns = 0;
for sind = 1:1:ns0
    if ~isnan(xxx(sind).actual)
        ns = ns+1;
    end
end

% reterive information
temp_id = [];
for sind = 1:1:ns
    temp_id  = [temp_id; xxx(sind).id_hist];
end

% total responses
disp(size(temp_id));

% total unique ids
disp(size(unique(temp_id)));

% unique id
id_unq = unique(temp_id);
id_n = length(id_unq);

id_res = zeros(ns,id_n);
for sind = 1:1:ns
    temp_id = xxx(sind).id_hist;
    for i=1:1:id_n
        if ismember(id_unq(i), temp_id)
            id_res(sind,i) = 1;
        end
    end
end

% -----------------------
% forecasters selection
% The filter is such that forecasters with more
% than four consecutive missing observations are excluded
% from the panel. (Genre, Kenny, Meyler, Timmermann)
KK  = 5;
KK0 = 0;
iniK = 4; %don't check this condition for the iniK quarters
id_sel = zeros(id_n,1);

for iind = 1:1:id_n
    temp_x = id_res(:,iind);
    
    % crit 1 = first KK0 exists
    valid1 = sum(temp_x(1:KK0,1)) == KK0;
    
    % crit 2 = keep only those who are not missing KK consecutive surveys
    valid2 = 1;
    temp_y = temp_x == 0;
    for i=(KK+iniK):1:size(temp_y)
        if sum(temp_y(i-KK+1:i,1)) == KK
            valid2 = 0;
            break;
        end
    end
    
    if valid1&&valid2
        id_sel(iind,1) = 1;
    end
end
disp('number of forecasters after filtering');
sum(id_sel)
id_in = id_unq(logical(id_sel));


% -------------------------------------
% construct a new panel
id_n = numel(id_in);
zzz = xxx; %copy
for sind = 1:1:ns
    temp_xxx = xxx(sind);
    temp_zzz = zzz(sind);
    
    % temp_xxx.hist
    % temp_xxx.id_hist
    
    temp_id = temp_xxx.id_hist; %id in the original data
    nbins = size(temp_xxx.hist,2);
    
    temp_zzz.hist = nan*ones(id_n,nbins);
    temp_zzz.id_hist = id_in;
    for i=1:1:id_n
        if ismember(id_in(i), temp_id)
            temp_zzz.hist(i,:) = temp_xxx.hist(temp_id == id_in(i),:);
        end
    end
    
    % add average hist (all forecasters)
    temp_zzz.hist_avg_a = mean(temp_xxx.hist,1,'omitnan');
    
    % add average hist (filtered forecasters)
    temp_zzz.hist_avg_s = mean(temp_zzz.hist,1,'omitnan');
    
    % store
    zzz(sind).hist       = temp_zzz.hist;
    zzz(sind).id_hist    = temp_zzz.id_hist;
    zzz(sind).hist_avg_a = temp_zzz.hist_avg_a;
    zzz(sind).hist_avg_s = temp_zzz.hist_avg_s;
end

% example figure
figure
plot(temp_zzz.hist', 'k', 'linewidth',1)
hold on
plot(temp_zzz.hist_avg_s', 'r', 'linewidth',3)
plot(temp_zzz.hist_avg_a', 'b', 'linewidth',3)
hold off

% -------------------------------------
% extrapolation
% 1) t=1, we replace with average forecasts
% 2) t>1, we replace with the similar forecasters (ordered by ranked score, by five group, based on previous survey performance)
sind = 1;
temp_zzz = zzz(sind);
temp_nan = isnan(temp_zzz.hist(:,1));
temp_zzz.hist(temp_nan,:) = repmat(temp_zzz.hist_avg_s, sum(temp_nan),1);
zzz(sind).hist = temp_zzz.hist;

% rps
temp_n = numel(temp_nan);
zzz(sind).rps = zeros(temp_n,1);
for i = 1:1:temp_n
    temp_p = temp_zzz.hist(i,:)';
    temp_x = temp_zzz.histx(:,:);
    temp_y = temp_zzz.actual;
    temp_rps = -sum( (cumsum(temp_p)/100 - (temp_y<=temp_x(:,2))).^2 ); %rps from Czado, Genieting ...
    zzz(sind).rps(i,1) = temp_rps;
end

% rps ranking
[a,b] = sort(zzz(sind).rps, 'descend');
zzz(sind).rps_ranking = zeros(temp_n,1);
zzz(sind).rps_ranking(b,1) = (1:1:temp_n)';
% (Ngr) groups by rps
Ngr = 5;
temp_nn = round(temp_n/Ngr);
% zzz(sind).rps_group(b,1) = [ones(temp_nn,1); 2*ones(temp_nn,1); 3*ones(temp_n-2*temp_nn,1)];
temp_rps_group = ones(temp_nn,1);
for igr = 2:(Ngr-1)
    temp_rps_group = [temp_rps_group; igr*ones(temp_nn,1)];
end
temp_rps_group = [temp_rps_group; Ngr*ones(temp_n-(Ngr-1)*temp_nn,1)];
zzz(sind).rps_group(b,1) = temp_rps_group;

% t>1
for sind = 2:1:ns
    temp_zzz = zzz(sind);
    temp_nan = isnan(temp_zzz.hist(:,1));
    temp_idloc = find(temp_nan);
    for i=1:1:numel(temp_idloc)
        %zzz(sind-1).rps(temp_idloc(i))
        temp_rg  = zzz(sind-1).rps_group;
        temp_rgi = temp_rg(temp_idloc(i));
        zzz(sind).hist(temp_idloc(i),:) = mean(zzz(sind).hist(temp_rg == temp_rgi,:),'omitnan');
    end
    
    % rps
    temp_n = numel(temp_nan);
    zzz(sind).rps = zeros(temp_n,1);
    for i = 1:1:temp_n
        temp_p = zzz(sind).hist(i,:)';
        temp_x = zzz(sind).histx(:,:);
        temp_y = zzz(sind).actual;
        temp_rps = -sum( (cumsum(temp_p)/100 - (temp_y<=temp_x(:,2))).^2 ); %rps from Czado, Genieting ...
        zzz(sind).rps(i,1) = temp_rps;
    end
    
    % rps ranking
    [a,b] = sort(zzz(sind).rps, 'descend');
    zzz(sind).rps_ranking = zeros(temp_n,1);
    zzz(sind).rps_ranking(b,1) = (1:1:temp_n)';
    
    % (Ngr) groups by rps
    Ngr = 3;
    temp_nn = round(temp_n/Ngr);
    % zzz(sind).rps_group(b,1) = [ones(temp_nn,1); 2*ones(temp_nn,1); 3*ones(temp_n-2*temp_nn,1)];
    temp_rps_group = ones(temp_nn,1);
    for igr = 2:(Ngr-1)
        temp_rps_group = [temp_rps_group; igr*ones(temp_nn,1)];
    end
    temp_rps_group = [temp_rps_group; Ngr*ones(temp_n-(Ngr-1)*temp_nn,1)];
    zzz(sind).rps_group(b,1) = temp_rps_group;
end

close all

%% save
eval([dname, '_bp = zzz;']);

cd(datapath)
save([dname,'_bp.mat'], [dname, '_bp']);
cd(workpath)
