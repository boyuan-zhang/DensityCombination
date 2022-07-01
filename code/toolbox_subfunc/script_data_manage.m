% script used to construct the data
% among others, it requires

% ==========================================
%     % location of inflation data
%     ind_head = indx_r(infl_loc(:,1)) + 1;
%     ind_0    = indx_r(infl_loc(:,1)) + 2;
%     ind_1    = indx_r(cinfl_loc(:,1)) -2;
% ==========================================

% Acutal script
ind_id    = 2; %column location of id
ind_point = 3; %column location of point forecast
temp      = indx_c(cellfun(@isstr, c(ind_head,:)));
ind_hist  = temp(4:end); %columan location of histgram

temp_data = c(ind_0:ind_1,:); %data portion
x = c(ind_0:ind_1,1);
xstr = make_num2str(x);
xindx = (1:1:length(x))';

% extracting target for
%"rolling-based" forecasts (1-year-ahead, 2-year-ahead)
%"calende-year-based" forecasts (1-year, 2-year, 5-year)
temp_target   = unique(xstr);
temp_target_l = cellfun(@length, temp_target);
temp_target_r = temp_target(temp_target_l~=4,1); %rolling
temp_target_c = temp_target(temp_target_l==4,1); %calendar

% temp_indx1 = ~cellfun(@isempty, strfind(xstr, temp_target_r(1))); %1-year-ahead, rolling
% temp_indx2 = ~cellfun(@isempty, strfind(xstr, temp_target_r(2))); %2-year-ahead, rolling
% temp_indx3 = ~cellfun(@isempty, strfind(xstr, temp_target_c(1))); %1-year-ahead, calendar
% temp_indx4 = ~cellfun(@isempty, strfind(xstr, temp_target_c(2))); %2-year-ahead, calendar
% temp_indx5 = ~cellfun(@isempty, strfind(xstr, temp_target_c(end))); %5-year-ahead, calendar

temp_indx1 = strcmp(xstr, temp_target_r(1)); %1-year-ahead, rolling
temp_indx2 = strcmp(xstr, temp_target_r(2)); %2-year-ahead, rolling
temp_indx3 = strcmp(xstr, temp_target_c(1)); %1-year-ahead, calendar
temp_indx4 = strcmp(xstr, temp_target_c(2)); %2-year-ahead, calendar
temp_indx5 = strcmp(xstr, temp_target_c(end)); %5-year-ahead, calendar


