function ecbspf00 = fill_data_survey(ecbspf00, sind, ...
    temp_data00, temp_target00, surveydate, c, ind_head, ind_hist, hist_xmode)
% fill in survey data
% 2020-7-27

% ---
% input
% ecbspf00
% sind
% temp_data00 = temp_data(temp_indx1,:);
% temp_target00 = temp_target_r{1}
% surveydate
% c
% ind_head
% ind_hist
% hist_xmode = 1 if infl, 2 if rgdp

% ---
% actual computation
ecbspf00(sind,1).sdate = surveydate{sind}; %survey date
ecbspf00(sind,1).tdate = temp_target00; %target date
ecbspf00(sind,1).id = cell2mat(temp_data00(:,2)); %forecaster id
ecbspf00(sind,1).point = cell2mat(temp_data00(:,3)); %point forecast

% prob fcst
temp_prob = cell2mat(temp_data00(:,ind_hist));
temp_prob(isnan(temp_prob)) = 0;

ecbspf00(sind,1).hist = temp_prob; %histogram for prob
ecbspf00(sind,1).histx_name = c(ind_head, ind_hist); %x-axis of the histogram
ecbspf00(sind,1).histx = ecb_histx(ecbspf00(sind,1).histx_name,hist_xmode); %x-axis of the histogram

% keep only who reported

% .. of point
temp_point = ~isnan(ecbspf00(sind,1).point);
ecbspf00(sind,1).id_point = ecbspf00(sind,1).id(temp_point,1);
ecbspf00(sind,1).point    = ecbspf00(sind,1).point(temp_point,1);

% .. of probability
temp_hist0 = (ecbspf00(sind,1).hist);
temp_hist  = sum(temp_hist0,2)>10; %keep ids only sum of prob >10, this is becuase we set NaN to 0, there could be spurious zeros
ecbspf00(sind,1).id_hist  = ecbspf00(sind,1).id(temp_hist,1);
ecbspf00(sind,1).hist  = ecbspf00(sind,1).hist(temp_hist,:);

















