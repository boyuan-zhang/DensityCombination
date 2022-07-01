function xxx = fill_data_actual(xxx,yyy,freq)
% Fill actual into the survey data
% freq = 12 (monthly)
% freq = 4 (quarterly)

% xxx = ecbspf_infl_1y;
% yyy = actual_infl;
zzz = [];
zzz2 = [];
nobs = size(xxx,1);
for sind = 1:1:nobs;
    
    if (length(xxx(sind).tdate) == 4)
        if freq == 4
            temp_tdate = [xxx(sind).tdate, 'Q4'];
        elseif freq == 12
            temp_tdate = [xxx(sind).tdate, 'Dec'];
        end
    else
        temp_tdate = xxx(sind).tdate;
    end
    temp_ind = strcmp(yyy(:,1), temp_tdate);
%     temp_ind = strcmp(yyy(:,1), xxx(sind).tdate);
    
    if sum(temp_ind) == 0
        xxx(sind).actual = nan;
        zzz = [zzz; nan];
        zzz2 = [zzz2; nan];
    else
        xxx(sind).actual = yyy{temp_ind,2};
        zzz = [zzz; yyy{temp_ind,2}];
        zzz2 = [zzz2; mean(xxx(sind).point)];
    end
end
% ecbspf_infl_1y = xxx;