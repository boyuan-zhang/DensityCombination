function tab_prob = ecb_histx(xname,xmode)
% dictionary for the intervals of histogram
% ecb prob forecasts
% last updated: 11/18/2017
% This function returns left and right cut-off of the xname

if nargin == 1
    xmode = 100;
end

% [name, left cut-off, right cut-off]
% HICP
infl_prob = {'TN2_0'     , -inf, -2.0; ...
'TN1_0'     , -inf, -1.0 ; ...
'FN2_0TN1_6', -2.0, -1.6 ; ...
'FN1_5TN1_1', -1.5, -1.1 ; ...
'FN1_0TN0_6', -1.0, -0.6 ; ...
'FN0_5TN0_1', -0.5, -0.1 ; ...
'T0_0'      , -inf , 0.0 ; ...
'F0_0T0_4'  ,   0.0, 0.4 ; ...
'F0_5T0_9'  ,   0.5, 0.9 ; ...
'F1_0T1_4'  ,   1.0, 1.4 ; ...
'F1_5T1_9'  ,   1.5, 1.9 ; ...
'F2_0T2_4'  ,   2.0, 2.4 ; ...
'F2_5T2_9'  ,   2.5, 2.9 ; ...
'F3_0T3_4'  ,   3.0, 3.4 ; ...
'F3_5T3_9'  ,   3.5, 3.9 ; ...
'F3_5'      ,   3.5, inf ; ...
'F4_0'      ,   4.0, inf ; ...
};

% RGDP
rgdp_prob = {'TN1_0', -inf, -1.0       ; ...
'TN6_0', -inf, -6.0       ; ...
'FN6_0TN5_6', -6.0, -5.6  ; ...
'FN5_5TN5_1', -5.5, -5.1  ; ...
'FN5_0TN4_6', -5.0, -4.6  ; ...
'FN4_5TN4_1', -4.5, -4.1  ; ...
'FN4_0TN3_6', -4.0, -3.6  ; ...
'FN3_5TN3_1', -3.5, -3.1  ; ...
'FN3_0TN2_6', -3.0, -2.6  ; ...
'FN2_5TN2_1', -2.5, -2.1  ; ...
'FN2_0TN1_6', -2.0, -1.6  ; ...
'FN1_5TN1_1', -1.5, -1.1  ; ...
'FN1_0TN0_6', -1.0, -0.6  ; ...
'FN0_5TN0_1', -0.5, -0.1  ; ...
'T0_0', -inf, 0.0         ; ...
'F0_0T0_4', 0.0, 0.4      ; ...
'F0_5T0_9', 0.5, 0.9      ; ...
'F1_0T1_4', 1.0, 1.4      ; ...
'F1_5T1_9', 1.5, 1.9      ; ...
'F2_0T2_4', 2.0, 2.4      ; ...
'F2_5T2_9', 2.5, 2.9      ; ...
'F3_0T3_4', 3.0, 3.4      ; ...
'F3_5T3_9', 3.5, 3.9      ; ...
'F4_0T4_4', 4.0, 4.4      ; ...
'F4_5T4_9', 4.5, 4.9      ; ...
'F4_0', 4.0, inf          ; ...
'F5_0', 5.0, inf          ; ...
};

if xmode == 1 %infl
    tab = infl_prob;
elseif xmode == 2 %rgdp
    tab = rgdp_prob;
else
    tab = [infl_prob; rgdp_prob];
end

tab_prob = [];
nx = length(xname);
for i=1:1:nx
%     indx = ~cellfun(@isempty, strfind(tab(:,1), xname{i}));
    indx = strcmp(tab(:,1), xname{i});
    temp = tab(indx,:);
    tab_prob = [tab_prob; [temp{2}, temp{3}]];
end

% there is a jump between two intervals ..., we fill this gap
for i=2:1:nx
    tab_prob(i-1,2) = tab_prob(i,1);
end


% UNEMP
% T5_5, -inf, 5.5
% T6_5, -inf, 6.5
% T7_0, -inf, 7.0
% T7_5, -inf, 7.5
% T9_0, -inf, 9.0
% F5_5T5_9, 5.5, 5.9
% F6_0T6_4, 6.0, 6.4
% F6_5T6_9, 6.5, 6.9 
% F7_0T7_4, 7.0, 7.4
% F7_5T7_9, 7.5, 7.9
% F8_0T8_4, 8.0, 8.4 
% F8_5T8_9, 8.5, 8.9 
% F9_0T9_4, 9.0, 9.4
% F9_5T9_9, 9.5, 9.9 
% F10_0T10_4, 10.0, 10.4
% F10_5T10_9, 10.5, 10.9 
% F11_0T11_4, 11.0, 11.4 
% F11_5T11_9, 11.5, 11.9 
% F12_0T12_4, 12.0, 12.4
% F12_5T12_9, 12.5, 12.9
% F13_0T13_4, 13.0, 13.4 
% F13_5T13_9, 13.5, 13.9 
% F14_0T14_4, 14.0, 14.4 
% F14_5T14_9, 14.5, 14.9
% F11_0, 11.0, inf
% F11_5, 11.5, inf
% F12_0, 12.0, inf
% F15_0, 15.0, inf
