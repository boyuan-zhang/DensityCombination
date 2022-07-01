function sss = xcaltrans(xxcal)
% transform xcal (number) to string 
% only supports monthly

% example: 1999M3
% xxcal = 1999 + (3-1)/12;

% year part
yyy = floor(xxcal);

% month part
mmm = (xxcal - yyy)*12 + 1;

% make string
sss = cell(length(xxcal),1);
for i=1:1:length(xxcal)
    sss{i,1} = [num2str(yyy(i,1)), 'M', num2str(mmm(i,1))];
end