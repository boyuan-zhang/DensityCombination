function []=shadenber(blah)
% shadenber.m
%  Routine to shade the nberdates in a figure
% Taken from: http://www.stanford.edu/~chadj/ChadMatlabFunctions.html
% Modified by Minchul Shin (make it self contained function)
% Last updates: 7/18/2013
% NOTE
% - x-axis should follow the date convention like:
%   date convention: year + (month-1)/12
% EXAMPLE
% xcal = ((1959+(3-1)/12):1/12:(2013+(7-1)/12))';
% x    = randn(length(xcal,1));
% plot(xcal,x);
% shadenber;

%% Define NBER dates
%  Official NBER recession dates, from
%   http://www.nber.org/cycles/cyclesmain.html  (see NBERcycles.html)
dates=[
% Peak Year,Month   Trough Year, Month
1857 6      1858 12;
1860 10     1861 6 ;
1865 4      1867 12;
1869 6      1870 12;
1873 10     1879 3 ;
1882 3      1885 5 ;
1887 3      1888 4 ;
1890 7      1891 5 ;
1893 1      1894 6 ;
1895 12     1897 6 ;
1899 6      1900 12;
1902 9      1904 8 ;
1907 5      1908 6 ;
1910 1      1912 1 ;
1913 1      1914 12;
1918 8      1919 3 ;
1920 1      1921 7 ;
1923 5      1924 7 ;
1926 10     1927 11;
1929 8      1933 3 ;
1937 5      1938 6 ;
1945 2      1945 10;
1948 11     1949 10;
1953 7      1954 5 ;
1957 8      1958 4 ;
1960 4      1961 2 ;
1969 12     1970 11;
1973 11     1975 3 ;
1980 1      1980 7 ;
1981 7      1982 11;
1990 7      1991 3 ;
2001 3      2001 11;
2007 12     2009 06;  % Starts Dec 2007, End=June 2009
      ];

start  = dates(:,1)+(dates(:,2)-1)/12;
finish = dates(:,3)+(dates(:,4)-1)/12;


%% Plot NBER recession shade
curax=axis;
indx1=find(finish>curax(1));  % First recession to include;
indx2=find(start<curax(2));  % Last recession to include;
indx1=indx1(1);
indx2=indx2(length(indx2));
if start(indx1)<curax(1);
  start(indx1)=curax(1);
end;
if finish(indx2)>curax(2);
  finish(indx2)=curax(2);
end;

% colorstr=[159 182 205]/256;
colorstr = rgb('silver');
% colorstr = [0.9,0.7,0.7];

shade(start(indx1:indx2),finish(indx1:indx2),colorstr);
end


% function used here in this file
function []=shade(start,finish,colorstr)

% function []=shade(start,finish,colorstr);
%
%  start and finish are Nx1 vectors of starting and ending years.
%  The function shades between the start and finish pairs using colorstr

if ~exist('colorstr'); colorstr='y'; end;  % default is yellow
curax=axis;
y=[curax(3) curax(4) curax(4) curax(3)];
hold on;
for i=1:length(start);
  x=[start(i) start(i) finish(i) finish(i)];
  fill(x,y,colorstr);
end;
  
% Now, prevent the shading from covering up the lines in the plot.  
h = findobj(gca,'Type','line');
% set(h,'EraseMode','xor');
set(h,'EraseMode','background');
 
h = findobj(gca,'Type','patch');
set(h,'EdgeColor','none');
alpha(h, 0.7);


% This last one makes the tick marks visible
set(gca, 'Layer', 'top')
end

