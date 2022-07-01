% script_heatmap_plot_diff.m
% -> generate heatmap version of diff (with uniform forecasters)

temp_x = xxx(1).histx_fixed;
temp_x(1,1) = temp_x(1,2)-0.5;
temp_x(end,2) = temp_x(end,1)+0.5;
temp_mid = temp_x(:,1) + (temp_x(:,2)-temp_x(:,1))/2;

big_grid = temp_mid;
ngrid = numel(big_grid);
p_grid = zeros(ns,ngrid);
m_grid = zeros(ns,1);
s_grid = zeros(ns,1);

for sind = ns0:1:ns1
    
    % Simple average
    temp_x = xxx(sind).histx_fixed;
    temp_x(1,1) = temp_x(1,2)-0.5;
    temp_x(end,2) = temp_x(end,1)+0.5;
    temp_mid = temp_x(:,1) + (temp_x(:,2)-temp_x(:,1))/2;
    temp_hist0 = xxx(sind).hist_fixed_nozero;
    temp_hist = mean(temp_hist0,1);
    
    temp_hist1 = mat_b1(sind,:) * temp_hist0;
    temp_hist2 = mat_b2(sind,:) * temp_hist0(1:end-1,:);
    
    % probability
    p_grid(sind,:) = temp_hist1 - temp_hist2;
    
    % mean and sd
    m_grid(sind,1) = (temp_hist/100)*temp_mid;
    s_grid(sind,1) = sqrt( ((temp_hist/100)*temp_mid.^2)-m_grid(sind,1).^2);    
end

% actual figure
mymap1 = (cbrewer('seq', 'Blues', 100));
mymap2 = (cbrewer('seq', 'Reds', 100));
mymap = [flipud(mymap1); [1,1,1]; mymap2];

imagesc(ns0:ns1, big_grid, p_grid(ns0:ns1,:)')

colormap(mymap)
cb = colorbar('eastoutside');
caxis([-25,25]) ; %change color axis
% xlabel('Time');

cb.Label.String = '(%)';
cb.Ticks = (-25:5:25);

xticks_id = 1:8:ns;
xticks_date = {};

for i = 1:numel(xticks_id)
   xticks_date{i} = xxx(xticks_id(i)).sdate;
end

set(gca,'YDir','normal')

xticks(xticks_id)
xticklabels(xticks_date)
yticks(0:4)
% set(gca,'TickLength',[0.1, 0.01], 'fontsize', 15)
set(gca,'fontsize', 16)

% view(-70,60) % fix the angles of the view

% xlabel('Survey date');
% ylabel('Inflation rate (%)')
% ylabel('Change in prob. and inflation rate (%)')
ylabel('Forecast Difference and Realized Inflation (%)')
zlabel('Probability');