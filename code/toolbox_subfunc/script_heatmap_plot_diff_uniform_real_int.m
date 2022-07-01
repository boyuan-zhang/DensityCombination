% script_heatmap_plot_diff_uniform_real_int.m
% -> generate heatmap version of diff (with uniform forecasters)
% -> for real interest rate.

temp_grid = [];

for sind=ns0:1:ns
    
    temp_x = xxx(sind).histx_fixed - mat_ni_data(sind);
    temp_x(1,1) = temp_x(1,2)-0.5;
    temp_x(end,2) = temp_x(end,1)+0.5;
    
    temp_x = - temp_x;
    
    temp_mid = temp_x(:,1) + (temp_x(:,2)-temp_x(:,1))/2;

    temp_grid = [temp_grid; temp_mid];
end

big_grid = unique(temp_grid);
ngrid = numel(big_grid);
p_grid = zeros(ns,ngrid);

for sind = ns0:1:ns
    
    % grid x
    temp_x = xxx(sind).histx_fixed - mat_ni_data(sind);
    temp_x(1,1) = temp_x(1,2)-0.5;
    temp_x(end,2) = temp_x(end,1)+0.5;
    temp_x = -temp_x;
    temp_mid = temp_x(:,1) + (temp_x(:,2)-temp_x(:,1))/2;
    temp_hist0 = xxx(sind).hist_fixed_nozero;
    
    temp_hist1 = mat_b1(sind,:) * temp_hist0;
    temp_hist2 = mat_b2(sind,:) * temp_hist0(1:end-1,:);
    
    temp_hist = temp_hist1 - temp_hist2;
    
    % avereage probability
    for j=1:1:numel(temp_mid)
        loc_j = find(big_grid == temp_mid(j));
        p_grid(sind,loc_j) = temp_hist(j);
    end
 
end


% Adjust big_grid and p_grid
grid_max = round(max(big_grid)*4)/4;
grid_min = round(min(big_grid)*4)/4;

if max(big_grid) > grid_max || min(big_grid) < grid_min
   grid_max = grid_max + 0.25;
   grid_min = grid_min + 0.25;
end

% big_grid_new = grid_min:0.5:grid_max;
big_grid_new = grid_min:0.5:grid_max;
ngrid_new = numel(big_grid_new);
p_grid_new = zeros(ns, ngrid_new);

% assign each midpoint in big_grid to the nearest point in the new and small grid 
grid_convert_id = zeros(ngrid,1);
for i = 1:1:ngrid
   temp_dist = abs( big_grid(i) - big_grid_new );
   [a, b] = min(temp_dist);
   grid_convert_id(i) = b;
end

% aggregate prob
for sind = ns0:1:ns
   for i = 1:1:ngrid_new
      temp_hist = p_grid(sind, grid_convert_id == i);
      p_grid_new(sind, i) = sum(temp_hist);
   end
end



% actual figure
mymap1 = (cbrewer('seq', 'Blues', 100));
mymap2 = (cbrewer('seq', 'Reds', 100));
mymap = [flipud(mymap1); [1,1,1]; mymap2];

imagesc(ns0:ns1, big_grid_new, p_grid_new(ns0:ns1,:)')

colormap(mymap)
cb = colorbar('eastoutside');
caxis([-30,30]) ; %change color axis
% xlabel('Time');

cb.Label.String = '(%)';
cb.Ticks = (-30:5:30);

xticks_id = 1:8:ns;
xticks_date = {};

for i = 1:numel(xticks_id)
   xticks_date{i} = xxx(xticks_id(i)).sdate;
end

set(gca,'YDir','normal')

xticks(xticks_id)
xticklabels(xticks_date)
yticks(-4:6)
% set(gca,'TickLength',[0.1, 0.01], 'fontsize', 15)
set(gca,'fontsize', 16)

% view(-70,60) % fix the angles of the view

% xlabel('Survey date');
% ylabel('Inflation rate (%)')
% ylabel('Change in prob. and real interest rate (%)')
ylabel('Forecast Difference and Realized Real Rate (%)')
zlabel('Probability');