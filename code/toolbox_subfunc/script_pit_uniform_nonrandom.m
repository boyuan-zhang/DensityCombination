% script_pit_uniform_nonrandom.m
% -> generate histogram for nonrandom version of PIT

% Figure 
% fig = figure(2001);
% setmyfig(fig, [2,2,6.5,5]);

mid_edges = (edges(1:end-1) + edges(2:end)) / 2;
bb = bar(mid_edges, avg_bar_nonrandom);
bb.EdgeColor = rgb('silver');
bb.FaceColor = rgb('silver');
hold on
p1 = plot([0,1], [(1/nb), (1/nb)], '-', 'linewidth', 1.5, 'color', rgb('black'));
p1.Color(4) = 0.5;
se_line = 2*sqrt(1/nb * (1-1/nb) / (ns1-ns0+1));
p1 = plot([0,1], (1/nb)+ [se_line, se_line], '--', 'linewidth', 3, 'color', rgb('red'));
p2 = plot([0,1], (1/nb)+ [-se_line, -se_line], '--', 'linewidth', 3, 'color', rgb('red'));
p1.Color(4) = 0.5;
p2.Color(4) = 0.5;
hold off

grid on
% ylim([0,0.4])
% ylim([0,0.45]) % only for the version for subperiods
ylim([0,0.5]) % only for the version for subperiods, second break point
xlim([0,1])

yticks(0:0.1:0.4)

set(gca, 'fontsize', 20, 'linewidth', 2);
