% script_pit_uniform.m
% -> generate histogram for PIT

nb = 5; %number of bins on [0,1]
edges = linspace(0,1,nb+1)';

mat_bar = zeros(nsim, nb, ns);
for sind = ns0:ns1
    temp_p = mat_pit_random(sind,:);
    for bind = 1:nb
        mat_bar(:,bind,sind) = edges(bind) < temp_p & edges(bind+1) >= temp_p;
    end
end
% hist(mat_pit_random(:))
avg_mat_bar = mean(mat_bar(:,:,ns0:ns1),3); %average over time
avg_bar = mean(avg_mat_bar,1); %average over simulation

mid_edges = (edges(1:end-1) + edges(2:end)) / 2;
bb = bar(mid_edges, avg_bar);
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
