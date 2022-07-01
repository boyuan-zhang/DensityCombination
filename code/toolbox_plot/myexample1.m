% example of figure
for iter_N = 1:1:3
    % 95% confidence band
    se = zeros(T,2);
    se(:,1) = pf_s_up(:, iter_N) - 1.96*sqrt(pf_P_up(:,iter_N));
    se(:,2) = pf_s_up(:, iter_N) + 1.96*sqrt(pf_P_up(:,iter_N));
    
    % plotting     
    f1 = figure(iter_N);
    set(f1, 'color', 'w');
    set(f1, 'units', 'inches');
    set(f1, 'outerposition', [2.5, 1.1, 9, 6]);
    set(f1, 'paperpositionmode', 'auto');
    ax = axes('FontSize', 15);
    
    plot(s, 'linewidth', 2, 'linestyle', '--', 'color', 'r')
    hold on
    plot(pf_s_up(:,iter_N), 'linewidth', 1.5, 'color', 'b');
    legend1 = legend('True', 'Particle filter');
    set(legend1, 'location', 'southeast', 'fontsize', 15, 'fontname', 'Tahoma')
    hold on
    [ph,msg]=jbfill((1:1:T),se(:,1)',se(:,2)',rgb('SlateBlue'),rgb('Navy'),0,0.2);
    hold off
    axis([1, T, -4, 4])
    title(['s_{t} from Particle Filter with ', num2str(set_N(iter_N)), ...
        ' particles'], 'FontSize', 15, 'FontName', 'Tahoma')
    set(ax, 'LineWidth', 1.5);
end