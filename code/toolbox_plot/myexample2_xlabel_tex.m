% writing xlabel, ylabel with latex interpreter


xxx = [1,2,3];

plot(xxx)
% xTicks = set(gca, 'xtick', [1,2,3], 'xticklabel', {'$\sigma$', 'b', 'c'}, 'interpreter', 'latex');
set(gca, 'xticklabel', [], 'yticklabel', []);
set(gca, 'xtick', [1,2.5,3]);
yTicks = get(gca,'ytick');
xTicks = get(gca, 'xtick');

ax = axis; %Get left most x-position
HorizontalOffset = 0.1;
% Reset the xtick labels in desired font 
minY = min(yTicks);
verticalOffset = 0.1;

myticklabel = {'$\sigma$', 'b', 'c'};
for xx=1:1:3
text(xTicks(xx), minY - verticalOffset, myticklabel{xx},...
'HorizontalAlignment','Right','interpreter', 'latex'); 
end