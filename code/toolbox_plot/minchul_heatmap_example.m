% heatmap example


%% Using external function
figure(12)
addpath('brewer');
m1 = flipud(cbrewer('seq', 'Blues', 4));
m2 = cbrewer('seq', 'Reds', 4);

pcolor(cvt_all')
colormap([m1; m2])
colorbar;
caxis([0,1]); %change color axis

%% Using my functions
figure(13)
mymap = [ ...
    rgb('FireBrick'); ...
    rgb('lightsalmon'); ...
    rgb('LightSkyBlue'); ...
    rgb('DarkBlue'); ...
    ];
pcolor(cvt_all')
colormap(mymap)
colorbar;
caxis([0,1]) ; %change color axis