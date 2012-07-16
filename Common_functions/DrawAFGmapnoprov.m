% Draw the AFG map with no provinces

ax = worldmap([29,39],[60,75]);
setm(ax,'Grid','off','Frame','off','meridianlabel','off','parallellabel','off')
set(gcf,'color',[1,1,1])
geoshow(ax, Countrybounds, 'FaceColor', [0.93 0.87 0.51])
