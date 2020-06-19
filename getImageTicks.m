function [tickLocs, tickLabels] = getImageTicks(n,image_width)


tickLocs = [1 n/2 n];                   % locs of the tick labels (in pixels)
tickLabels = (0:image_width/2:image_width)*1e3;     % labels of ticks
tickLabels = cellstr(num2str(tickLabels(:)));