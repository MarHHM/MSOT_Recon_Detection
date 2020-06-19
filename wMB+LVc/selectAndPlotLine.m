function linePosOut = selectAndPlotLine(linePosIn)
% - opens a saved .fig file & selects a line on it to draw its profile
%
% - if no line coords are input, you have to draw the line
% - only vertical & horizontal lines are allowed, inclined lines are approximated according to thier slope

% % choose & plot the image on which you will select the linePlot (e.g. the ROI from the prev section)
% imHandel = figure, imagesc(roi.Im), colormap(gray), axis image off;
[fName, pName, ~] = uigetfile({'*.fig'}, 'please select a .fig file', 'C:\Users\marwan.muhammad\Desktop');
imHandel = openfig([pName fName]);

if ~exist('linePosIn','var')
    h = imline;
    linePos = wait(h);
    delete(h);
else
    linePos = linePosIn;
end
linePlot.x0 = floor(linePos(1,1));
linePlot.y0 = floor(linePos(1,2));
linePlot.xf = floor(linePos(2,1));
linePlot.yf = floor(linePos(2,2));
    
line([linePos(1,1) linePos(2,1)], [linePos(1,2) linePos(2,2)], ...
     'Color', 'y',...
     'LineStyle', '--');
temp = getimage(imHandel);
figure;
if atan(abs(linePlot.yf-linePlot.y0)/abs(linePlot.xf-linePlot.x0)) < pi/4           % detect if line is vertical or horizontal
    plot([linePlot.x0:linePlot.xf], temp(linePlot.y0,linePlot.x0:linePlot.xf)), xlabel('pixel number'), ylabel('OA intensity (au)');
else
    plot([linePlot.y0:linePlot.yf], temp(linePlot.y0:linePlot.yf,linePlot.x0)), xlabel('pixel number'), ylabel('OA intensity (au)');
end

linePosOut = linePos;