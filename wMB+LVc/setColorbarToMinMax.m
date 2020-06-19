function [] = setColorbarToMinMax()
% opens a .fig then sets its color bar to min & max

[fName, pName, ~] = uigetfile({'*.fig'}, 'please select a .fig file', 'C:\Users\marwan.muhammad\Desktop');
openfig([pName fName]);

h = colorbar;
colorbar('XTick', [h.Limits(1) h.Limits(2)] ,'XTickLabel', {'min','max'});