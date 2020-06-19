% prerequisites:
%   - load the weighted recon "ReconW" from "tryLuisWeighting.m"

%% PARAMS
PROCESSING_METHOD = 'adapt HistEq';     % 'non' - 'adapt HistEq'
gamma = 1.3;            % > 1 darker, < 1 brighter
logShift = 1;         % to avoid log(0)
RECON_LOADED_IN_RAM = 1;               % if recon is already loaded in the RAM

IS_BRAIN_IM = 0;        % to flip the brain vertically (according to Saak)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% loading conv MB
if ~RECON_LOADED_IN_RAM
    Recon_MB = load('D:\DATASETS\MSOT 256\Brain (Ivan)\study 2\Scan_65 (Ivan - brain data to test reflections)\recons\reconMB_Tik_nonNeg_0_imSz_250_slcs_1_wls_3_reps_1.mat');
    Recon_MB = Recon_MB.Recon_MB;
end

% I = squeeze(Recon_MB(:,:,1,slcIn,repIn,wlIn));
I = ReconW(:,:);
if IS_BRAIN_IM
    I = imrotate(I, 180);
end

switch PROCESSING_METHOD
    case 'non'
        Im_processed = I;
    case 'gamma correction'
        Im_processed = imadjust(mat2gray(I),[],[],gamma);      % not too much diffrence for the moment
        figure, imshowpair(I,Im_processed,'montage'), colormap('gray'), axis image off;
    case 'DR compression'
        Im_processed = log(mat2gray(I)+logShift);
        figure, imshowpair(I,Im_processed,'montage'), colormap('gray'), axis image off;
        figure, subplot(1,2, 1), imagesc(I), colormap(gray), axis image off;
                subplot(1,2, 2), imagesc(Im_processed), colormap(gray), axis image off;
    case 'adapt HistEq'
        I_norm = mat2gray(I);
        Im_processed = adapthisteq(I_norm);
        % figure, imhist(I_norm);
        % figure, imhist(J);
        % figure, imshowpair(I,J, 'montage');
end

Im_processed = truncateBgd(Im_processed, targetPath);         % remove empty areas of the image
% figure, imagesc(I), colormap(gray), colorbar, axis image off; title('before processing');
figure, imagesc(Im_processed), colormap(gray), colorbar, axis image off;

%% ROI selection
if ~exist('pos','var')
    h = imrect;
    pos = wait(h);          % wait till the rec is double-clicked & get pos
    delete(h);
    
    roi.x0 = floor(pos(1));
    roi.y0 = floor(pos(2));
    roi.xf = floor(pos(1)+pos(3));
    roi.yf = floor(pos(2)+pos(4));
end

rectangle('position', [pos(1) pos(2) pos(3) pos(4)],...
        'EdgeColor', 'y',...
        'LineStyle', '--');
roi.Im = Im_processed(roi.y0:roi.yf,roi.x0:roi.xf);
figure, imagesc(roi.Im), colormap(gray), axis image off;
% DON'T forget to save the figure!!

%% select & plot line
linePosOut = selectAndPlotLine(linePosOut);
% % % choose & plot the image on which you will select the linePlot (e.g. the ROI from the prev section)
% imHandel = figure, imagesc(roi.Im), colormap(gray), axis image off;
% 
% if ~exist('linePos','var')
%     h = imline;
%     linePos = wait(h);
%     delete(h);
%     
%     linePlot.x0 = floor(linePos(1,1));
%     linePlot.y0 = floor(linePos(1,2));
%     linePlot.xf = floor(linePos(2,1));
%     linePlot.yf = floor(linePos(2,2));
% end
%  
% myLine = line([linePos(1,1) linePos(2,1)], [linePos(1,2) linePos(2,2)], ...
%      'Color', 'y',...
%      'LineStyle', '--');
% temp = getimage(imHandel);
% figure;
% if atan(abs(linePlot.yf-linePlot.y0)/abs(linePlot.xf-linePlot.x0)) < pi/4           % detect if line is vertical or horizontal
%     plot([linePlot.x0:linePlot.xf], temp(linePlot.y0,linePlot.x0:linePlot.xf)), xlabel('pixel number'), ylabel('OA intensity (au)');
% else
%     plot([linePlot.y0:linePlot.yf], temp(linePlot.y0:linePlot.yf,linePlot.x0)), xlabel('pixel number'), ylabel('OA intensity (au)');
% end
% 
% linePosXX = selectAndPlotLine(linePosIn);

%% BGD removal
% mask = roipoly;
% figure, imagesc(J.*mask), colormap(gray), axis image off;
