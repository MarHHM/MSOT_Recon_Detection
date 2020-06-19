%%% try Luis weighting (wBP - wBP with apriori - wMB) - run after 'myRecon_single.m'

%% VIP
% - make sure that folder 'masks' exist in the dataset path
% - input range (slc, rep, wl)
% - for weighted methods+apriori --> choose w higher than 1 
% - for each dataset, try different 'w' till u find the best & save it
% - nPairs --> if too large, can fill memory!!

%% PATHS & PARAMS

%%% PARAMS
SIGS_LOADED_IN_RAM = 1;               % if sigMat_truncated (from the convMB_tik) is already loaded in the RAM
LOAD_A_DIMS = 1;                      % if the area containing the target has been delineated before
LOAD_B_DIMS = 1;                     % if the area containing the reflector has been delineated before
slcToChooseMask = Recon_MB(:, :, 1, 1, 1, 4);

%%% weighted recon params
LUIS_WEIGHT_METHOD = 'wMB';    % 'wBP+apriori'     'wMB'     'wMB+apriori'
w = 1;                                 % weighting parameter (heuristic (change for optimal choice with each dataset) - use higher value when acoustic heterogenities (i.e. reflectors or scatterers) are more) (typically 1 - more than 2 begins to distort) (prev called WP)
%%% wMB+apriori params
% P_r.nHist = 40;                 % originially in Luis code 40
% P_r.nPairs = n^2;              % originially in Luis code 500 (used also 50e3 in his paper) --> max theoritical should be '(n^2)^2' if both A & B cover the whole image
%%% wMB+LVc params
% LIMITED_VIEW_CORRCTN = 'spatio-temporal weighting';     % 'spatio-temporal weighting'   'pixel reverse-weighting in Regularization'    'pixel reverse-normalization after recon'
CORRCT_FOR_LV = 1;
% note: for "pixel reverse-weighting in Regularization", deactivate the compGradient part (use the old    "P_d(j,i) = compnsFactor*w;")
t_undercompnsated_ratio = 0.5;     % further part of the record length (ratio of the total record length) (typically 0.5)
compnsFactor = 0.75;       % should be < 1 (1 deactivates the LV correction (i.e. is the standard Luis algo where degradation happens at upper part)) (typically 0.75)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = floor(im_w/reconRes);   % reconstructed im size (pixels)
disp( ['--recon method: ' LUIS_WEIGHT_METHOD ' with w = ' num2str(w) '--'] );

nDet = datainfo.HWDesc.NumDetectors;

if ~SIGS_LOADED_IN_RAM
  [datainfo_recon, Recon_MB] = loadRecon_iThera(sigMat_pathName);
end

% mask containing all optical absorpers & acoustic reflectors (A in the paper)
[xc_A, yc_A, Rc_A, logicalMask_A] = checkMask( LOAD_A_DIMS, 'A', newDataPath, slcToChooseMask, im_w );
mask_absRef = rsmak( 'circle', Rc_A*n/im_w, [(n/2+xc_A*n/im_w) (n/2+yc_A*n/im_w)] );       % just for drawing

% mask containing acoustic reflectors (B in the paper)
[xc_B, yc_B, Rc_B, logicalMask_B] = checkMask( LOAD_B_DIMS, 'B', newDataPath, slcToChooseMask, im_w );
mask_ref = rsmak( 'circle', Rc_B*n/im_w, [(n/2+xc_B*n/im_w) (n/2+yc_B*n/im_w)] );


%%% draw input slice & overlay masks A & B to confirm right dims
% % [tickLocs, tickLabels] = getImageTicks(n, im_w);
% figure('name', 'to confirm that masks are correct!!'),...
%   imagesc( slcToChooseMask ), title( ['conv ' RECON_METHOD ' - arbitrary slice'], 'interpreter', 'none' ), colormap( bone ), axis image off;
% % set(gca, 'YTick', tickLocs, 'XTick', tickLocs, 'XTickLabel', tickLabels, 'YTickLabel', tickLabels);
% hold on; fnplt( mask_absRef ), axis image off;
% hold on; fnplt( mask_ref ), axis image off;

%% extract input frame & do weighted recon
ReconW = zeros(n,n,1,slc_end,rep_end,wl_end);
reconItr = 0;
totNumRecons = size(sigMat_truncated,3)*size(sigMat_truncated,4)*size(sigMat_truncated,5)*size(sigMat_truncated,6);


for slc_idx = 1 : slc_end-slc_begin+1
  for rep_idx = 1 : rep_end-rep_begin+1
    for wl_idx = 1 : wl_end-wl_begin+1
        
      tic;
      
      sigMat_current = sigMat_truncated(:, :, 1, slc_idx, rep_idx, wl_idx);
      r_sensor = datainfo.HWDesc.Radius;
      
      switch LUIS_WEIGHT_METHOD
        case 'wBP'
          image_select = 'derivative';    % type of back-projection in Luis wBP
          [ReconBP,ReconWBP,X,Y] = backproject_luis_prc_circle( sigMat_current, n, r_sensor,angle_sensor,c,image_select,ts,im_w,0,0,xc_A,yc_A,Rc_A,w,1);
          figure; imagesc(ReconWBP); title(['weighted BP (\omega = ' num2str(w) ' )']), colormap(bone); axis image off;
        case 'wBP+apriori'
          image_select = 'derivative';    % type of back-projection in Luis wBP
          %             %%% calculate masks for absorption and reflection
          %             % note: X & Y --> mesh coming from the conv BP
          %             x_abs = xc_A;
          %             y_abs = yc_A;
          %             R_abs = Rc_A;
          %             Dxy = im_w/(n-1); % spatial sampling distance of the reconstructed image
          %             x = [(-1)*(n/2-0.5)*Dxy:Dxy:(n/2-0.5)*Dxy]; x = x+xc; % position of pixels of the reconstructed image in the x direction
          %             y = [(-1)*(n/2-0.5)*Dxy:Dxy:(n/2-0.5)*Dxy]; y = y+yc; % position of pixels of the reconstructed image in the y direction
          %             [X,Y] = meshgrid(x,y); % position of the pixels of the reconstructed image (mesh)
          %             mask_abs = zeros(size(X,1),size(X,2));
          %             mask_abs(((X-x_abs).^2+(Y-y_abs).^2)<R_abs^2) = 1;
          %             circleAbs = rsmak('circle', R_abs*n/im_w, [(n/2+x_abs*n/im_w) (n/2+y_abs*n/im_w)]);
          %             figure, imagesc(Recon_MB(:, :, 1, slc_idx, rep_idx, wl_idx)), title('conv BP'), colormap(bone), axis image off;
          %             hold on; fnplt(circleAbs);
          %
          %             %%% reflector dims (ink up)
          %             %         x_ref = 0e-3;
          %             %         y_ref = 7.2e-3;
          %             %         R_ref = 2.1e-3;     % physical r of the straw
          %             %%% reflector dims (straw up)
          %             xc_ref = 2.8e-3;
          %             yc_ref = -7e-3;
          %             Rc_ref = 2.1e-3;
          %             %%% reflector dims (straw up)
          %             %         x_ref = 7.1e-3;
          %             %         y_ref = -1.5e-3;
          %             %         R_ref = 2.1e-3;
          %
          %             mask_ref = zeros(size(X,1),size(X,2));
          %             mask_ref(((X-xc_ref).^2+(Y-yc_ref).^2)<Rc_ref^2) = 1;
          %             circleRef = rsmak('circle', Rc_ref*n/im_w, [(n/2+xc_ref*n/im_w) (n/2+yc_ref*n/im_w)]);
          %             hold on; fnplt(circleRef);
          
          % MAYBE IT'S NOT USEFUL!!
%           [Recon, X, Y] = backproject_luis( sigMat_current, n, r_sensor, angle_sensor, c0, image_select, ts, image_width, 0, 0 );
          
          %%% calculate histogram
          T = datainfo.AverageTemperature;
          c = 12 + round(1.402385 * 1e3 + 5.038813 * T - 5.799136 * 1e-2 * T^2 + 3.287156 * 1e-4 * T^3 - 1.398845 * 1e-6 * T^4 + 2.787860 * 1e-9 * T^5 );
          Dxy = im_w/(n-1); % spatial sampling distance of the reconstructed image
          x = -1*(n/2-0.5)*Dxy : Dxy : (n/2-0.5)*Dxy; x = x+xc_A; % position of pixels of the reconstructed image in the x direction
          y = -1*(n/2-0.5)*Dxy : Dxy : (n/2-0.5)*Dxy; y = y+yc_A; % position of pixels of the reconstructed image in the y direction
          [X,Y] = meshgrid(x,y); % position of the pixels of the reconstructed image (mesh)
          %             histogram = calculate_time_histograms_refs(mask_abs, mask_ref, X, Y, ts, n_points, r_sensor, angle_sensor, c, hist.nHist);
          histogram = calculate_time_histograms_refs( logicalMask_A, logicalMask_B, X, Y, ts, P_r.nPairs, r_sensor, angle_sensor, c, P_r.nHist );
          %         figure, imagesc(histogram);
          
          % weighted reconstruction
          WT = 1;
          [ReconW, X, Y] = backproject_luis_ref_priors( sigMat_current, n, r_sensor, angle_sensor, c, image_select, ts,...
                                                        im_w, histogram, WT, w );
          %           figure; imagesc(ReconW); title(['weighted BP with apriori (\omega = ' num2str(w) ' - n points = ' num2str(n_points) ')']), colormap(bone); axis image off;
          figure; imagesc(ReconW); title(['weighted BP with apriori (\omega = ' num2str(w) ' - nHist = ' num2str(P_r.nHist) ')']), colormap(bone); axis image off;
        
        case 'wMB'
          Dxy = im_w/(n-1); % sampling distance in x and y
          x = (-1)*(n/2-0.5)*Dxy:Dxy:(n/2-0.5)*Dxy; % position of pixels in the x direction
          y = (-1)*(n/2-0.5)*Dxy:Dxy:(n/2-0.5)*Dxy; % position of pixels in the y direction
          [X,Y] = meshgrid(x,y); % grid of points in x and y
          nn = n*n;
          lsqr_iter = 50;
          
          
          if CORRCT_FOR_LV
            %%% accessing specific detector signals causing probs
            detAngleBegin = pi/4;
            detAngleEnd = 3*pi/4;
            detIdxAtAngleBegin = floor((detAngleBegin-datainfo.HWDesc.StartAngle)/datainfo.HWDesc.StepAngle) + 1;
            detIdxAtAngleEnd = floor((detAngleEnd-datainfo.HWDesc.StartAngle)/datainfo.HWDesc.StepAngle) + 1;
            % detIdx = 256;
            %         detSigRangeCausingProb = sigMat_truncated(:,detIdxAtAngleBegin:detIdxAtAngleEnd,1,1,1,1);
            %         figure, plot(detSigRangeCausingProb);
            
            %%% form spatio-temporal variant weight factor w (increasing weight at t_ij making undercompensation at pixels at upper part of im)
            %           if strcmp( LIMITED_VIEW_CORRCTN, 'pixel reverse-weighting in Regularization') || strcmp(LIMITED_VIEW_CORRCTN, 'pixel reverse-normalization after recon' )
            %             compnsFactor = 1;       % deactivates this correction
            %           end
            wSpatiotemp_mat = w*ones(length(t),nDet);
            for i = 1:nDet
              for j = 1:length(t)
                % if at underweighted zone
                if (i >= detIdxAtAngleBegin) && (i <= detIdxAtAngleEnd) && (j > floor(t_undercompnsated_ratio*size(wSpatiotemp_mat,1)))
                  %                     P_d(j,i) = compnsFactor*w;        % lower the penalization for bigger coverage at these transducers
                  compGradient = (length(t)-j)/(length(t)*(1-t_undercompnsated_ratio));   % from 1 to 0 at the culprit part of the signal (starting at "t_undercompnsated_ratio*size(P_d,1)")
                  wSpatiotemp_mat(j,i) = compnsFactor*w*compGradient;
                end
              end
            end
            %%% TEST
            figure, imagesc(wSpatiotemp_mat), title('wSpatiotemp_{mat}'), xlabel('detectors'), ylabel('t samples');
            figure, plot(wSpatiotemp_mat(:,128));
          end
          
          
          %%% form weighting matrix P_d(t_ij) as a col-major vector (probablity of direct propagation)(to be multiplied by the model matrix)
          P_d = zeros(1,nDet*length(t));
          for i = 1:nDet                                        % for each transducer i
            theta = angle_sensor(i);                             % angle of the transducer
            x_sensor = r_sensor*cos(theta);                     % horizontal position of the transducer
            y_sensor = r_sensor*sin(theta);                     % vertical position of the transducer
            
            % calculate the part of A covered by the current transducer for different time instants (A_ij for all t_j values for the current trans i)
            A_ij = area_covered_detector_circle_distance( c, t, xc_A, yc_A, Rc_A, r_sensor, theta );
            AreaNorm = A_ij / (pi*Rc_A^2);      % normalize A_ij to the area of the target (circle)

            if CORRCT_FOR_LV
              P_d(((i-1)*length(t)+1):(i*length(t))) = max(0,1-(wSpatiotemp_mat(:,i)'.*AreaNorm));   % eq 6 & 7 in paper (time instants corresponding to the current detector)
            else
              P_d(((i-1)*length(t)+1):(i*length(t))) = max(0,1-w*AreaNorm);   % eq 6 & 7 in paper (org Luis weighting)
            end

          end
          %%% TEST
          figure, plot(t, A_ij), title('A_{ij} for last transducer');
          
          pos = linspace( 1, nDet*length(t), nDet*length(t) );
          W = sparse( pos, pos, P_d, nDet*length(t), nDet*length(t) );   % just putting P_d as a diagonal matrix
          
          %%% TEST
          P_d_mat = reshape(P_d, [length(t) nDet]);
          figure, imagesc(P_d_mat);
          
          
          %%% weighting A_mat & b_vec
          A_matW = W*A_mat;
          
          b_vec_In = prepareMeasuredPressureVec( sigMat_current, ts, t, nDet );
          b_vecW = W*b_vec_In;
          
          %%% do recon
          ReconW(:, :, 1, slc_idx, rep_idx, wl_idx) = reconstruction( A_matW, b_vecW, n, RECON_METHOD, MB_regu, IMPOSE_NONNEGATIVITY );
          
%           switch LIMITED_VIEW_CORRCTN
%             case 'spatio-temporal weighting'
%               %%% Inversion
%               % ReconW = reconstruct_MBLuisWeighted(A_matW, b_vecW, w_vec_pixels_diag, n, reguForLuisMthd, 0);
%               ReconW(:, :, 1, slc_idx, rep_idx, wl_idx) = reconstruction( A_matW, b_vecW, n, RECON_METHOD, MB_regu, IMPOSE_NONNEGATIVITY );
%               
%             case 'pixel reverse-weighting in Regularization'      %%% [MHH] pixel reverse weightining for limited view case
%               w_vec_mat = reshape(P_d', length(t), nDet);     % reshape to matrix to input to BP algo
%               % VIP: don't use "full" here as we don't need the derivatives for the pixel weighs reconstruction
%               [w_vec_pixels, X, Y] = backproject_luis(w_vec_mat, n, datainfo.HWDesc.Radius, angle_sensor, c, 'direct', t, im_w, 0, 0);
%               figure, imagesc(w_vec_pixels), title('reconstructed pixel weights (using BP)'), colormap(bone), axis image off;
%               w_vec_pixels_colMajor = reshape(w_vec_pixels, n*n,1);     % reshape back to vector
%               w_vec_pixels_colMajor = w_vec_pixels_colMajor + 1;      % to avoid division by zero in the next command
%               w_vec_pixels_colMajor = 1./w_vec_pixels_colMajor;         % Q: is this correct?
%               w_vec_pixels_diag = diag(w_vec_pixels_colMajor);
%               
%               %%% Inversion
%               ReconW = reconstruct_MBLuisWeighted(A_matW, b_vecW, w_vec_pixels_diag, n, MB_regu, IMPOSE_NONNEGATIVITY);
%               % ReconW = reconstruction(A_matW, b_vecW, n, RECON_METHOD, reguForLuisMthd, IMPOSE_NONNEGATIVITY);
%             case 'pixel reverse-normalization after recon'
%               disp('!! conceptually wrong as the signal in the distorted upper part is already lost after the recon!!');
%               %                 w_vec_mat = reshape(P_d', length(t), nDet);     % reshape to matrix to input to BP algo
%               %                 [w_vec_pixels, X, Y] = backproject_luis(w_vec_mat, n, datainfo.HWDesc.Radius, angle_sensor, c, 'full', t, datainfo.HWDesc.SamplingFrequency, im_w, 0, 0);
%               %                 w_vec_pixels = w_vec_pixels + 1;      % to avoid div by zero
%               %
%               %                 %%% Inversion then normalization by pixel weights
%               %                 ReconW = reconstruction(A_matW, b_vecW, n, RECON_METHOD, reguForLuisMthd, IMPOSE_NONNEGATIVITY);
%               %                 ReconW = ReconW./w_vec_pixels;
%           end
             
        case 'wMB+apriori'
          wMB_apriori;
      end
      
      reconItr = reconItr+1;
      disp(['recon ' num2str(reconItr) ' of ' num2str(totNumRecons) ' done..']);
      toc;
    end
  end
end

%%% draw result
% [tickLocs, tickLabels] = getImageTicks(n, im_w);
figure,...
  imagesc(ReconW(:, :, 1, 1, 1, 4)); colormap(bone), title([LUIS_WEIGHT_METHOD ' (\omega = ' num2str(w) ')']), colorbar, axis image off;    % OR is f(t_{ij})
% set(gca, 'YTick', tickLocs, 'XTick', tickLocs, 'XTickLabel', tickLabels, 'YTickLabel', tickLabels);

if SAVE_RECON
  % don't forget to create the folder "recons" in this path
  note = 'sigMat_truncated attached with this recon is laser-energy corrected & filtered';
  save( strcat( sigMat_pathName, 'recons\', LUIS_WEIGHT_METHOD, ' - nonNeg_', num2str(IMPOSE_NONNEGATIVITY), '-slcs_', num2str(slc_end-slc_begin+1), '-reps_',...
    num2str(rep_end-rep_begin+1), '-wls_', num2str(wl_end-wl_begin+1), '-reconRes_', num2str(reconRes),...
    '-imW_', num2str(im_w) , '-w_', num2str(w), '.mat'), 'ReconW', 'sigMat_truncated', 'note', ...
    'datainfo', 'im_w', 'reconRes', 'f_min', 'f_max', 't_res', 'MB_regu', 'slc_begin', 'slc_end', 'wl_begin',...
    'wl_end', 'rep_begin', 'rep_end', 'IMPOSE_NONNEGATIVITY', 'P_r', '-v7.3' );
  
end

%% sharpness calculation (Hong code) --> don't work well with my algo!!
% % Recon_temppp = Recon_MB(imLims.y0:imLims.yf,imLims.x0:imLims.xf, 1, slc_idx, rep_idx, wl_idx);
% Recon_temppp = ReconW(imLims.y0:imLims.yf,imLims.x0:imLims.xf);
% [Gx, Gy]=gradient(Recon_temppp);
% S1=sqrt(Gx.*Gx+Gy.*Gy);
% sharpness1=sum(sum(S1));
% sharpness11=sum(sum(S1))./(numel(Gx));
