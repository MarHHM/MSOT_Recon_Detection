function recon_wMB(scan_path, nonNeg, LVc, RECON_ALL, shift)
%%% try Luis weighting (wBP - wBP with apriori - wMB)
% LVc: false -> original wMB (Luis PMB 2013)

%% VIP
% - make sure that folder "masks" exist in the dataset path
% - input range (zpos, rep, wl)
% - for weighted methods+apriori --> choose w higher than 1
% - for each dataset, try different "w" till u find the best & save it
% - nPairs --> if too large, can fill memory!!


if nargin == 4
    shift = 0;
end
%% PATHS & PARAMS

%%% PARAMS
% scan_path = "S:\PA_DATASETS\MSOT 256\Qutaiba phantom\Scan_8";
disp("Loading conv MB recon structure..")
parent_convMB_recon = load(scan_path + "\recons\MB_Tik - nonNeg_" + int2str(nonNeg) + "-zpos_1-reps_1-wls_11-reconRes_7.5e-05-imW_0.03.mat");      % load an aribtrary recon to test the masks
% parent_convMB_recon = load([scan_path "\recons\MB_Tik - nonNeg_" int2str(nonNeg) "-zpos_8-reps_1-wls_13-reconRes_7.5e-05-imW_0.03.mat"]);      % load an aribtrary recon to test the masks
disp("done.")

SAVE_RECON = true; 
% SIGS_LOADED_IN_RAM = 1;               % if sigMat_truncated (from the convMB_tik) is already loaded in the RAM
LOAD_A_DIMS = 1;                      % if the area containing the whole target has been delineated before
LOAD_B_DIMS = 0;                      % if the area containing the reflector has been delineated before

%%% weighted recon params
LUIS_WEIGHT_METHOD = "wMB";    %    "wMB"     "wMB+apriori"
w = 1;                                   % weighting parameter (heuristic (change for optimal choice with each dataset) - use higher value when acoustic heterogenities (i.e. reflectors or scatterers) are more) (typically 1 - more than 2 begins to distort) (prev called WP) - 0->should produce like conventional MB (MB_Tik)
t_undercompnsated_ratio = 0.5;     % further part of the record length (ratio of the total record length) (typically 0.5) (choosing the latter half of the center of rotation is the same idea by Wang in thier Half-time approach)
% compnsFactor = 0.75;       % should be < 1 (1 deactivates the LV correction (i.e. is the standard Luis algo where degradation happens at upper part)) (typically 0.75)
% LIMITED_VIEW_CORRCTN = "spatio-temporal weighting";     % "spatio-temporal weighting"   "pixel reverse-weighting in Regularization"    "pixel reverse-normalization after recon" % note: for "pixel reverse-weighting in Regularization", deactivate the compGradient part (use the old    "P_d(j,i) = compnsFactor*w;")

%%% wMB+apriori params (as in the Luis APL 2011 paper (wBP+apriori))
P_r.nHist = 40;                 % originially in Luis code 40
P_r.nPairs = parent_convMB_recon.n^2;              % originially in Luis code 500 (used also 50e3 in his paper) --> max theoritical should be "(n^2)^2" if both A & B cover the whole image

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datainfo = parent_convMB_recon.datainfo;
im_w = parent_convMB_recon.im_w;
t = parent_convMB_recon.t;
angle_sensor = parent_convMB_recon.angle_sensor;
ts = parent_convMB_recon.ts;
sigMat_truncated = parent_convMB_recon.sigMat_truncated;
reconRes = parent_convMB_recon.reconRes;
f_min = parent_convMB_recon.f_min;
f_max = parent_convMB_recon.f_max;
t_res = parent_convMB_recon.t_res;
MB_regu = parent_convMB_recon.MB_regu;

% n_det = datainfo.HWDesc.n_det;
% r_sensor = datainfo.HWDesc.Radius;
n_det = 256;
r_sensor = 0.0405;
T = parent_convMB_recon.datainfo.AverageTemperature;
c = 12 + round(1.402385 * 1e3 + 5.038813 * T - 5.799136 * 1e-2 * T^2 + 3.287156 * 1e-4 * T^3 - 1.398845 * 1e-6 * T^4 + 2.787860 * 1e-9 * T^5 );
% [t, ts] = formInterpolationVec(datainfo, n, t_res, im_w, c);

A_mat = load(parent_convMB_recon.A_matPath).A_mat;
run_idx = 1;
if RECON_ALL
        disp("The full dataset will be reconstructed!");
        zpos_begin = 1;       zpos_end = length(parent_convMB_recon.datainfo.ZPositions);
        rep_begin = 1;        rep_end = parent_convMB_recon.datainfo.RepNum;
        wl_begin = 1;         wl_end = length(parent_convMB_recon.datainfo.Wavelengths);
    else
        zpos_begin = 1;                            % 5
        zpos_end = 1;                              % 5
        rep_begin = 1;
        rep_end = 1;
        wl_begin = 1;                               % 6 is 750nm for Hong phantom (scan_69 -> AlexaFl_750)
        wl_end = 1;
end
    

n = floor(im_w/parent_convMB_recon.reconRes);   % reconstructed im size (pixels)
disp("--recon method: " + LUIS_WEIGHT_METHOD + " with w = " + num2str(w) + " - LVc = " + num2str(LVc) + "--");

% if ~SIGS_LOADED_IN_RAM
%     [datainfo_recon, Recon_MB] = loadRecon_iThera(sigMat_pathName);
% end
imToChooseMask = parent_convMB_recon.Recon_MB(:, :, 1, 1, 1, 1);   % arbitrary slice shown to draw (or check) the mask

% load mask containing the whole target (all optical absorpers & acoustic reflectors - A in the 2011 paper)
[xc_A, yc_A, Rc_A, logicalMask_A] = checkMask( LOAD_A_DIMS, "A", scan_path, imToChooseMask, im_w );
if strcmp(LUIS_WEIGHT_METHOD, "wMB+apriori")
    [~, ~, ~, logicalMask_B] = checkMask( LOAD_B_DIMS, "B", scan_path, imToChooseMask, im_w ); % also load mask containing acoustic reflectors (B in the 2001 APL apriori paper)
end

%%% TEST - draw input slice & overlay masks A & B to check them
mask_wholeTarget = rsmak( "circle", Rc_A*n/im_w, [(n/2+xc_A*n/im_w) (n/2+yc_A*n/im_w)] );
figure("name", "just to confirm that masks are correct!!");
    imagesc( imToChooseMask ), title( "conv MB - arbitrary slice" ), colormap( bone ), colorbar, axis image off;
    hold on, fnplt(mask_wholeTarget), axis image off, hold off;
%     print(h);       % forces the figure to be drawn, but sends it to the printer as well

% mask_reflector = rsmak( "circle", Rc_B*n/im_w, [(n/2+xc_B*n/im_w) (n/2+yc_B*n/im_w)] );
% hold on; fnplt(mask_reflector), axis image off;

%% extract input frame & do weighted recon
ReconW = zeros(n, n, run_idx, zpos_end-zpos_begin+1, rep_end-rep_begin+1, wl_end-wl_begin+1);
reconItr = 0;
totNumRecons = size(ReconW,3)*size(ReconW,4)*size(ReconW,5)*size(ReconW,6);

disp("Reconstruction of " + int2str(totNumRecons) + " slices has begun..");
for zpos_idx = 1 : zpos_end-zpos_begin+1
    for rep_idx = 1 : rep_end-rep_begin+1
        for wl_idx = 1 : wl_end-wl_begin+1
            tic;
            
            sigMat_current = sigMat_truncated(:, :, 1, zpos_idx, rep_idx, wl_idx);
            
            switch LUIS_WEIGHT_METHOD
                case "wMB"
%                     Dxy = im_w/(n-1); % sampling distance in x and y
%                     x = (-1)*(n/2-0.5)*Dxy:Dxy:(n/2-0.5)*Dxy; % position of pixels in the x direction
%                     y = (-1)*(n/2-0.5)*Dxy:Dxy:(n/2-0.5)*Dxy; % position of pixels in the y direction
%                     [X,Y] = meshgrid(x,y); % grid of points in x and y
%                     nn = n*n;
%                     lsqr_iter = 50;
                    
                    if LVc
                        %%% accessing specific detector signals causing probs
                        j_c = floor( t_undercompnsated_ratio * length(t) );
                        
                        % see "angles illustration.jpg" 
                        theta_start = -pi/4;            % start of detection ring (counting clock-wise !)
                        total_detection_ring_coverage = 3*pi/2;         % inVision 256 coverage is 270 deg
                        theta_end = theta_start + total_detection_ring_coverage;     % end of detection ring
                        theta_step = (theta_end - theta_start)/(n_det);
                        
                        theta_1 = pi/2;             % start of detection part affected by Luis weighting method ("partial D_u" in my paper)
                        theta_2 = pi;               % end of detection part affected by Luis weighting method
                        det_idx_at_theta_1 = floor(theta_1/theta_step);
                        det_idx_at_theta_2 = floor(theta_2/theta_step);
                        
                        %             % detIdx = 256;
                        %             %         detSigRangeCausingProb = sigMat_truncated(:,det_idx_at_theta_1:det_idx_at_theta_2,1,1,1,1);
                        %             %         figure, plot(detSigRangeCausingProb);
                        %
                        %             %%% form spatio-temporal variant weight factor w (increasing weight at t_ij making undercompensation at pixels at upper part of im)
                        %             %           if strcmp( LIMITED_VIEW_CORRCTN, "pixel reverse-weighting in Regularization") || strcmp(LIMITED_VIEW_CORRCTN, "pixel reverse-normalization after recon" )
                        %             %             compnsFactor = 1;       % deactivates this correction
                        %             %           end
                        %             wSpatiotemp_mat = w*ones(length(t),n_det);
                        %             for i = 1:n_det
                        %               for j = 1:length(t)
                        %                 % if at underweighted zone
                        %                 if (i >= det_idx_at_theta_1) && (i <= det_idx_at_theta_2) && (j > floor(t_undercompnsated_ratio*size(wSpatiotemp_mat,1)))
                        %                   %                     P_d(j,i) = compnsFactor*w;        % lower the penalization for bigger coverage at these transducers
                        %                   compGradient = (length(t)-j)/(length(t)*(1-t_undercompnsated_ratio));   % from 1 to 0 at the culprit part of the signal (starting at "t_undercompnsated_ratio*size(P_d,1)")
                        %                   wSpatiotemp_mat(j,i) = compnsFactor*w*compGradient;
                        %                 end
                        %               end
                        %             end
                        %             %%% TEST
                        %             figure, imagesc(wSpatiotemp_mat), title("wSpatiotemp_{mat}"), xlabel("detectors"), ylabel("t samples");
                        %             figure, plot(wSpatiotemp_mat(:,128));
                    end
                    
                    
                    %%% form weighting matrix P_d(t_ij) as a col-major vector (probablity of direct propagation)(to be multiplied by the model matrix)
                    P_d = zeros(1,n_det*length(t));
                    for i = 1:n_det                                        % for each transducer i
                        theta = angle_sensor(i);                             % angle of the transducer
%                         x_sensor = r_sensor*cos(theta);                     % horizontal position of the transducer
%                         y_sensor = r_sensor*sin(theta);                     % vertical position of the transducer
                        
                        % calculate the part of A covered by the current transducer for different time instants (A_ij for all t_j values for the current trans i)
                        A_ij = area_covered_detector_circle_distance( c, t, xc_A, yc_A, Rc_A, r_sensor, theta );
                        
                        if LVc && (i >= det_idx_at_theta_1+shift) && (i <= det_idx_at_theta_2-shift)  % access detectors in red part (partial(D_u) in fig. 2a) & replace the latter part of A_ij(j) by the const A_ij(j_c) (i.e. fig. 2c)
                            for j = 1:length(t)
                                if j > j_c
                                    A_ij(j) = A_ij(j_c);
                                end
                            end
                            %%% TEST
                            %                 if i == 128
                            %                   figure, plot(t, A_ij), title("A_{ij} for a transducer in the red part after fixing A_{ij}(t_{j}) after j_c"), ylim([0 pi*Rc_A^2]);
                            %                 end
                        end
                        
                        AreaNorm = A_ij / (pi*Rc_A^2);      % normalize A_ij to the area of the target (circle)
                        P_d(((i-1)*length(t)+1):(i*length(t))) = max(0,1-w*AreaNorm);   % eq 6 & 7 in paper (org Luis weighting)
                        
                    end
                    
                    %%% TEST
                    %           P_d_mat = reshape(P_d, [length(t) n_det]);
                    %           figure, plot(t, P_d_mat(:,128)), ylim([0 1]), title("P_{d}(t_{ij}) for a transducer in the red part after fixing A_{ij}(t_{j}) after j_c");
                    %           figure, plot(t, P_d_mat(:,end)), ylim([0 1]), title("P_{d}(t_{ij}) for last transducer");
                    
                    pos = linspace( 1, n_det*length(t), n_det*length(t) );
                    W = sparse( pos, pos, P_d, n_det*length(t), n_det*length(t) );   % just putting P_d as a diagonal matrix
                    
                    %%% weighting A_mat & b_vec
%                     disp("weighting A_mat & b_vec");
                    A_matW = W*A_mat;
                    
                    b_vec_In = prepareMeasuredPressureVec( sigMat_current, ts, t, n_det );
                    b_vecW = W*b_vec_In;
                    
                    %%% TEST
                    %           goodDet = 1;
                    %           badDet = 128;
                    % %           figure, plot( t, b_vec_In( (goodDet-1)*length(t) + 1 : goodDet*length(t) ) ), title("good det - no weighting"), ylim([-3 4]);
                    % %           figure, plot( t, b_vecW( (goodDet-1)*length(t) + 1 : goodDet*length(t) ) ), title("good det - weigthing + LVc"), ylim([-3 4]);
                    %           figure, plot( t, b_vec_In( (badDet-1)*length(t) + 1 : badDet*length(t) ) ), title("bad det - no weighting"), ylim([-3 4]);
                    %           figure, plot( t, b_vecW( (badDet-1)*length(t) + 1 : badDet*length(t) ) ), title("bad det - weigthing ORG"), ylim([-3 4]);
                    
                    
                    %%% do recon
                    ReconW(:, :, 1, zpos_idx, rep_idx, wl_idx) = reconstruction( A_matW, b_vecW, n, parent_convMB_recon.RECON_METHOD, parent_convMB_recon.MB_regu, nonNeg );
                    
                    %           switch LIMITED_VIEW_CORRCTN
                    %             case "spatio-temporal weighting"
                    %               %%% Inversion
                    %               % ReconW = reconstruct_MBLuisWeighted(A_matW, b_vecW, w_vec_pixels_diag, n, reguForLuisMthd, 0);
                    %               ReconW(:, :, 1, zpos_idx, rep_idx, wl_idx) = reconstruction( A_matW, b_vecW, n, RECON_METHOD, MB_regu, IMPOSE_NONNEGATIVITY );
                    %
                    %             case "pixel reverse-weighting in Regularization"      %%% [MHH] pixel reverse weightining for limited view case
                    %               w_vec_mat = reshape(P_d", length(t), n_det);     % reshape to matrix to input to BP algo
                    %               % VIP: don"t use "full" here as we don"t need the derivatives for the pixel weighs reconstruction
                    %               [w_vec_pixels, X, Y] = backproject_luis(w_vec_mat, n, Radius, angle_sensor, c, "direct", t, im_w, 0, 0);
                    %               figure, imagesc(w_vec_pixels), title("reconstructed pixel weights (using BP)"), colormap(bone), axis image off;
                    %               w_vec_pixels_colMajor = reshape(w_vec_pixels, n*n,1);     % reshape back to vector
                    %               w_vec_pixels_colMajor = w_vec_pixels_colMajor + 1;      % to avoid division by zero in the next command
                    %               w_vec_pixels_colMajor = 1./w_vec_pixels_colMajor;         % Q: is this correct?
                    %               w_vec_pixels_diag = diag(w_vec_pixels_colMajor);
                    %
                    %               %%% Inversion
                    %               ReconW = reconstruct_MBLuisWeighted(A_matW, b_vecW, w_vec_pixels_diag, n, MB_regu, IMPOSE_NONNEGATIVITY);
                    %               % ReconW = reconstruction(A_matW, b_vecW, n, RECON_METHOD, reguForLuisMthd, IMPOSE_NONNEGATIVITY);
                    %             case "pixel reverse-normalization after recon"
                    %               disp("!! conceptually wrong as the signal in the distorted upper part is already lost after the recon!!");
                    %               %                 w_vec_mat = reshape(P_d", length(t), n_det);     % reshape to matrix to input to BP algo
                    %               %                 [w_vec_pixels, X, Y] = backproject_luis(w_vec_mat, n, Radius, angle_sensor, c, "full", t, SamplingFrequency, im_w, 0, 0);
                    %               %                 w_vec_pixels = w_vec_pixels + 1;      % to avoid div by zero
                    %               %
                    %               %                 %%% Inversion then normalization by pixel weights
                    %               %                 ReconW = reconstruction(A_matW, b_vecW, n, RECON_METHOD, reguForLuisMthd, IMPOSE_NONNEGATIVITY);
                    %               %                 ReconW = ReconW./w_vec_pixels;
                    %           end
                case "wMB+apriori"
                    wMB_apriori;
                case "wBP"
                    image_select = "derivative";    % type of back-projection in Luis wBP
                    [~,ReconWBP,~,~] = backproject_luis_prc_circle( sigMat_current, n, r_sensor,angle_sensor,c,image_select,ts,im_w,0,0,xc_A,yc_A,Rc_A,w,1);
                    figure; imagesc(ReconWBP); title(["weighted BP (\omega = " num2str(w) " )"]), colormap(bone); axis image off;
                case "wBP+apriori"
                    image_select = "derivative";    % type of back-projection in Luis wBP
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
                    %             circleAbs = rsmak("circle", R_abs*n/im_w, [(n/2+x_abs*n/im_w) (n/2+y_abs*n/im_w)]);
                    %             figure, imagesc(Recon_MB(:, :, 1, zpos_idx, rep_idx, wl_idx)), title("conv BP"), colormap(bone), axis image off;
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
                    %             circleRef = rsmak("circle", Rc_ref*n/im_w, [(n/2+xc_ref*n/im_w) (n/2+yc_ref*n/im_w)]);
                    %             hold on; fnplt(circleRef);
                    
                    % MAYBE IT"S NOT USEFUL!!
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
                    [ReconW, ~, ~] = backproject_luis_ref_priors( sigMat_current, n, r_sensor, angle_sensor, c, image_select, ts,...
                        im_w, histogram, WT, w );
                    %           figure; imagesc(ReconW); title(["weighted BP with apriori (\omega = " num2str(w) " - n points = " num2str(n_points) ")"]), colormap(bone); axis image off;
                    figure; imagesc(ReconW); title(["weighted BP with apriori (\omega = " num2str(w) " - nHist = " num2str(P_r.nHist) ")"]), colormap(bone); axis image off;
            end
            
            reconItr = reconItr+1;
            disp("recon " + num2str(reconItr) + " of " + num2str(totNumRecons) + " done..");
            toc;
        end
    end
end

%%% draw result
% [tickLocs, tickLabels] = getImageTicks(n, im_w);
figure,...
    imagesc(ReconW(:, :, 1, 1, 1, 1)); colormap(bone), title(LUIS_WEIGHT_METHOD+" (\omega = "+num2str(w)+") - LVc = "+num2str(LVc), "Interpreter", "tex"), colorbar, axis image off;    % OR is f(t_{ij})
% set(gca, "YTick", tickLocs, "XTick", tickLocs, "XTickLabel", tickLabels, "YTickLabel", tickLabels);

if SAVE_RECON
    if ~exist((scan_path + "\recons"), "dir")
        mkdir((scan_path + "\recons"))
    end
    reconStruct__fname = ( scan_path + "\recons\"+ LUIS_WEIGHT_METHOD+ " - nonNeg_" + num2str(nonNeg) + "-zpos_" +...
                           num2str(zpos_end-zpos_begin+1) + "-reps_" + num2str(rep_end-rep_begin+1) + "-wls_" + num2str(wl_end-wl_begin+1) +...
                           "-reconRes_" + num2str(reconRes) + "-imW_" + num2str(im_w)  + "-w_" + num2str(w) +"-LVc_" + num2str(LVc) + "shift" + num2str(shift) + ".mat"...
                          );
    if exist(reconStruct__fname, "file")                % to avoid overwriting an existing recon
       [filepath, name, ext] = fileparts(reconStruct__fname);
        reconStruct__fname = filepath+"/"+name+"__"+ext;
    end
%     note = "sigMat_truncated attached with this recon is laser-energy corrected & filtered";
    save(reconStruct__fname, "ReconW", "sigMat_truncated", "datainfo", "im_w", "reconRes", "f_min", "f_max", "t_res", "MB_regu",...
        "zpos_begin", "zpos_end", "wl_begin", "wl_end", "rep_begin", "rep_end", "nonNeg", "LVc", "P_r", "-v7.3" );
    disp("-> Reconstruction saved to disk.");
    
end

%% sharpness calculation (Hong code) --> doesn"t work well with my algo!!
% % Recon_temppp = Recon_MB(imLims.y0:imLims.yf,imLims.x0:imLims.xf, 1, zpos_idx, rep_idx, wl_idx);
% Recon_temppp = ReconW(imLims.y0:imLims.yf,imLims.x0:imLims.xf);
% [Gx, Gy]=gradient(Recon_temppp);
% S1=sqrt(Gx.*Gx+Gy.*Gy);
% sharpness1=sum(sum(S1));
% sharpness11=sum(sum(S1))./(numel(Gx));
