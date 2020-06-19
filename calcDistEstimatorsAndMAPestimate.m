% % VIP
% %     - change 'path_priorDataset' to the path of the dataset from which u want to construct the prior estimators (mean & cov)
% %     - make sure that "sosVec" is the same as for the dataset
% %     - make sure u have folder "recons\MAP recons" in the dataset path & change the naming of the saved results (to reflect changes from last time)
% %     - update "path_datasetIn" to the desired input dataset on which u want to test inversion on
% 
% addReconPaths;
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% prior dist param estimators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%PARAMS
% 
% % load any recon just to read the dims from it
% path_priorDataset = 'E:\Marwan\Datasets_offline\MSOT 256\(2017-02-04) uSpheres for basis cov matrix (JP idea)\depth50mm_res200um\';
% tempRecon = load([path_priorDataset 'recons\recons to build prior dist\BPn\no added noise\reconBackProjection_imSz100_sos1480_rep1to1000.mat']);
% time_res = 2;
% tempRecon = tempRecon.reconStack;
% Zpos_begin = 1;
% Zpos_end = size(tempRecon, 4);
% % Zpos_end = 1;
% slcNum = Zpos_end-Zpos_begin+1;
% wl_begin = 1;
% wl_end = size(tempRecon, 3);
% wlNum = wl_end-wl_begin+1;
% rep_begin = 1;
% rep_end = size(tempRecon, 5);
% repNum = rep_end-rep_begin+1;
% im_ht = size(tempRecon,1);
% im_wd = size(tempRecon,2);
% imSz = im_ht*im_wd;
% 
% E_MODE = 'white standard';      % noise mode
% newNoiseDataPath = 'E:\Marwan\Datasets_offline\MSOT 256\(2016-12-01) central rod with uspheres around\no avg per acqstn\water only\';
% 
% %% concatenate reps for each pixel from diffenret SOS recons (to increase the effective #obs to be = #RVs for the Cov Mat to be pos def)
% sosVec = linspace(1480,1540,11);
% reconStack_imColMajor = zeros(slcNum, wlNum, repNum*length(sosVec), imSz);
% pxlVsRep_currSOS = zeros(repNum,imSz);
% for sosIdx = 1:length(sosVec)
%     sos_current = sosVec(sosIdx);
%     reconStack_currSOS = load([path_priorDataset 'recons\recons to build prior dist\BPn\no added noise\reconBackProjection_imSz100_sos' num2str(sos_current) '_rep1to' num2str(repNum) '.mat' ]);
%     reconStack_currSOS = reconStack_currSOS.reconStack;
%     
%     for slcIdx = 1 : slcNum
%         for wlIdx = 1 : wlNum
%             for repIdx = 1 : repNum
%                 % convert recon images to column-major vecs
%                 pxlVsRep_currSOS(repIdx,:) = reshape(reconStack_currSOS(:,:,wlIdx,slcIdx,repIdx), imSz, 1);
%             end
%             % rearrange reconstack dims to a logical order (descending --> Zpos, wl, rep, im)
%             reconStack_imColMajor(slcIdx,wlIdx,((sosIdx-1)*repIdx)+1:sosIdx*repIdx,:) = pxlVsRep_currSOS;
%         end
%     end
% end
% 
% %% form priors (per wl per Zpos)
% p0mean = zeros(slcNum, wlNum, imSz);
% 
% % take the mean over the "reps" dim
% for slcIdx = 1 : slcNum
%     for wlIdx = 1 : wlNum
%         p0mean(slcIdx,wlIdx,:) = mean(reconStack_imColMajor(slcIdx, wlIdx, :, :), 3);
%     end
% end
% p0mean = squeeze(p0mean);
% 
% % form the Cov matrix (per Zpos, per wl)
% P0_COV_TYPE = 'sample estimator';
% p0cov = zeros(slcNum, wlNum, imSz, imSz);
% for slcIdx = 1 : slcNum
%     for wlIdx = 1 : wlNum
%         disp('calculating cov mat of the prior p0..');
%         tic;
%         switch P0_COV_TYPE
%             case 'white noise prior'
%             case 'Matern prior'
% %                                 % to-do: continue filling these params
% %                                 sig_m = ;
% %                                 v = ;
% %                                 r_i = ;
% %                                 r_j = ;
% %                                 l = ;
% %                 
% %                                 term = sqrt(2*v)*abs(r_i-r_j)/l;
% %                                 for i = 1: size(p0cov,3)
% %                                     for j = 1:size(p0cov,4)
% %                                         p0cov(slcIdx,wlIdx,i,j) = sig_m^2*(2^(1-v)/gamma(v)) * term^v * besselk(v, term);
% %                                     end
% %                                 end
%             case 'sample estimator'
%                 p0cov(slcIdx,wlIdx,:,:) = cov(squeeze(reconStack_imColMajor(slcIdx, wlIdx, :, :)));
%         end
%         disp(['prior cov mat calculation using ' P0_COV_TYPE ' model took ' num2str(toc) ' s']);
%     end
% end
% 
% p0cov = squeeze(p0cov);
% 
% 
% % SUM COV MATRICES OVER SLC DIM THEN SEE HOW TO WEIGHT IT USING INPUT HEURISTIC RECON
% % (Q) LINEARITY prob suspected!!
% p0meanAvg = mean(p0mean, 1);
% p0covAvg = squeeze(mean(p0cov, 1));
% 
% % % calc cholesky decomp & its inv of every cov mat (per slc per wl)
% % p0covInvCholsky = zeros(slcNum,wlNum,imSz,imSz);
% % for slcIdx = 1 : slcNum
% %     for wlIdx = 1 : wlNum
% %         tic;
% %         temp = chol(squeeze(p0cov(slcIdx,wlIdx,:,:)), 'lower');
% %         p0covInvCholsky(slcIdx,wlIdx,:,:) = inv(temp);
% %         toc;
% %     end
% % end
% 
% 
% %% %%%%% noise dist param estimation      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % calc noise dist param estimators (per wl per Zpos)
% 
% % load sigMat to read dims & params
% if ~exist('sigMatNoise_pathName', 'var') || ~strcmp(sigMatNoise_pathName, newNoiseDataPath)
%     [datainfo_noise, sigMat_noise, sigMatNoise_pathName] = loadSigMat_iThera(newNoiseDataPath);
% else
%     disp(['scan "' sigMatNoise_pathName '" is already loaded']);
% end
% 
% Zpos_begin = 1;
% Zpos_end = size(sigMat_noise, 4);
% % Zpos_end = 1;
% slcNum = Zpos_end-Zpos_begin+1;
% wl_begin = 1;
% wl_end = size(sigMat_noise, 3);
% % wl_end = 1;
% wlNum = wl_end-wl_begin+1;
% rep_begin = 1;
% rep_end = size(sigMat_noise, 5);
% repNum = rep_end-rep_begin+1;
% tSamplesPerDetNum = size(sigMat_noise,1);       %% TO BE TRUNCATED!!
% detNum = size(sigMat_noise,2);
% 
% sigMat_noise_truncated = sigMat_noise(:, :, run_idx, Zpos_begin:Zpos_end, rep_begin:rep_end, wl_begin:wl_end);
% % pre-filtration (is not optimized for large datasets (e.g. high rep num))
% sigMat_noise_truncated = filter_function(sigMat_noise_truncated, filter_f, datainfo_noise.HWDesc.SamplingFrequency);
% time_res_noise = 4;     % just to match the t_res of the numerical phantom input record length
% [t_noise, ~] = formInterpolationVec(datainfo_noise, n, time_res_noise, image_width);
% % ts_noise = 0 : 1/datainfo_noise.HWDesc.SamplingFrequency : (datainfo_noise.MeasurementDesc.RecordLength-1)/datainfo_noise.HWDesc.SamplingFrequency; % sampling instants
% tSamplesTotNum = detNum*length(t_noise);
% 
% eMean = zeros(tSamplesTotNum,1);
% switch E_MODE
%     % hardcoding estimators (assuming uncorrelated Gaussian (i.e. white) noise as in Tick 2016)
%     case 'white standard'
%         % eCov = eye(tSamplesTotNum);
%         eCov = speye(tSamplesTotNum);     % sparse
%         
%     case 'white with sd prop to pDet peaks'         % like Tick 2016
%      % STILL TO CALC
%         
%     case 'sample estimator'
%     
%         repVec_intrpltd = linspace(1,repNum, (tSamplesTotNum)+1);
%         rep_vs_tSamples_intrpltd = zeros(slcNum,wlNum,length(repVec_intrpltd),tSamplesTotNum);
%         
%         for run_idx = 1:1
%             for slc_idx = 1 : slcNum
%                 for wl_idx = 1 : wlNum
%         
%                     rep_vs_tSamples = zeros(repNum,tSamplesTotNum);
%                     for rep_idx = 1 : repNum
%                         % TO-DO: use noise params!! 
%                         sigMat_currentRep = sigMat_noise(:, :, run_idx, slc_idx, rep_idx, wl_idx);
%         
%                         % interpolate raw data to smaller resolution (downsampling) & truncate (same as in recon code)
%                         sigMat_interpTrunc = zeros(length(t), size(sigMat_currentRep, 2));
%                         for j = 1:detNum
%                             sigMat_interpTrunc(:,j) = interp1(ts, sigMat_currentRep(:,j), t);
%                         end
%         
%                         % reshape rep to a row-vec (i.e. [det1 det2 det3 ....]) &
%                         % fill in samplesVsRep mat
%                         sigMat_currRepRowVec = reshape(sigMat_interpTrunc, [1 tSamplesTotNum]);
%                         rep_vs_tSamples(rep_idx,:) = sigMat_currRepRowVec;
%                     end
%         
%                     % interpolate rep_vs_tSamples mat along col dim (i.e. reps) to
%                     % reach #obs > #RVs (detNum*samplesNum)
%                     for j = 1:tSamplesTotNum
%                         rep_vs_tSamples_intrpltd(slc_idx,wl_idx,:,j) = interp1(1:repNum, rep_vs_tSamples(:,j), repVec_intrpltd);
%                     end
%                     % Q: I'm concerned that the cov matrix resulting from this
%                     % won't be full rank as the rows will be lin depndent on
%                     % the original rows used in this interpolation (truly acquired)
%                 end
%             end
%         end
%         
%         
%         % clear sigMat rep_vs_tSamples sigMat_currentRep sigMat_currRepRowVec sigMat_interpTrunc repVec_intrpltd;
%         
%         % take the mean along the "reps" dim
%         eCov = zeros(slcNum, wlNum, tSamplesTotNum, tSamplesTotNum);
%         for slc_idx = 1 : slcNum
%             for wl_idx = 1 : wlNum
%                 eMean(slc_idx,wl_idx,:) = mean(rep_vs_tSamples_intrpltd(slc_idx, wl_idx, :, :), 3);
%             end
%         end
%         
%         % form the Cov matrix (per Zpos, per wl)
%         for slc_idx = 1 : slcNum
%             for wl_idx = 1 : wlNum
%                 tic;
%                 eCov(slc_idx,wl_idx,:,:) = cov(squeeze(rep_vs_tSamples_intrpltd(slc_idx, wl_idx, :, :)));
%                 toc;
%             end
%         end
%         
% end
% 
% 
% 
% 
% % % calc cholesky decomp & its inv of every cov mat (per slc per wl)
% % e_cholInv = zeros(size(eCov));
% % for slc_idx = 1 : slcNum
% %     for wl_idx = 1 : wlNum
% %         tic;
% %         [temp,p] = chol(squeeze(eCov(slc_idx,wl_idx,:,:)), 'lower');
% %         e_cholInv(slc_idx,wl_idx,:,:) = inv(temp);
% %         toc;
% %     end
% % end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INVERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% %%%% calc MAP estimate of p0 (analytical way - eq 10 in Tick 2016)
% For averaged input (to test MAP inversion in noise), set repEnd > 1
% This part is not effective when working with numerical data
path_datasetIn = 'E:\Marwan\Datasets_offline\MSOT 256\(2017-02-04) uSpheres for basis cov matrix (JP idea)\depth50mm_res200um\';     % don't forget the backslash in the end
slcInIdx = 6;
wlInIdx = 1;
repBegin = 1;            % bet 1 & 500, when changed --> only need to recalc b in eq.10
repEnd = 1;

n_repIn = 100;          % im size of input
ADD_NOISE = 1;          % if add synth noise to pDetIn
snr = 1;
BPF_INPUT_Pdet = 1;
filter_min = 0.02e6;
filter_max = 7e6;
P0_MODE = 'weighted';

filter_f = [filter_min filter_max];      % a band-pass filter (originlay between 100kHz and 8MHz), or maybe try the low pass a little lower (it should be depending on the BW of the transducer)
MB_MAT_FOLDER = 'E:\Marwan\MB_matrices';
A_mat = load([MB_MAT_FOLDER '\A_mat_t_res_4_100x100_width_25_c_1530_nDet_256_t_566']);                 % (Q): is this correct??
A_mat = A_mat.A_mat;

% extract input pDet rep as an input to the inversion algo
IP_SRC = 'numerical_MouseCart';                  % 'MouseCart' - 'Acuity' - 'numerical_MouseCart'
switch IP_SRC
    case 'MouseCart'
        if ~exist('sigMat_pathName_pDetIn', 'var') || ~strcmp(sigMat_pathName_pDetIn, path_datasetIn) 
            [datainfo_pDetIn, sigMat_pDetIn, sigMat_pathName_pDetIn] = loadSigMat_iThera(path_datasetIn);
        end
        disp(['input to inv algo is loaded from dataset "' sigMat_pathName_pDetIn '"']);
        
        [t, ts] = formInterpolationVec(datainfo_pDetIn, n_repIn);
        
        pDetInStack = 0;
        for repIdx = repBegin:repEnd
            pDetIn = squeeze(sigMat_pDetIn(:,:,wlInIdx,slcInIdx,repIdx));
            
            if ADD_NOISE
                % add synth noise
                pDetIn = awgn(pDetIn, snr, 'measured');
            end
            
            if BPF_INPUT_Pdet
                pDetIn = filter_function(pDetIn, filter_f, datainfo_pDetIn.HWDesc.SamplingFrequency);              % (Q): should I BPF input pDet (also, I should filter only selected frame)
            end
            
            % interpolate raw data to smaller resolution (downsampling) & truncate
            temp = zeros(length(t), size(pDetIn, 2));
            for j = 1:datainfo_pDetIn.HWDesc.NumDetectors
                temp(:,j) = interp1(ts, pDetIn(:,j), t);
            end
            pDetIn = temp;
            
            % reshape interpolated measures to column-major format
            pDetIn = reshape(pDetIn, length(t)*datainfo_pDetIn.HWDesc.NumDetectors, 1);
            
            % accumulate to calc average
            pDetInStack = pDetInStack + pDetIn;
        end
        % calc avg
        pDetIn = pDetInStack/(repEnd-repBegin+1);     
    case 'Acuity'
        % should load pDetIn from hh dataset (still to complete!! --> BPF done but still adding wt noise & truncating to t vec like MouseCart above)
        myReadAcuityData;
    case 'numerical_MouseCart'
        pDetIn_struct = load('E:\Marwan\Dropbox\PA_imaging\MSOT_Recon_Detection\numerical phantoms\pDetIn_derenzo_timeRes_4.mat');
        pDetIn = pDetIn_struct.sigMat_numerical;
        if ADD_NOISE
            % add synth noise
            pDetIn = awgn(pDetIn, snr, 'measured');
        end
        
end

%% do Bayesian inversion (MAP)
% calc A & b expression for the MAP estimate (eq 10)
p0meanRaw = p0mean;
p0covRaw = p0cov;       % save it before trying weighting to avoid recalcing it
if(strcmp(P0_MODE,'weighted'))
    %%% do heuristic recon (BPn) for input to use this heuristic as a weighting for both mean vec & cov mat 
    disp('prior mode: WEIGHTED (JP IDEA)');
    nDet = length(angle_sensor);
    pDetIn = reshape(pDetIn, [length(pDetIn_struct.t) nDet]);     % resize to work with BP code
    [pDetInWeighting, X, Y] = backproject_luis(pDetIn, pDetIn_struct.n , datainfo.HWDesc.Radius,...
                                            angle_sensor, c, 'full', pDetIn_struct.t, datainfo.HWDesc.SamplingFrequency, pDetIn_struct.image_width, 0, 0);
    pDetIn = reshape(pDetIn, [length(pDetIn_struct.t)*nDet 1]);     % resize back                                    
    if NONNEG
        pDetInWeighting = max(pDetInWeighting,0);
    end
    % if BPn already done on the input before, just load it
%     pDetInWeighting = load([path_priorDataset '\recons\recons to build prior dist\BPn\no added noise\reconBackProjection_imSz100_sos1510_rep1to1000.mat']);
%     pDetInWeighting = pDetInWeighting.reconStack;
%     pDetInWeighting = squeeze(pDetInWeighting(:,:,wlInIdx,slcInIdx,repBegin:repEnd));
    pDetInWeighting = reshape(pDetInWeighting, [size(pDetInWeighting,1)^2 1]);
    
    %%% weight p0 basis estimators
    p0mean = pDetInWeighting.*p0meanAvg';
%     p0cov = diag(pDetInWeighting)*p0covAvg;
    p0cov = diag(pDetInWeighting)*p0covAvg*diag(pDetInWeighting);     %% now it's symmetric
end
p0covInv = inv(p0cov);

%%% calc A (WARNING: % not sparse (takes LONG TIME 608 s))
% %u don't need to calc it when u change the input pDetIn (it's dependent on noise covariance, conventional A_mat & prior covariance)
% tic;
% A = A_mat'*(eCov\A_mat) + p0covInv;
% disp(['elapsed time for A calc is ' num2str(toc) ' s']);

%%% calc B
% WARNING: make sure that A_mat, eCov, eMean & pDetIn have the same sizes concerning the recordLength (length of a single channel)
tic;
b = A_mat'*(eCov\(pDetIn-eMean)) + (p0cov\p0mean);            % takes 9 s
disp(['elapsed time for b calc is ' num2str(toc) ' s']);

% do MAP recon
p0_MAP = A\b;                % takes only 4 secs
p0_MAPvecMajor = p0_MAP;
p0_MAP = reshape(p0_MAPvecMajor, [n_repIn n_repIn]);

figure, imagesc(p0_MAP), title('MAP recon'), colormap(bone);
% saveName = ['p0_MAP_addNoiseSnr' num2str(snr) '_avg' num2str(repEnd-repBegin+1)];
% saveName = ['p0_MAP_phntm_slc1_noNoise'];
saveName = ['p0_MAP_weightedPriorBasis_weightingWithDADNumPhatnom_slc' num2str(slcInIdx) '_snr_' num2str(snr) ];
savefig(['E:\Marwan\Dropbox\PA_imaging\MAP recons\' saveName '.fig']);
save([path_priorDataset '\recons\MAP recons\' saveName '.mat'], 'p0_MAP');

%%% Calc CNR to compare
% targetROI.ul.x = 36;
% targetROI.ul.y = 84;
% targetROI.lr.x = 41;
% targetROI.lr.y = 87;
% bgdROI.ul.x = ;
% bgdROI.ul.y = ;
% bgdROI.lr.x = ;
% bgdROI.lr.y = ;



%%% Compare to normal MB
NONNEG = 0;                     % 0 --> normal lsqr
regu = 20;
MB_regu = 1e6*(2^((regu-20)));
p0convMB = reconstruction(A_mat, pDetIn, n_repIn, 'MB_Tik', MB_regu, NONNEG);
figure, imagesc(p0convMB), title('conventional MB Tik recon'), colormap(bone);
saveName = ['p0_convMBtik_weightedPriorBasis_weightingWithDADNumPhatnom_slc' num2str(slcInIdx) '_snr_' num2str(snr) ];
savefig(['E:\Marwan\Dropbox\PA_imaging\MAP recons\' saveName '.fig']);
save([path_priorDataset '\recons\MAP recons\' saveName '.mat'], 'p0convMB');



