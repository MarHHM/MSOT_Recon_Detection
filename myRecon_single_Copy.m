% Recon code
% VIPs:
%     - don't forget to create the folder "recons" in the dataset path for saving recons
%     - don't forget to change the name of the saved recon to reflect the last changes u did (e.g. experimented with)
%     - make sure that the dataset u will load is not huge if you use a normal PC (limited RAM)
%     - if mouse data, 800nm is the important wl to show (isosbestic, best total blood signal)

addReconPaths;
% addpath(genpath('C:\Users\marwan.muhammad\Dropbox\PA_imaging\Luis code\'));
newDataPath = 'D:\DATASETS\MSOT 256\Evangelos (high reflections)\Scan_90\';
% 'D:\DATASETS\MSOT 256\Evangelos (high reflections)\Scan_87\'
% 'D:\DATASETS\MSOT 256\Scan_96 (Ivan - cleanest data with tumor)\'
% 'D:\DATASETS\MSOT 256\Scan_70 (Ivan)\'
% 'E:\Marwan\Dropbox\PA_imaging\MSOT_Recon_Detection\numerical phantoms\'
% 'D:\DATASETS\MSOT 256\Marwan\(2017-03-03) testing phantom ink+straw\both horizontal\'
% 'D:\DATASETS\MSOT 256\Marwan\(2017-03-03) testing phantom ink+straw\straw up\'
% 'D:\DATASETS\MSOT 256\Marwan\(2017-02-22) testing phantom ink+straw (FALSE SPECTRUM!!)\averaging 10 frames\'
% 'D:\DATASETS\MSOT 512\(2017-02-20) Mouse data by Steven\'


%% params
% if sigMat was generated using numerical phantom propagated through model matrix, use the same dims (n - image_width) that were used there
n = 250;        % reonstructed im size (pixels) (org: better to recon at 100um resolution (e.g. if width = 25mm, recon at n = 250 pixles))
image_width = 25e-3;        % physical im width at acquisition (m)

filter_min = 0.1e6;             % typically bet 0.1 to 7 MHz
filter_max = 7e6;
time_res = 2;                                       % time resolution for model-based (bigger -> finer res but bigger A & more recon time) (originally 2)
RECON_METHOD = 'MB_Tik';        % 'BP'  'MB_Tik'  'MB_TVL1'
NONNEG = 0;                     % 0 --> normal lsqr
SAVE_RECON = 1;

MB_MAT_FOLDER = 'D:\MB_matrices';           % 'D:\MB_matrices'


%% %%%%%%%%%%%%%%%%%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% load sigMat to read dims & params
if ~exist('sigMat_pathName', 'var') || ~strcmp(sigMat_pathName,newDataPath) || ~exist('sigMat', 'var')
    [datainfo, sigMat, sigMat_pathName] = loadSigMat_iThera(newDataPath);
else
    disp(['scan "' sigMat_pathName '" is already loaded']);
end

%%% measurements to recon (limit #reps due to memory limitation of prefilt fun (no need if on WS))
if exist('sigMat', 'var')
    slc_begin = 26;
%     slc_end = size(sigMat, 4);
    slc_end = 29;
    wl_begin = 10;
%     wl_end = size(sigMat, 6);
    wl_end = 12;
    rep_begin = 1;
    % rep_end = 1;
    rep_end = size(sigMat, 5);
    run_idx = 1; 
    
    %%% extract a stack from the whole measures to recon
    sigMat_truncated = sigMat(:, :, run_idx, slc_begin:slc_end, rep_begin:rep_end, wl_begin:wl_end);
    clear sigMat;
    
    % if ADD_NOISE_TO_RAW_DATA
    %     for run_idx = 1:1
    %         for slc_idx = 1 : slc_end-slc_begin+1
    %             for wl_idx = 1 : wl_end-wl_begin+1
    %                 for rep_idx = 1 : rep_end-rep_begin+1
    %                     for detIdx = 1:size(sigMat_truncated,2)
    %                         sigMat_truncated(:,detIdx,run_idx,slc_idx,rep_idx,wl_idx) = awgn(sigMat_truncated(:,detIdx,run_idx,slc_idx,rep_idx,wl_idx),...
    %                                                                                             snr, 'measured');
    %                     end
    %                 end
    %             end
    %         end
    %     end
    % end
    
    % pre-filtration (is not optimized for large datasets (e.g. high rep num))
    filter_f = [filter_min filter_max];      % a band-pass filter (originally between 100kHz and 8MHz), or maybe try the low pass a little lower (it should be depending on the BW of the transducer)
    sigMat_truncated = filter_function(sigMat_truncated, filter_f, datainfo.HWDesc.SamplingFrequency);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 1:length(sos_vec)
%     msgbox_h = msgbox(['currently at sos ' num2str(i) ' of ' num2str(length(sos_vec))] , 'current SOS');

% reading general params from datainfo
try
    angle_sensor = datainfo.HWDesc.StartAngle : datainfo.HWDesc.StepAngle : datainfo.HWDesc.EndAngle;
catch     % error due to MSOT 512, try to fill them using same values from MSOT 256
    datainfo.HWDesc.StartAngle = -0.7762;
    datainfo.HWDesc.EndAngle = 3.9178;
    datainfo.HWDesc.NumDetectors = 512;
    datainfo.HWDesc.StepAngle = (datainfo.HWDesc.EndAngle - datainfo.HWDesc.StartAngle)/(datainfo.HWDesc.NumDetectors-1);
    angle_sensor = datainfo.HWDesc.StartAngle : datainfo.HWDesc.StepAngle : datainfo.HWDesc.EndAngle;
    
    datainfo.HWDesc.Radius = 0.0405;
end

T = datainfo.AverageTemperature;
c = 12 + round(1.402385 * 1e3 + 5.038813 * T - 5.799136 * 1e-2 * T^2 + 3.287156 * 1e-4 * T^3 - 1.398845 * 1e-6 * T^4 + 2.787860 * 1e-9 * T^5 );
%     c = sos_vec(i);
[t, ts] = formInterpolationVec(datainfo, n, time_res, image_width, c);

if strcmp(RECON_METHOD,'MB_Tik')||strcmp(RECON_METHOD,'MB_TVL1')
    % set MB recon params
    n_angles = 2*n;                                     % number of points for discretizing the curve
    sizeT = length(t);
    regu = 20;
    MB_regu = 1e6*(2^((regu-20)));
    
    % check for model matrix if saved, otherwise build & save it
    A_matPath = [MB_MAT_FOLDER '\A_mat_t_res_' num2str(time_res) '_' num2str(n) 'x' num2str(n) '_width_' num2str(image_width*1e3)...
                '_c_' num2str(c) '_nDet_' num2str(length(angle_sensor)) '_t_' num2str(length(t)) '.mat'];
    A_mat = compOrLoadA_mat(A_matPath, c, n, image_width, t, datainfo.HWDesc.Radius, angle_sensor, n_angles);
end

%% recon the extracted stack
Recon_MB = zeros(n, n, run_idx, slc_end-slc_begin+1, rep_end-rep_begin+1, wl_end-wl_begin+1);
reconItr = 0;
totNumRecons = size(Recon_MB,3)*size(Recon_MB,4)*size(Recon_MB,5)*size(Recon_MB,6);

for run_idx = 1:1
    for slc_idx = 1 : slc_end-slc_begin+1
        for wl_idx = 1 : wl_end-wl_begin+1
            for rep_idx = 1 : rep_end-rep_begin+1
                tic;
                % extract current measure (2d) to recon it
                sigMat_current = sigMat_truncated(:, :, run_idx, slc_idx, rep_idx, wl_idx);
                
                if strcmp(RECON_METHOD,'BP')
                    [Recon_tmp, X, Y] = backproject_luis(sigMat_current, n, datainfo.HWDesc.Radius, angle_sensor, c, 'full', ts, image_width, 0, 0);
                    if NONNEG
                        Recon_tmp = max(Recon_tmp,0);
                    end
                end
                
                if strcmp(RECON_METHOD,'MB_Tik')||strcmp(RECON_METHOD,'MB_TVL1')
                    % interpolate raw data to smaller resolution (downsampling) & truncate
                    sigMat2 = zeros(length(t), size(sigMat_current, 2));
                    for j = 1:datainfo.HWDesc.NumDetectors
                        sigMat2(:,j) = interp1(ts, sigMat_current(:,j), t);
                    end
                    clear sigMat_current;
                    
                    % reshape interpolated measures to column-major format
                    b_vec = reshape(sigMat2, sizeT*datainfo.HWDesc.NumDetectors, 1);
                    clear sigMat2
                    
                    % do recon
                    Recon_tmp  = reconstruction(A_mat, b_vec, n, RECON_METHOD, MB_regu, NONNEG);
                end
                
                Recon_MB(:, :, run_idx, slc_idx, rep_idx, wl_idx) = Recon_tmp;
                clear Recon_tmp;
                reconItr = reconItr+1;
                disp(['recon ' num2str(reconItr) ' of ' num2str(totNumRecons) ' done..']);
                toc;
            end
        end
    end
end

if SAVE_RECON
    % don't forget to create the folder "recons" in this path
    save(strcat(sigMat_pathName, 'recons\recon', RECON_METHOD, '_imSz_', num2str(n), '_slcs_', num2str(slc_end-slc_begin+1),...
                                                               '_wls_', num2str(wl_end-wl_begin+1), '_reps_', num2str(rep_end-rep_begin+1), '.mat'), 'Recon_MB',...
                                                               'datainfo', 'image_width','slc_begin', 'slc_end', 'wl_begin', 'wl_end', 'rep_begin', 'rep_end');

    % save(strcat(sigMat_pathName, 'recons\recon', RECON_METHOD, '_imSz', num2str(n), '_sos', num2str(c), '_rep', num2str(rep_begin), 'to', num2str(rep_end), '.mat'), 'Recon_MB');
end
% delete(msgbox_h);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

