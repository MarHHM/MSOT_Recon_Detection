% Recon code
%% VIPs:
%     - choose correct path & desired recon ranges (for zpos, wls & reps)
%     - don't forget to create the folder "recons" & "masks" in the dataset path for saving recons
%     - don't forget to change the name of the saved recon to reflect the last changes u did (e.g. experimented with)
%     - make sure that the dataset u will load is not huge if you use a normal PC (limited RAM)
%     - if mouse data, 800nm is the important wl to show (isosbestic, best total blood signal)

%% PATHS & PARAMS
% addReconPaths;        % no need, it's in the startup.m
newDataPath =  'S:\PA_DATASETS\MSOT 256\Hong\phantom\Scan_69';      % \Hong\phantom\Scan_69  ||  \Marwan\(2017-03-03) testing phantom ink+straw\both horizontal
im_w = 25e-3;                % physical im width at acquisition (m) - org: 25e-3
reconRes = 75e-6;                  % resolution of reconstruction (m/pixel) (org: 100 - preferred at 100 or 75 um)

%%% PARAMS
RECON_METHOD = 'MB_Tik';        % 'BP'  'MB_Tik'  'MB_TVL1'
RECON_ALL = false;
SAVE_RECON = true;


%% %%%%%%%%%%%%%%%%%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MB_MAT_FOLDER = 'S:\MB_matrices';
IMPOSE_NONNEGATIVITY = false;                     % 0 --> allow negatives
f_min = 0.04e6;             % typically 0.04 (for unmixing, better 0.25) to 8 MHz
f_max = 8e6;                
t_res = 2;                  % time resolution for model-based (bigger -> finer res but bigger A & more recon time) (originally 2)
MB_regu = 1e6;              % regulariztaion param for numerical solution (typically 1e6)

n = floor(im_w/reconRes);   % reconstructed im size (pixels)

%%% load sigMat to read dims & params
if ~exist('sigMat_pathName', 'var') || ~strcmp(sigMat_pathName,newDataPath) || ~exist('sigMat', 'var')
    [datainfo, sigMat, sigMat_pathName] = loadSigMat_iThera( newDataPath );
else
    disp(['--scan "' sigMat_pathName '" is already loaded--']);
end

%%% measurements to recon (limit #reps due to memory limitation of prefilt fun (no need if on WS))
if exist('sigMat', 'var')
    run_idx = 1;
    if RECON_ALL
        disp('the full dataset will be reconstructed..');
        zpos_begin = 1;       zpos_end = size(sigMat, 4);
        rep_begin = 1;        rep_end = size(sigMat, 5);
        wl_begin = 1;         wl_end = size(sigMat, 6);
    else
        zpos_begin = 1;                            % 5
        zpos_end = 1;                              % 5
        rep_begin = 1;
        rep_end = 1;
        wl_begin = 6;                               % 6 is 750nm for Hong phantom (scan_69 -> AlexaFl_750)
        wl_end = 6;
    end
    
    %%% extract a stack from the whole measures to recon
    sigMat_truncated = sigMat(:, :, run_idx, zpos_begin:zpos_end, rep_begin:rep_end, wl_begin:wl_end);
    
    %%% correct laser energy (NO NEED; done auto before loading the sigs by iThera function)
    % sigMat_truncated = normByLaserEnergy(sigMat_truncated, datainfo);
    
    %%% pre-filtration (is not optimized for large datasets (e.g. high rep num))
    filter_f = [f_min f_max];      % a band-pass filter (originally between 100kHz and 8MHz), or maybe try the low pass a little lower (it should be depending on the BW of the transducer)
    sigMat_truncated = filter_function( sigMat_truncated, filter_f, datainfo.HWDesc.SamplingFrequency );
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
[t, ts] = formInterpolationVec(datainfo, n, t_res, im_w, c);

if strcmp(RECON_METHOD,'MB_Tik')||strcmp(RECON_METHOD,'MB_TVL1')
    % set MB recon params
    n_angles = 2*n;                                     % number of points for discretizing the curve
    sizeT = length(t);
    
    % check for model matrix if saved, otherwise build & save it
    A_matPath = [MB_MAT_FOLDER '\A_mat_t_res_' num2str(t_res) '_' num2str(n) 'x' num2str(n) '_width_' num2str(im_w*1e3)...
        '_c_' num2str(c) '_nDet_' num2str(length(angle_sensor)) '_t_' num2str(length(t)) '.mat'];
    A_mat = compOrLoadA_mat( A_matPath, c, n, im_w, t, datainfo.HWDesc.Radius, angle_sensor, n_angles );
end

%% Recon the extracted stack
Recon_MB = zeros(n, n, run_idx, zpos_end-zpos_begin+1, rep_end-rep_begin+1, wl_end-wl_begin+1);
reconItr = 0;
totNumRecons = size(Recon_MB,3)*size(Recon_MB,4)*size(Recon_MB,5)*size(Recon_MB,6);


disp(['Reconstruction started (Method: ' RECON_METHOD ')..'])
for run_idx = 1:1
    for zpos_idx = 1 : zpos_end-zpos_begin+1
        for wl_idx = 1 : wl_end-wl_begin+1
            for rep_idx = 1 : rep_end-rep_begin+1
                tic;
                % extract current measure (2d) to recon it
                sigMat_current = sigMat_truncated(:, :, run_idx, zpos_idx, rep_idx, wl_idx);
                
                if strcmp(RECON_METHOD,'BP')
                    [Recon_tmp, X, Y] = backproject_luis(sigMat_current, n, datainfo.HWDesc.Radius, angle_sensor, c, 'full', ts, im_w, 0, 0);
                    if IMPOSE_NONNEGATIVITY
                        Recon_tmp = max(Recon_tmp,0);
                    end
                end
                
                if strcmp(RECON_METHOD,'MB_Tik')||strcmp(RECON_METHOD,'MB_TVL1')
                    
                    b_vec = prepareMeasuredPressureVec( sigMat_current, ts, t, datainfo.HWDesc.NumDetectors );
                    
                    % do recon
                    Recon_tmp  = reconstruction(A_mat, b_vec, n, RECON_METHOD, MB_regu, IMPOSE_NONNEGATIVITY);
                end
                
                Recon_MB(:, :, run_idx, zpos_idx, rep_idx, wl_idx) = Recon_tmp;
                reconItr = reconItr+1;
                disp(['recon ' num2str(reconItr) ' of ' num2str(totNumRecons) ' done..']);
                toc;
            end
        end
    end
end

if SAVE_RECON
    % don't forget to create the folder "recons" in this path
    % note = 'sigMat_truncated attached with this recon is laser-energy corrected & filtered';
    save(strcat(sigMat_pathName, 'recons\', RECON_METHOD, ' - nonNeg_', num2str(IMPOSE_NONNEGATIVITY), '-zpos_', num2str(zpos_end-zpos_begin+1),...
        '-reps_', num2str(rep_end-rep_begin+1), '-wls_', num2str(wl_end-wl_begin+1), '-reconRes_', num2str(reconRes),...
        '-imW_', num2str(im_w) ,'.mat'), 'Recon_MB',...
        'datainfo', 'im_w', 'reconRes', 'f_min', 'f_max', 't_res', 'ts', 't', 'angle_sensor',...
        'MB_regu', 'zpos_begin', 'zpos_end', 'wl_begin', 'wl_end', 'rep_begin', 'rep_end', 'IMPOSE_NONNEGATIVITY', '-v7.3');
    disp('reconstruction saved to disk.');
    
end
