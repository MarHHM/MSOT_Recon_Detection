function recon_single_msot_scan(scan_path, NON_NEG, recon_lims, RECON_ALL, n_det__custom, addtnl_note)
%% VIPs:
%     - choose correct path & desired recon ranges (for zpos, wls & reps)
%     - don't forget to change the name of the saved recon to reflect the last changes u did (e.g. experimented with)
%     - make sure that the dataset u will load is not huge if you use a normal PC (limited RAM)
%     - if mouse data, 800nm is the important wl to show (isosbestic, best total blood signal)
random_temp__path = tempname;   % creates a temporary file with different name at each run
diary(random_temp__path);       % activate logging of command window
%% PATHS & PARAMS
MB_MAT__path = 'C:\PA_local\MB_matrices';
RECON_METHOD = "MB_Tik";        % 'BP'  'MB_Tik'  'MB_TVL1'
im_w = 25*1e-3;         % physical im width at acquisition (m) - org: 25e-3
res = 75*1e-6;          % resolution of reconstruction (m/pixel) (org: 100 - preferred at 100 or 75 um)
f_min = 0.04e6;         % typically 0.04 (for unmixing, better 0.25)
f_max = 8e6;            % typically 8M
MB_regu = 1e6;          % regulariztaion param for numerical solution (typically 1e6)
t_res = 2;              % time resolution for model-based (bigger -> finer res but bigger A & more recon time) (originally 2)
% f_min = 0;
% f_max = 5e6;          % typically 8M
% MB_regu = 0;          % regulariztaion param for numerical solution (typically 1e6)
SAVE_RECON = true;


%% %%%%%%%%%%%%%%%%%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['---------------- Reconstructing scan : ' scan_path ' - nonneg = ' num2str(NON_NEG) '----------------']);

n = floor(im_w/res);   % reconstructed im size (pixels)
%%% load sigMat to read dims & params
if ~exist('sigMat_pathName', 'var') || ~strcmp(sigMat_pathName,scan_path) || ~exist('sigMat', 'var')
    [datainfo, sigMat, ~] = loadSigMat_iThera( scan_path );
else
    disp(['--scan "' sigMat_pathName '" is already loaded--']);
end

% disp("replacing numDetectors in datainfo.HWDesc.NumDetectors with a custom value..");
% datainfo.HWDesc.NumDetectors = custom_numDetectors;

%%% measurements to recon (limit #reps due to memory limitation of prefilt fun (no need if on WS))
run_idx = 1;
rep_begin = 1;        rep_end = datainfo.RepNum;  %usually only 1 rep
if RECON_ALL
    disp('!! THE FULL SCAN WILL BE RECONSTRUCTED !!');
    zpos_begin = 1;       zpos_end = length(datainfo.ZPositions);
    wl_begin = 1;         wl_end = length(datainfo.Wavelengths);
else
    zpos_begin = recon_lims.zpos(1);      zpos_end = recon_lims.zpos(2);
    wl_begin = recon_lims.wls(1);         wl_end = recon_lims.wls(2);
end

%%% extract a stack from the whole measures to recon
sigMat_truncated = sigMat(:,:,run_idx,zpos_begin:zpos_end,rep_begin:rep_end,wl_begin:wl_end);



disp("[TEST] n_det__custom = "+int2str(n_det__custom)+" (dets at both ends are zeroed)");
n_det__zeroed = datainfo.HWDesc.NumDetectors - n_det__custom;
sigMat_truncated(:,(1:floor(n_det__zeroed/2)),:,:,:,:) = 0;
sigMat_truncated(:,(end-(floor(n_det__zeroed/2))+1:end),:,:,:,:) = 0;


%%% correct laser energy (NO NEED; done auto before loading the sigs by iThera function)
% sigMat_truncated = normByLaserEnergy(sigMat_truncated, datainfo);

%%% pre-filtration (is not optimized for large datasets (e.g. high rep num))
filter_f = [f_min f_max];      % a band-pass filter (originally between 100kHz and 8MHz), or maybe try the low pass a little lower (it should be depending on the BW of the transducer)
sigMat_truncated = filter_function( sigMat_truncated, filter_f, datainfo.HWDesc.SamplingFrequency );

%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 1:length(sos_vec)
%     msgbox_h = msgbox(['currently at sos ' num2str(i) ' of ' num2str(length(sos_vec))] , 'current SOS');

% reading general params from datainfo
try
    angle_sensor = datainfo.HWDesc.StartAngle : datainfo.HWDesc.StepAngle : datainfo.HWDesc.EndAngle;
catch     % error due to different MSOT format - check if new reader is required from iThera
    disp("!! different MSOT format - check if new reader is required from iThera..");
    datainfo.HWDesc.Radius = 0.0405;
    datainfo.HWDesc.StartAngle = -0.7762;
    datainfo.HWDesc.EndAngle = 3.9178;
    datainfo.HWDesc.NumDetectors = 256;
    datainfo.HWDesc.StepAngle = (datainfo.HWDesc.EndAngle - datainfo.HWDesc.StartAngle)/(datainfo.HWDesc.NumDetectors-1);
    angle_sensor = datainfo.HWDesc.StartAngle : datainfo.HWDesc.StepAngle : datainfo.HWDesc.EndAngle;
end

%%%
% disp("using sos for heavy water instead of calulating it from the metadat of 'datainfo'!!");
% c = 1465;
%%%
T = datainfo.AverageTemperature;
c = 12 + round(1.402385 * 1e3 + 5.038813 * T - 5.799136 * 1e-2 * T^2 + 3.287156 * 1e-4 * T^3 - 1.398845 * 1e-6 * T^4 + 2.787860 * 1e-9 * T^5 );
[t, ts] = formInterpolationVec(datainfo, n, t_res, im_w, c);

if RECON_METHOD == "MB_Tik" || RECON_METHOD == "MB_TVL1"
    % set MB recon params
    n_angles = 2*n;                                     % number of points for discretizing the curve
    %     sizeT = length(t);
    
    % check for model matrix if saved, otherwise build & save it
    A_matPath = [MB_MAT__path '\A_mat_t_res_' num2str(t_res) '_' num2str(n) 'x' num2str(n) '_width_' num2str(im_w*1e3)...
        '_c_' num2str(c) '_nDet_' num2str(length(angle_sensor)) '_t_' num2str(length(t)) '.mat'];
    A_mat = compOrLoadA_mat( A_matPath, c, n, im_w, t, datainfo.HWDesc.Radius, angle_sensor, n_angles );
else
    A_matPath = "";
end

%% Recon the extracted stack
Recon = zeros(n, n, run_idx, zpos_end-zpos_begin+1, rep_end-rep_begin+1, wl_end-wl_begin+1);
reconItr = 0;
totNumRecons = size(Recon,3)*size(Recon,4)*size(Recon,5)*size(Recon,6);


disp(['Reconstruction started (Method: ' RECON_METHOD ')..'])
for run_idx = 1:1
    for zpos_idx = 1 : zpos_end-zpos_begin+1
        for wl_idx = 1 : wl_end-wl_begin+1
            for rep_idx = 1 : rep_end-rep_begin+1
                tic;
                % extract current slice (2d) to recon it
                sigMat_current = sigMat_truncated(:, :, run_idx, zpos_idx, rep_idx, wl_idx);
                
                if RECON_METHOD == "BP"
                    [Recon_tmp, ~, ~] = backproject_luis(sigMat_current, n, datainfo.HWDesc.Radius, angle_sensor, c, 'full', ts, datainfo.HWDesc.SamplingFrequency,...
                        im_w, 0, 0);
                    if NON_NEG
                        Recon_tmp = max(Recon_tmp,0);
                    end
                elseif RECON_METHOD == "MB_Tik" || RECON_METHOD == "MB_TVL1"
                    b_vec = prepareMeasuredPressureVec( sigMat_current, ts, t, datainfo.HWDesc.NumDetectors );
                    Recon_tmp  = reconstruction(A_mat, b_vec, n, RECON_METHOD, MB_regu, NON_NEG);
                end
                
                Recon(:, :, run_idx, zpos_idx, rep_idx, wl_idx) = Recon_tmp;
                reconItr = reconItr+1;
                
                elpsd_t = toc;
                disp(['recon ' num2str(reconItr) ' of ' num2str(totNumRecons) ' done (' num2str(elpsd_t) 's)..']);
            end
        end
    end
end

if SAVE_RECON
    if ~exist([scan_path '\recons'], 'dir')
        mkdir([scan_path '\recons']);
    end
    reconStruct__fName = strcat(RECON_METHOD, '-nonNeg_', num2str(NON_NEG), '-zposs_',num2str(zpos_end-zpos_begin+1), '-wls_', num2str(wl_end-wl_begin+1),...
        addtnl_note);
    scan_path = convertCharsToStrings(scan_path);
    if exist((scan_path+"\recons\"+reconStruct__fName), 'dir')                % to avoid overwriting an existing recon
        reconStruct__fName = (reconStruct__fName+"__");
    end
    reconStruct__path = (scan_path+"\recons\"+reconStruct__fName);
    mkdir(reconStruct__path);
    save((reconStruct__path+"\"+reconStruct__fName+".mat"),...
        'scan_path', 'RECON_METHOD', 'Recon', 'sigMat_truncated', 'datainfo', 'im_w', 'res', 'n', 'f_min', 'f_max','t_res', 'ts', 't',...
        'angle_sensor','MB_regu', 'run_idx', 'zpos_begin', 'zpos_end', 'wl_begin', 'wl_end', 'rep_begin', 'rep_end', 'NON_NEG',...
        'A_matPath', '-v7.3');       % sigMat_truncated attached with this recon is laser-energy corrected & filtered
    disp(("-> Reconstruction saved to "+reconStruct__path));
    
    % save a representative image (for ease of browsing later)
    rprsnttv_im = Recon(:,:,1,1,1,1);
    figure('Position', [1.3082e+03 549.8000 469.6000 400]), imagesc(rprsnttv_im),...
        title(("recon: "+reconStruct__fName+" (rprsnttv_im)"), 'Interpreter', 'none'), colormap(bone), colorbar, axis image off;
    export_fig((reconStruct__path+"\"+reconStruct__fName+".png"), '-nocrop');
end

diary off;
movefile(random_temp__path, (reconStruct__path+"\log.txt"));    % rename & relocate log file to recon folder
end