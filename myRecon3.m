javaaddpath MSOTBeans\msotbeans.jar
javaaddpath java_class\recon.jar
javaaddpath java_class\patbeans.jar
javaaddpath java_class\xmlbeans-2.5.0\lib\xbean.jar
addpath(genpath('D:\Marwan\Dropbox\PA_imaging\MSOT_Recon_Detection'));

% clearvars;

% params
n = 75;        % im width (pixels) (org: 200)
image_width = 20e-3;        % im width (m)
filter_min = 0.02e6;
filter_max = 7e6;
RECON_METHOD = 'MB_Tik';        % 'BackProjection'  'MB_Tik'  'MB_TVL1'
NONNEG = 0;
time_res = 2;                                       % time resolution for model-based (bigger -> finer res but bigger A & more recon time) (originally 2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loading into RAM (if not the same last data)
[FileNamesigMat, temp, FilterIndexsigMat] = uigetfile({'*.msot'}, 'mytitle', 'D:\Marwan\Datasets_offline\Marwan_1\Scan_5');
if ~strcmp(temp, PathNamesigMat)
    PathNamesigMat = temp;
    [datainfo] = loadMSOT([PathNamesigMat '\' FileNamesigMat]);
    [sigMat] = loadMSOTSignals(datainfo);
end

%% extract a stack from the whole measures to recon
% limit #reps due to memory limitiation of prefilt fun (no need if on WS)
Zpos_begin = 1;
Zpos_end = size(sigMat, 4);
wl_begin = 1;
wl_end = size(sigMat, 6);
rep_begin = 1;
rep_end = size(sigMat, 5);
% rep_end = 100;
run_idx = 1;

sigMat_truncated = sigMat(:, :, run_idx, Zpos_begin:Zpos_end, rep_begin:rep_end, wl_begin:wl_end);

% reading general params from datainfo
angle_sensor = datainfo.HWDesc.StartAngle : datainfo.HWDesc.StepAngle : datainfo.HWDesc.EndAngle;
ts = 0 : 1/datainfo.HWDesc.SamplingFrequency : (datainfo.MeasurementDesc.RecordLength-1)/datainfo.HWDesc.SamplingFrequency; % sampling instants
T = datainfo.AverageTemperature;
c = 12 + round(1.402385 * 1e3 + 5.038813 * T - 5.799136 * 1e-2 * T^2 + 3.287156 * 1e-4 * T^3 - 1.398845 * 1e-6 * T^4 + 2.787860 * 1e-9 * T^5 );
avgs = datainfo.MeasurementDesc.Averages;
runs = datainfo.RunNum;         % # full studies (usually 1)
repetitions = datainfo.RepNum;
slices = numel(datainfo.ZPositions);
lambdas = numel(datainfo.Wavelengths);

if strcmp(RECON_METHOD,'MB_Tik')||strcmp(RECON_METHOD,'MB_TVL1')
    % set MB recon params
    n_angles = 2*n;                                     % number of points for discretizing the curve
    limits(1) = datainfo.HWDesc.Radius-(image_width)*sqrt(2)/2;       % limits for the signal
    limits(2) = datainfo.HWDesc.Radius+(image_width)*sqrt(2)/2;       % limits for the signal
    dx = image_width/n;                                 % increment in x
    dt = dx/(time_res*c);                              % increment in t employed to make the model-based reconstruction
    fac = datainfo.HWDesc.SamplingFrequency/c;
    pos_start = max(1,int32((limits(1))*fac));
    pos_end = min(datainfo.MeasurementDesc.RecordLength,int32((limits(2))*fac));
    t = ts(pos_start):dt:ts(pos_end);           % downsampled (& cut) time vector  (less than ts)
    sizeT = length(t);
    regu = 20;
    MB_regu = 1e6*(2^((regu-20)));
    
    % check for model matrix if saved, otherwise build & save it
    init_folder_name = pwd;
    sub_folder_name = 'MB matrices';
    folder_name = [init_folder_name '\' sub_folder_name];
    if exist([folder_name '\' 'A_mat_t_res_',num2str(time_res),'_',num2str(n),'x',num2str(n),'_width_',num2str(image_width*1e3),'_c_',num2str(c),'.mat'], 'file')
        load ([folder_name '\' 'A_mat_t_res_',num2str(time_res),'_',num2str(n),'x',num2str(n),'_width_',num2str(image_width*1e3),'_c_',num2str(c),'.mat'])
    else
        A_mat = Calculate_MatrixMB_Luis(c,n,image_width,t,datainfo.HWDesc.Radius,angle_sensor,n_angles);
        save(([folder_name '\' 'A_mat_t_res_',num2str(time_res),'_',num2str(n),'x',num2str(n),'_width_',num2str(image_width*1e3),'_c_',num2str(c),'.mat']) , 'A_mat','-v7.3')
    end
end

% pre-filtration (is not optimized for large datasets (e.g. high rep num))
filter_f = [filter_min filter_max];      % a band-pass filter (originlay between 100kHz and 8MHz), or maybe try the low pass a little lower (it should be depending on the BW of the transducer)
sigMat_truncated = filter_function(sigMat_truncated, filter_f, datainfo.HWDesc.SamplingFrequency);

%% recon the extracted stack
reconStack = zeros(n, n, runs, Zpos_end-Zpos_begin+1, rep_end-rep_begin+1, wl_end-wl_begin+1);

for run_idx = 1:1
    for slc_idx = 1 : Zpos_end-Zpos_begin+1
        for wl_idx = 1 : wl_end-wl_begin+1
            for rep_idx = 1 : rep_end-rep_begin+1
                % extract current measure (2d) to recon it
                sigMat_current = sigMat_truncated(:, :, run_idx, slc_idx, rep_idx, wl_idx);
                
                switch RECON_METHOD
                    case 'BackProjection'
                        Recon_tmp = backproject_luis(sigMat_current, n, datainfo.HWDesc.Radius, angle_sensor, c, 'full', ts, datainfo.HWDesc.SamplingFrequency, image_width, 0, 0);
                        if NONNEG
                            Recon_tmp = max(Recon_tmp,0);
                        end
                end
                
                if strcmp(RECON_METHOD,'MB_Tik')||strcmp(RECON_METHOD,'MB_TVL1')                    
                    % interpolate raw data to smaller resolution (downsampling)
                    sigMat2 = zeros(length(t), size(sigMat_current, 2));
                    for j = 1:datainfo.HWDesc.NumDetectors
                        sigMat2(:,j) = interp1(ts, sigMat_current(:,j), t);
                    end
                    clear sigMat_current;
                    
                    % reshape interpolated measures to column-major format
                    b_vec = reshape(sigMat2, sizeT*datainfo.HWDesc.NumDetectors, 1);
                    clear sigMat2
                    
                    % do recon
                    tic;    Recon_tmp  = reconstruction(A_mat, b_vec, n, RECON_METHOD, MB_regu, NONNEG);    toc;
                end
                
                reconStack(:, :, run_idx, slc_idx, rep_idx, wl_idx) = Recon_tmp;
                clear Recon_tmp;
                
            end
        end
    end
end

save(strcat(PathNamesigMat, 'recon', RECON_METHOD, '_imSz', num2str(n), '_rep', num2str(rep_begin), 'to', num2str(rep_end)), 'reconStack');

%% show specific im
slcToShow = 1;
wlToShow = 1;
repToShow = 1;
reconIm = reconStack(:, :, 1, slcToShow, repToShow, wlToShow);
if strcmp(RECON_METHOD, 'BackProjection')
    MB_regu = 0;
end
reconIm = addtext(reconIm, n, 1, slcToShow, datainfo.Wavelengths, wlToShow, image_width, c, RECON_METHOD, MB_regu, 0, NONNEG);
figure, imagesc(reconIm); colormap('gray'); colorbar ;



