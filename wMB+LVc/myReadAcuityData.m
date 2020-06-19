% %%%% NOT WORKING!!
% "loading" data is correct (by comparing "sigMat_TOT" coming from this code & from the test script "tests\test_loadnprocessModelBased.m"),
% so most probab the prob is in the recon!!

% %% add paths
% % handHeldPath = 'D:\Dropbox\PA_imaging\msotrecon_handheld';
% handHeldPath = 'C:\Users\marwan.muhammad\Dropbox\PA_imaging\msotrecon_handheld';
% scanFolder = 'D:\DATASETS\Acuity\from Hong\Scan_17 (3 tubes)';
% 
% javaaddpath( [ handHeldPath '\Recon Acuity - Basic\msotbeans.jar' ] ) ;
% %addpath( [ handHeldPath '\reconstruction code' ] ) ;
% addpath( [ handHeldPath '\preprocessing' ] ) ;
% addpath( [ handHeldPath '\reconstruction' ] ) ;
% addpath( [ handHeldPath '\reconstruction\modelbased' ] ) ;
% addpath( [ handHeldPath '\reconstruction\BackProj_MexCore' ] ) ;
% addpath( [ handHeldPath '\reconstruction\ReconPolar' ] ) ;
% addpath( [ handHeldPath '\corelib' ] ) ;
% setenv('KMP_DUPLICATE_LIB_OK', 'TRUE')
% addpath( [ handHeldPath '\corelib\extension' ] ) ;
% addpath( [ handHeldPath '\system_parameters' ] ) ;
% javaaddpath ( [ handHeldPath '\corelib\gui\java_class\xbean.jar' ] ) ;
% javaaddpath ( [ handHeldPath '\corelib\gui\java_class\patbeans.jar' ] ) ;
% javaaddpath ( [ handHeldPath '\corelib\gui\java_class\msotbeans.jar' ] ) ;
% %  javaaddpath ( [ handHeldPath '\corelib\gui\java_class\recon_gui_4_handheld.jar' ] ) ;
% 
% 
% FID = fopen([scanFolder '\Scan_17.bin'], 'r', 'ieee-le.l64' ) ;      % supply instead the .bin file
% % fseek(FID, 0, 'eof') ;
% % frames = ftell(FID);
% % channels = 256;
% % NumFrames = frames / 2 / (channels*2030) ;
% % Frame = NumFrames;
% 
% 
% dataInfoHH = loadMSOTAcuity([scanFolder '\Scan_17.msot']);
% 
% % taken as is from Recon_GUI_Handheld.m when loaded Scan_17 (3 tubes)
% handles.index_selected = 1;
% handles.FileInfo.Frame(handles.index_selected) = 237;
% handles.FileInfo.ScanInfo(handles.index_selected).Data = dataInfoHH.ScanFrames;
% 
% % read sigMat (BUG!!!!   sigMat_raw is all ZEROS!!)
% sigMat_raw = FileRead(FID, handles.FileInfo.Frame(handles.index_selected), 1, handles.FileInfo.ScanInfo(handles.index_selected).Data) ;
% 
% 
% % handles.sigMat_raw = sigMat_raw;
% % preprocess data
% sigMat = PreProcessing(sigMat_raw);         % BPF (equivalent to the 2 lines under from gui code)
% pDetIn = sigMat;
% % processdata = handles.ProcessData(handles.index_selected);
% % sigMat = processdata.preProcessSingleWaveSlide(sigMat_raw);
% % handles.sigMat = sigMat;


%% reading (taken from 'tests\test_loadnprocessModelBased.m')

handHeldPath = 'C:\Users\marwan.muhammad\Dropbox\PA_imaging\msotrecon_handheld\';
javaaddpath([handHeldPath 'corelib\gui\java_class\xbean.jar']);
javaaddpath([handHeldPath 'corelib\gui\java_class\patbeans.jar']);
javaaddpath([handHeldPath 'corelib\gui\java_class\recon.jar']);
javaaddpath([handHeldPath 'Recon Acuity - Basic\msotbeans.jar']) ;

addpath(strcat(handHeldPath,'corelib'));
addpath(strcat(handHeldPath,'corelib\extension\com.itheramedical.msotlib'));
addpath(strcat(handHeldPath,'Recon Acuity - Basic'));
addpath(strcat(handHeldPath,'preprocessing'));
addpath(strcat(handHeldPath,'reconstruction'));
addpath(strcat(handHeldPath,'reconstruction\modelbased'));

scanFolder = 'D:\DATASETS\Acuity\from Hong\Scan_17 (3 tubes)\';
% 'D:\DATASETS\Acuity\(2016-12-15) Study_16 (Murad head)\Scan_5\'
% 'D:\DATASETS\Acuity\from Hong\Scan_17 (3 tubes)\'
metaFile = 'Scan_17.msot';
binFile = 'Scan_17.bin';

datainfo  = loadMSOTAcuity( [scanFolder '\' metaFile] );

imp_resp = zeros( 2030, 1 ) ;
imp_resp( 1015 ) = 1 ;
radiusSensor = .061;
rparams = struct( 'n', 600, ...
                  'proj', 256, ...
                  'time_res', 4, ...
                  'r_sensor', radiusSensor , ...
                  'angle_sensor', ( -3 : -174/255 : -177 ) * pi / 180, ...
                  'c', 1500, ...
                  'limits', radiusSensor + [ -.03 .03 ], ...
                  'imp_resp', imp_resp, ...
                  'fs', datainfo.HWDesc.SamplingFrequency, ...
                  'filter_f', [0.05e6 10e6 ], ...
                  'lsqr_iter', 20 ) ;

SampleShift = 460;
image_width = 0.06; %unit:m;    
channels = rparams.proj ;
wavelengths = datainfo.Wavelengths;   
LaserEnergy = [datainfo.ScanFrames(1:end).LaserEnergy];
slices = datainfo.RepNum;
Nwav = length(wavelengths);
LaserEnergy_cal = LaserEnergy./(max(LaserEnergy));

%% read data
FID  = fopen( [scanFolder binFile], 'r', 'ieee-le.l64' ) ;
sigMat_TOT = zeros ( 2030, channels,slices,Nwav) ;
sigMat_TOTcal = zeros ( 2030, channels,slices,Nwav) ;
for sl = 1:slices
    sl
    for wav = 1:Nwav
        sig_temp = reshape( fread( FID, channels*2030, 'uint16' ), 2030, channels ) ;       % positions the file pointer after the last value read, so no need for indexing
        sigMat_TOT (:, :, sl,  wav) = sig_temp;
        sigMat_TOTcal(:, :, sl,  wav) = sig_temp./LaserEnergy_cal(1,wav+Nwav*(sl-1));
    end
end
fclose(FID);


%% reconstruc images TODO: Bug: Recon=0?
obj = ReconMethods(DeviceInfo('handheld_acuity'));
obj.method = 6;     % 3: MB_tik
obj.nonneg = 0;

switch obj.method
    
    case 1 %'Backprojection_Gael' --> works at highest resolution possible (native resoultion: 1600*1600 (FOV/pxlSz))
        Recon = BackProjectGmex (obj.devInfo.samples, obj.devInfo.proj, 1, double(sigMat_TOT), ...
            double(obj.devInfo.angle_sensor), obj.devInfo.fsample, obj.devInfo.delayNumber, ...
            obj.devInfo.c0, obj.devInfo.r_sensor, obj.devInfo.image_width);
        
        Recon = double(-Recon);
        [x, y] = size(Recon);
        Recon = equalise(x, y, 1, 1, Recon, 2, 6);
        if obj.nonneg
            Recon=max(Recon,0);
        end
        
    case 2 %'Backprojection_Luis'
        %account for delay by adding zeros
        sigMat_TOT=cat(1, zeros([obj.devInfo.delayNumber,size(sigMat_TOT,2), size(sigMat_TOT,3),size(sigMat_TOT,4),size(sigMat_TOT,5),size(sigMat_TOT,6)]),sigMat_TOT);
        sigMat_TOT=new_preprocessing(sigMat_TOT,4e7,[obj.filter_min, obj.filter_max],1,760);
        angle_sensor = obj.devInfo.angle_sensor(end:-1:1);
        Recon = backproject_luis(sigMat_TOT,obj.devInfo.n,obj.devInfo.r_sensor,angle_sensor,obj.devInfo.c0,'full',...
            (0:1/40000000:(obj.devInfo.delayNumber+2030-1)/40000000),...%(0:1/obj.devInfo.fsample:(obj.devInfo.delayNumber+obj.devInfo.samples-1)/obj.devInfo.fsample),...
            obj.devInfo.fsample,obj.devInfo.image_width,0,0);
        if obj.nonneg
            Recon=max(Recon,0);
        end
    case 3 %'MB_Tik'
        %account for delay by adding zeros
        sigMat_TOT=cat(1, zeros([obj.devInfo.delayNumber,size(sigMat_TOT,2), size(sigMat_TOT,3),size(sigMat_TOT,4),size(sigMat_TOT,5),size(sigMat_TOT,6)]),sigMat_TOT);
        
        [A_mat, ts, t, sizeT] = obj.computeOrLoadModelMatrix();
        b_vec = obj.prepareSignalVec_b(sigMat_TOT,ts,t,sizeT);
        
        nn = obj.devInfo.n*obj.devInfo.n;
        FILT = obj.MB_regu*calculate_matrix_highpass_1(obj.devInfo.n);
        L = sparse(FILT);
        cell_mat{1,1} = A_mat;
        cell_mat{2,1} = L;
        A_mat = cell2mat(cell_mat);
        b_vec = [b_vec; zeros(nn,1)];
        if obj.nonneg
            R=nnls_conjgrad_armijo(A_mat,b_vec,zeros(obj.devInfo.n*obj.devInfo.n,1),0.001,5,3);
        else
            R=lsqr_b(A_mat,b_vec,50);
        end
        Recon=reshape(R(:,end),obj.devInfo.n,obj.devInfo.n);
        disp('MB_Tik done!')
    case 4 %'MB_TVL1'
        %account for delay by adding zeros
        sigMat_TOT=cat(1, zeros([obj.devInfo.delayNumber,size(sigMat_TOT,2), size(sigMat_TOT,3),size(sigMat_TOT,4),size(sigMat_TOT,5),size(sigMat_TOT,6)]),sigMat_TOT);
        
        [A_mat, ts, t, sizeT] = obj.computeOrLoadModelMatrix();
        b_vec = obj.prepareSignalVec_b(sigMat_TOT,ts,t,sizeT);
        
        TVWeight = 0e5; 	% Weight for TV penalty
        xfmWeight = obj.MB_regu;	% Weight for Transform L1 penalty
        Itnlim = 50;		% Number of iterations
        Recon  = TVL1_yiyong(A_mat,b_vec,TVWeight,xfmWeight,Itnlim,obj.nonneg);
        disp('MB_TVL1 done!')
    case 5 % 'WaveletPacket'
        %account for delay by adding zeros
        sigMat_TOT=cat(1, zeros([obj.devInfo.delayNumber,size(sigMat_TOT,2), size(sigMat_TOT,3),size(sigMat_TOT,4),size(sigMat_TOT,5),size(sigMat_TOT,6)]),sigMat_TOT);
        
        [A_mat, ts, t, sizeT] = obj.computeOrLoadModelMatrix();
        b_vec = obj.prepareSignalVec_b(sigMat_TOT,ts,t,sizeT);
        
        [Ainv_A,places,BB_part_temp] = obj.computeWaveletMBmatrices(A_mat,sizeT);
        Recon = reconstruction_WP(A_mat,Ainv_A, b_vec,places,BB_part_temp,obj.devInfo.n, t,obj.devInfo.proj,obj.nonneg);
    case 6 %'ReconPolar'
        [systemParams,reconParams,interpolParams] = obj.prepareReconPolarParams();
        
        A_inv = [] ;
        A_stack = [] ;
        sigMat_TOT=cat(1, zeros([obj.devInfo.delayNumber,size(sigMat_TOT,2), size(sigMat_TOT,3),size(sigMat_TOT,4)]),sigMat_TOT);
        [Recon,systemParams,reconParams,interpolParams,Res,A_stack,A_inv]=reconPolar(sigMat_TOT,...
            systemParams,reconParams,interpolParams,A_stack,A_inv);
    otherwise
        error('reconstruction method unknown::Choose from: (1)Backprojection_Gael, (2)Backprojection_Luis, (3)MB_Tik, (4)MB_TVL1, (5)WaveletPacket');
end

% Recon = processdata.resonstructSingleWaveSlide(sigMat_TOT);
figure, imagesc(Recon), colormap('gray'); 