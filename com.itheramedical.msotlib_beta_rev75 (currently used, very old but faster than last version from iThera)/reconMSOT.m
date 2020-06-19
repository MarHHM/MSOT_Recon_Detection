function [R wls zpos ts] = reconMSOT(fn,varargin)
% [R wls zpos ts] = reconMSOT(fn,parameters[,verbose]);
% 
% Reconstruction function to reconstruct a selected dataset using the
% supplied parameters
%
% Input Parameters:
% -----------------
% - Filename: Either filename of .msot file or datainfo object
% - Parameters: (struct)
%   n            -> Resolution (default: 200)
%   proj         -> Number of Projections (default: 2*number of detectors)
%   c            -> Speed of Sound (either absolute or offset, default: 0)
%   roi          -> Region of Interest (in cm, default: 20)
%   image_select -> Reconstruction algorithm (*direct, model_lin, wavelet)
%   filter_f     -> Bandpass frequencies (default: 50kHz - 7Mhz)
%   iter         -> LSQR Iterations for 'model_lin'
%   timeres      -> Timeresolution for 'model_lin' and 'wavelet'
%   selMat       -> Selection Matrix based on datainfo.ScanStructure
% - Verbose: Print verbose logging information (true/false)
%
% Return Values:
% --------------
% @return R          An array containing the reconstructed images in 8D
%                    1) x
%                    2) y
%                    3) z (only if 3D detector, currently unused)
%                    4) RUN (repetition of the whole acquisition step, time
%                       dimension (timepoints as separate vector)
%                    5) z-position (if 2D system with translation stage)
%                    6) REPETITION (currently unused)
%                    7) Wavelength (see wls vector)
%                    8) Individual Images (only if not averaged)
% @return wls        Wavelength Vector
% @return zpos       Array of z-stage positions
% @return ts         Multidimensional array of timestamps in s, same 
%                    dimension as R

%% Input Parameters
% Filename
if (isstruct(fn))
    datainfo = fn;
    fn = datainfo.XMLFileName;
else
    if (~exist(fn,'file'))
      error(['MSOT File not found: ' fn]);
    else 
        datainfo = loadMSOT(fn);
    end
end

% Parameter default set
par = getWLdefaults;
par.n = 200;
par.proj = datainfo.HWDesc.NumDetectors*2-1;
par.r_sensor = datainfo.HWDesc.Radius;
par.c = 0;
par.roi = 20e-3;
par.image_select = 'direct';
par.filter_f = [10 0]*1e3;
par.iter = 50;
par.timeres = 3;
par.n_angles = par.n*2;
par.limitSensors = [];
par.selMat = datainfo.ScanStructure;
par.Amat_dir = 'E:\_Amatrices';
par.f_HPFOrder = 4;
par.f_HPFRipple = 0.01;
par.f_HPFzp = 1;
par.save = 0;       % 2: save in matlab file



% Copy parameters from input struct
if numel(varargin) >= 1
    cpar = varargin{1};
    % transfer all fields to parameter array
    fx = fieldnames(cpar);
    for j = 1:numel(fx)
        par = setfield(par,fx{j},getfield(cpar,fx{j}));
    end
    clear cpar j fx;
end

if numel(varargin) >= 2
    verbose = varargin{2};
end
    

%% derived parameters
if (par.c < 1000)
    T = datainfo.AverageTemperature;
    par.c = round(1.402385 * 1e3 + 5.038813 * T - 5.799136 * 1e-2 * T^2 + 3.287156 * 1e-4 * T^3 - 1.398845 * 1e-6 * T^4 + 2.787860 * 1e-9 * T^5 + par.c); 
else 
    par.c = round(par.c);
end


% put selMat to 1D 
selMat = reshape(par.selMat,numel(par.selMat),1);

%% load signal Data
% to be safe, cast as double if usePower is 0...
sigMat = double(loadMSOTSignals(datainfo,selMat,par));
% ts = datainfo.RelTime(selMat);


%% preprocess signals
[sigMat, t, angle_sensor] = preprocessSignals(sigMat,datainfo,par);


%% preparation
% backprojection
if strcmp(par.image_select,'direct')
    r_sensor = [ cos( angle_sensor' ) sin( angle_sensor' ) zeros( size( angle_sensor' ) ) ] * par.r_sensor ;
    fs = datainfo.HWDesc.SamplingFrequency;
% model based
elseif (strcmp(par.image_select,'model_lin'))
    A_mat = getAMat(par,angle_sensor,t);
% wavelet model based
elseif (strcmp(par.image_select,'wavelet'))
    startupcl('vendor','amd','type', 'cpu');
    wl_filename = getWLfilename(par);
    if (~exist(wl_filename,'file'))
        A_mat = getAMat(par,angle_sensor,t);
        wl_filename = invertWL(A_mat,par);   
    end   
else
    error(['Invalid Reconstruction selected: ' par.image_select ]);
end

%% reconstruct
R = zeros(par.n,par.n,size(selMat,1));
coreNum = feature( 'numCores' ) ;
% progress bar
wbar = waitbar(0,'Estimated Time Left: Unknown','Name','Reconstruction','CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');

for jj = 1:size(selMat,1)
    id = selMat(jj);
    tic;
    
    if getappdata(wbar,'canceling')
        close(wbar);
        delete(wbar);
        error('User cancelled the reconstruction');
    end


    % Reconstruction (lsqr or backprojection)    
    if strcmp(par.image_select,'direct')
%         fprintf('Backprojection: Image %i of %i (Frame %i)\n',jj,size(selMat,1),id);    
        R(:,:,jj) = backproject( sigMat(:,:,jj), par.n, r_sensor(:, 1), r_sensor(:, 2), r_sensor(:, 3), par.c, 1, t', fs, par.roi, coreNum ) ;
    elseif (strcmp(par.image_select,'model_lin'))   
%         fprintf('Model Based: Image %i of %i (Frame %i), %i iterations\n',jj,size(selMat,1),id,par.iter);    
        [recon, flag] = lsqr(A_mat, ...
            reshape(sigMat(:,:,jj),size(sigMat,1)*size(sigMat,2),1),...
            1e-6, par.iter);
        R(:,:,jj) = -reshape(recon,par.n,par.n);
        clear recon;
    end
    
    % time estimation
    perimg = toc;
    est = round(perimg*(size(selMat,1)-jj)); est_unit = 's';
    if (est > 120) est = round(est / 60); est_unit = 'min'; end;    
    wbar = waitbar(jj/size(selMat,1),wbar,sprintf( 'Estimated Time Left: %i%s', est, est_unit )) ;
    
    
    % Meta information
    wls(jj) = datainfo.ScanFrames(id).Wavelength;
    zpos(jj) = datainfo.ScanFrames(id).ZPos;
    ts(jj) = datainfo.ScanFrames(id).RelTime;
end
close(wbar);delete(wbar);

if (strcmp(par.image_select,'wavelet'))   
    fprintf('Wavelet Reconstruction: Batch processing %i images...\n',size(selMat,1));    
    [Lo_D,Hi_D] = wfilters(par.wl_name);
    SetWaveletFilters_cl(Lo_D, Hi_D);

    ok = loadWLModel_cl(wl_filename);
    sigMat = reshape(sigMat,size(sigMat,1)*size(sigMat,2),size(sigMat,3));

    tic
    R = ReconWLcl(sigMat, numel(t),par.proj,par.depth_proj,par.n);
    recont = toc;
    fprintf('  Execution time: %.1fs (%.3fs per image)\n',recont,recont/size(sigMat,2));
    
    cleanupcl;
end   

R = reshape(R,[par.n par.n size(par.selMat)]);
ts = reshape(ts,size(par.selMat));
wls = unique(wls);
zpos = unique(zpos);

if (par.save == 2)
    savefile = [datainfo.FolderName '\' datainfo.FriendlyName ...
        '_' par.image_select...
        '_roi' num2str(par.roi*1e3,'%.1f')...
        '_n' num2str(par.n)...
        '_c' num2str(par.c)...
        '_proj' num2str(par.proj)...
        '_hp' num2str(par.filter_f(1)*1e-3,'%.0f') 'khz'...
        '.mat'];
    fprintf('Saving Recon to %s\n',savefile);
    save(savefile,'-v7.3','R','ts','wls','zpos','datainfo','par');
end
