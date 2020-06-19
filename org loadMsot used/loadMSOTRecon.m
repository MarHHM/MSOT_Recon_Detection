function [R wls zpos ts datainfo] = loadMSOTRecon(datainfo,varargin);
%% Recon Image Loader - loads images from a Recon set described in datainfo
%
% @param  datainfo   Data Description as retrieved from loadMSOT
%                    alternatively: Filename of .msot File with Info
% @opt    rNode      Number of the ReconNode to load
%                    in the GUI Dataloader dialog
% @opt    selMat     Selection Matrix based on
%                    datainfo.ReconNode(x).ReconStructure
%
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

%% parameters
progressBar = [];
rNode = 0;
selMat = [];

% in case first argument is not a struct, try to load data from file
if ~isstruct(datainfo)
    if (exist(datainfo,'file'))
        datainfo = loadMSOT(datainfo);
    else
        error('Wrong input parameters, first argument must be data description');
    end
end

% if no reconNode numner is given, use the first
if nargin >= 2
    rNode = varargin{1};
else
    rNode = 1;
end

% if no selMat is given, use all
if nargin >= 3
    selMat = varargin{2};
else
    selMat = datainfo.ReconNode(rNode).ReconStructure;
end

%% Initialise Progress Bar
dlg = [] ;
if( isempty( progressBar ) )
    if( isempty( getappdata( 0, 'loading_dialog' ) ) )
        dlg = com.helmholtz.pat.PATDaqProgressDialog( [], 0 ) ;
        dlg.getContentPane().getComponent(1).setText( 'Hide' ) ;
        dlg.setTitle( 'Load data progress' ) ;
        setappdata( 0, 'loading_dialog', dlg ) ;
    else
        dlg = getappdata( 0, 'loading_dialog' ) ;
    end ;
    progressBar = dlg.getContentPane().getComponent(0) ;
    dlg.setVisible(1) ;
end ;
progressBar.setIndeterminate(1) ;
progressBar.setString( 'Initialising, please wait.' ) ;


%% initialise return array
sel = ~isnan(selMat);   % this will retain structure, marking only not-nans
svec = size(sel);       % retain structural information
selMat = selMat(sel);   % remove nans and make selMat 1D

% image return vector
n = datainfo.ReconNode(rNode).Resolution;
R = nan([n n prod(svec)],'double');
% initial vectors for wavelengths and zpositions and timestamps
wls = zeros(size(selMat));
zpos = zeros(size(selMat));
ts = zeros(size(selMat));

%% Load Images
fname = [datainfo.RealPath '\RECONs\' datainfo.ReconNode(rNode).GUID '.bin'];
[FID str] = fopen(fname,'r','n');
if (FID == -1)
    error(['Error opening binary Recon File ' fname ': ' str]);
    return;
end
% handle progress bar
progressBar.setMaximum(numel(selMat)-1);
progressBar.setValue( 0 ) ;
progressBar.setIndeterminate(0) ;
progressBar.setString( '' ) ;

% load all selected images sequentially
for j = 1:numel(selMat)
    progressBar.setValue( j ) ;
    id = selMat(j);
    
    % copy Meta-Information from Scan Frame
    wls(j) = datainfo.ReconNode(rNode).Frames(id).Wavelength;
    zpos(j) = datainfo.ReconNode(rNode).Frames(id).ZPos;
    ts(j) = datainfo.ReconNode(rNode).Frames(id).RelTime;
    try
        startbyte = (id-1)*n*n*8;    %determine start position (double prec)
        fseek(FID,startbyte,-1);     % move towards position
        R(:,:,j) = fread(FID,[n n],'double');
    catch ex
        warning(['Cannot Read ReconFrame ' num2str(id) ', skipping...']);
    end
   
end
fclose(FID);

%% reshape data to fit requirements
R = reshape(R,[n n svec]);
% ts = reshape(ts,svec);
wls = unique(wls);
zpos = unique(zpos);

dlg.setVisible(0) ;
