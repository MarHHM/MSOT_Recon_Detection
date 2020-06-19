function varargout = evalFilters(varargin)
% EVALFILTERS MATLAB code for evalFilters.fig
%      EVALFILTERS, by itself, creates a new EVALFILTERS or raises the existing
%      singleton*.
%
%      H = EVALFILTERS returns the handle to a new EVALFILTERS or the handle to
%      the existing singleton*.
%
%      EVALFILTERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EVALFILTERS.M with the given input arguments.
%
%      EVALFILTERS('Property','Value',...) creates a new EVALFILTERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before evalFilters_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to evalFilters_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help evalFilters

% Last Modified by GUIDE v2.5 14-Feb-2013 17:19:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @evalFilters_OpeningFcn, ...
                   'gui_OutputFcn',  @evalFilters_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before evalFilters is made visible.
function evalFilters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to evalFilters (see VARARGIN)

% Choose default command line output for evalFilters
handles.output = hObject;

if numel(varargin) >= 1
    handles.sig = varargin{1};
else
    error('Insufficient arguments');
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes evalFilters wait for user response (see UIRESUME)
% uiwait(handles.fig);


% --- Outputs from this function are returned to the command line.
function varargout = evalFilters_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;






function sl_sos_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function cb_type_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function sl_orderLP_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function sl_lcf_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function sl_hcf_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function sl_orderHP_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- change c
function sl_sos_Callback(hObject, eventdata, handles)
set(handles.l_sos,'String',num2str(get(hObject,'Value'),'%.1f'));

% --- change zero phase
function cb_zerophase_Callback(hObject, eventdata, handles)

% --- change type
function cb_type_Callback(hObject, eventdata, handles)

% --- change order.
function sl_orderLP_Callback(hObject, eventdata, handles)
set(handles.l_orderLP,'String',num2str(get(hObject,'Value'),'%d'));

% --- change lcf
function sl_lcf_Callback(hObject, eventdata, handles)
set(handles.l_lcf,'String',num2str(get(hObject,'Value'),'%.1f kHz'));

% --- change order HP
function sl_orderHP_Callback(hObject, eventdata, handles)
set(handles.l_orderHP,'String',num2str(get(hObject,'Value'),'%d'));

% --- change hcf
function sl_hcf_Callback(hObject, eventdata, handles)
set(handles.l_hcf,'String',num2str(get(hObject,'Value'),'%.1f Mhz'));

function b_go_Callback(hObject, eventdata, handles)
    zerophase = get(handles.cb_zerophase,'Value');
    n = 200;
    proj = 256;
    fs=4e7;
    load E:\transducers\tdata_256A101.mat;
    r_sensor = 0.0405;
    c = get(handles.sl_sos,'Value');
    limits = [r_sensor-0.0125 r_sensor+0.0125];
%     filter_f = [get(handles.sl_lcf,'Value')*1e3 get(handles.sl_hcf,'Value')*1e6 ];
    filter_f = [get(handles.sl_lcf,'Value')*1e3 get(handles.sl_hcf,'Value')*1e6 ];
    if (~zerophase)
        filter_f = [filter_f 1];
    end

    global filter_struct ;
    f_LPF = filter_f(2) ;
    filter_struct.passband(2) = f_LPF ;
    filter_struct.a_LPF = [] ;
    filter_struct.b_LPF = [] ;
    if (f_LPF)
        [b_LPF,a_LPF] = cheby1( get(handles.sl_orderLP,'Value'), .01, 2 * f_LPF/fs/3 * .9 ) ;
        filter_struct.a_LPF = a_LPF ;
        filter_struct.b_LPF = b_LPF ;
    end

    f_HPF = filter_f(1) ;
    filter_struct.passband(1) = filter_f(1) ;
    filter_struct.a_HPF = [];
    filter_struct.b_HPF = [] ;
    if (f_HPF)
        [b_HPF,a_HPF] = cheby1( get(handles.sl_orderHP,'Value'), .01, 2 * f_HPF/fs/3 * 1.46, 'high' ) ;
        filter_struct.a_HPF = a_HPF ;
        filter_struct.b_HPF = b_HPF ;
    end

    [ BP AM sigf ] = evalFilters_reconWrapper([],handles.sig,[],n,proj,r_sensor,angle_sensor,c,...
        'direct',fs,limits,3,filter_f,0);
    axes(handles.ax_display);
    imagesc(BP);
    caxis([0 max(BP(:))]);
    axis image;
    colorbar;
    
%     figure;
    axes(handles.ax_sino);
    imagesc(sigf);
    xlabel('channels');
    ylabel('samples (time)');
    
    axes(handles.ax_plot);
    plot(sigf(:,[1:16:256]));
    xlabel('samples (time)');