function varargout = MSOT_RECON_DETECTION(varargin)
% MSOT_RECON_DETECTION MATLAB code for MSOT_RECON_DETECTION.fig
%      MSOT_RECON_DETECTION, by itself, creates a new MSOT_RECON_DETECTION or raises the existing
%      singleton*.
%
%      H = MSOT_RECON_DETECTION returns the handle to a new MSOT_RECON_DETECTION or the handle to
%      the existing singleton*.
%
%      MSOT_RECON_DETECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MSOT_RECON_DETECTION.M with the given input arguments.
%
%      MSOT_RECON_DETECTION('Property','Value',...) creates a new MSOT_RECON_DETECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MSOT_RECON_DETECTION_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MSOT_RECON_DETECTION_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MSOT_RECON_DETECTION

% Last Modified by GUIDE v2.5 16-Nov-2016 18:46:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @MSOT_RECON_DETECTION_OpeningFcn, ...
    'gui_OutputFcn',  @MSOT_RECON_DETECTION_OutputFcn, ...
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


% --- Executes just before MSOT_RECON_DETECTION is made visible.
function MSOT_RECON_DETECTION_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MSOT_RECON_DETECTION (see VARARGIN)

% Choose default command line output for MSOT_RECON_DETECTION
handles.output = hObject;
cd('C:\Users\ibn_a\OneDrive\PA_imaging\wMB+LVc\MSOT_Recon_Detection\');
addpath('unmixing code');
addpath('reconstruction code');
addpath('MB matrices');
axis off; axis equal;%colorbar ;

handles = guidata(hObject);
handles.output = 1;
clear handles.listbox_sliceSelect;
clear handles.listbox2;
clear handles.listbox3;
guidata(hObject, handles);

handles = guidata(hObject);
guidata(hObject, handles);

set(handles.pushbuttonSpec,'enable','off');
set(handles.pushbuttonAnalyze,'enable','off');
set(handles.pushbutton6,'enable','off');

set(handles.pushbuttonRecon,'enable','off');
set(handles.pushbuttonReconCurrSlice,'enable','off');

set(handles.pushbuttonSelectROI,'enable','off');
set(handles.pushbuttonExpSpec,'enable','off');

handles.recon_method='MB_Tik';
handles.n=200;
handles.image_width=20e-3;
handles.regu=20;
handles.Recon=zeros(handles.n,handles.n);
handles.WP_global_threshold=1;
handles.noneg=1;
handles.filter_min=0.02e6;
handles.filter_max=7e6;

% Update handles structure

set(gcf, 'WindowButtonDownFcn', @getMousePositionOnImage);
pan off % Panning will interfere with this code

guidata(hObject, handles);

% UIWAIT makes MSOT_RECON_DETECTION wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MSOT_RECON_DETECTION_OutputFcn(hObject, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in listbox_sliceSelect.
function listbox_sliceSelect_Callback(hObject, ~, handles)
% hObject    handle to listbox_sliceSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_sliceSelect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_sliceSelect
slice = get(handles.listbox_sliceSelect,'value');
handles.CurrSlice = slice;
guidata(hObject, handles);
if sum(strcmp(fieldnames(handles), 'Analysis')) == 0
    R = squeeze(handles.Recon(:,:,handles.CurrRun,handles.CurrSlice,1,handles.CurrWav));
    R=addtext(R,handles.n,handles.CurrRun,handles.CurrSlice,handles.datainfo.Wavelengths,handles.CurrWav,handles.image_width,handles.c,handles.recon_method,handles.MB_regu,handles.truncated,handles.noneg);
    axes(handles.axes1);
    imagesc(R(:,:)); colormap('gray'); axis off; axis equal;colorbar ;
else
    axes(handles.axes1);
    image(uint8(255*squeeze(handles.Analysis(:,:,:,handles.CurrRun,handles.CurrSlice)))); axis off; axis equal;colorbar ;
end

% --- Executes during object creation, after setting all properties.
function listbox_sliceSelect_CreateFcn(hObject, ~, handles)
% hObject    handle to listbox_sliceSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2
run = get(handles.listbox2,'value');
handles.CurrRun = run;
guidata(hObject, handles);
if sum(strcmp(fieldnames(handles), 'Analysis')) == 0
    R = squeeze(handles.Recon(:,:,handles.CurrRun,handles.CurrSlice,1,handles.CurrWav));
    R=addtext(R,handles.n,handles.CurrRun,handles.CurrSlice,handles.datainfo.Wavelengths,handles.CurrWav,handles.image_width,handles.c,handles.recon_method,handles.MB_regu,handles.truncated,handles.noneg);
    axes(handles.axes1);
    imagesc(R(:,:)); colormap('gray'); axis off; axis equal;colorbar ;
else
    axes(handles.axes1);
    image(squeeze(uint8(255*handles.Analysis(:,:,:,handles.CurrRun,handles.CurrSlice)))); axis off; axis equal;colorbar ;
end

% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in listbox3.
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox3
wav = get(handles.listbox3,'value');
handles.CurrWav = wav;
guidata(hObject, handles);
if sum(strcmp(fieldnames(handles), 'Analysis')) == 0
    R = squeeze(handles.Recon(:,:,handles.CurrRun,handles.CurrSlice,1,handles.CurrWav));
    R=addtext(R,handles.n,handles.CurrRun,handles.CurrSlice,handles.datainfo.Wavelengths,handles.CurrWav,handles.image_width,handles.c,handles.recon_method,handles.MB_regu,handles.truncated,handles.noneg);
    axes(handles.axes1);
    imagesc(R(:,:)); colormap('gray'); axis off; axis equal;colorbar ;
else
    axes(handles.axes1);
    image(uint8(255*squeeze(handles.Analysis(:,:,:,handles.CurrRun,handles.CurrSlice)))); axis off; axis equal;colorbar ;
end

% --- Executes during object creation, after setting all properties.
function listbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
algo = get(handles.popupmenu1, 'value');
algo
if algo == 3 | algo == 4
    set(handles.popupmenu4,'String','no thresh');
else
    set(handles.popupmenu4,'String',{'no thresh','PFA 2.5e-3 (100 pix)', 'PFA 5e-4 (20 pix)','PFA 1.25e-4 (5 pix)'});
end



% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
% Get value of popup
selectedIndex = get(handles.popupmenu2, 'value');
disp(selectedIndex)
% Take action based upon selection
if selectedIndex  == 1 & sum(strcmp(fieldnames(handles), 'spectra')) == 0
    set(handles.pushbuttonAnalyze,'enable','off');
elseif sum(strcmp(fieldnames(handles), 'datainfo')) == 0
    warndlg('Please load MSOT data before providing spectra','Warning!');
    set(handles.popupmenu2,'value',1);
else
    switch selectedIndex
        %ICG
        case 2
            folder_name = 'spectra\SpectralSpecifications_iCG';
            handles.spectra = LoadSpectra(folder_name,handles.datainfo.Wavelengths);
            guidata(hObject, handles);
            %Alexa 750
        case 3
            folder_name = 'spectra\SpectralSpecifications_AF750';
            handles.spectra = LoadSpectra(folder_name,handles.datainfo.Wavelengths);
            guidata(hObject, handles);
            %Alexa 790
        case 4
            folder_name = 'spectra\SpectralSpecifications_AF790';
            handles.spectra = LoadSpectra(folder_name,handles.datainfo.Wavelengths);
            guidata(hObject, handles);
            %IntSense
        case 5
            folder_name = 'spectra\SpectralSpecifications_IntSense';
            handles.spectra = LoadSpectra(folder_name,handles.datainfo.Wavelengths);
            guidata(hObject, handles);
        case 6
            folder_name = 'spectra\Spectral_Specifications_DiR';
            handles.spectra = LoadSpectra(folder_name,handles.datainfo.Wavelengths);
            guidata(hObject, handles);
        case 7
            set(handles.pushbuttonSpec,'enable','on');
            guidata(hObject, handles);
    end
    set(handles.pushbuttonAnalyze,'enable','on');
    axes(handles.axes2);
    plot(handles.datainfo.Wavelengths,handles.spectra(1,:));
    
    grid on;
    axis tight;
end

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure wguide
ith handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on button press in pushbutton_LoadsigMat.
function pushbutton_LoadsigMat_Callback(hObject, ~, handles)
% hObject    handle to pushbutton_LoadsigMat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileNamesigMat,PathNamesigMat,~] = uigetfile({'*.msot'});

if FileNamesigMat
    if FileNamesigMat(end-3:end) == "msot"
        [datainfo] = loadMSOT([PathNamesigMat '\' FileNamesigMat]);
        [sigMat] = loadMSOTSignals(datainfo);
        handles.sigMat=sigMat;
        handles.datainfo = datainfo;
        handles.PathNamesigMat=PathNamesigMat;
        %         handles.datainfo.Wavelengths
        clear datainfo sigMat
        clear handles.listbox_sliceSelect;
        clear handles.listbox2;
        clear handles.listbox3;
        set(handles.listbox_sliceSelect,'string',1:length(handles.datainfo.ZPositions));
        set(handles.listbox2,'string',1:handles.datainfo.RunNum);
        set(handles.listbox3,'string',handles.datainfo.Wavelengths);
        set(handles.text9,'String',handles.datainfo.Name);
%         T = handles.datainfo.AverageTemperature;
%         handles.c= 12+round(1.402385 * 1e3 + 5.038813 * T - 5.799136 * 1e-2 * T^2 + 3.287156 * 1e-4 * T^3 - 1.398845 * 1e-6 * T^4 + 2.787860 * 1e-9 * T^5 );
        handles.c = 1465;       % for heavy water (Qutaiba\Study_4)
        set(handles.text18,'String',handles.c);
        
        handles.CurrSlice = 1;
        handles.CurrRun = 1;
        handles.CurrWav = 1;
        
        set(handles.pushbuttonRecon,'enable','on');
        set(handles.pushbuttonReconCurrSlice,'enable','on');
        
        init_folder_name = pwd;
        sub_folder_name = 'MB matrices';
        folder_name = [init_folder_name '\' sub_folder_name];
        
        
        
        set(handles.edit7,'String',handles.filter_min/1e6);
        set(handles.edit8,'String',handles.filter_max/1e6);
        
        set(handles.slider7,'Max', 1580);
        set(handles.slider7,'Min', 1420);
        set(handles.slider7,'Value',handles.c);
        set(handles.slider7,'SliderStep',[3/39,5/39]);
        
        set(handles.slider8,'Max',handles.regu+20);
        set(handles.slider8,'Min',handles.regu-20);
        set(handles.slider8,'Value',handles.regu);
        set(handles.slider8,'SliderStep',[1/39,3/39]);
        
        handles.truncated = 0.1*(1.1^(handles.regu-20));
        handles.MB_regu=1e6*(2^((handles.regu-20)));
        
        set(handles.text16,'String',handles.MB_regu);
        
        try
            handles.angle_sensor = handles.datainfo.HWDesc.StartAngle : handles.datainfo.HWDesc.StepAngle : handles.datainfo.HWDesc.EndAngle;
        catch     % error due to different MSOT format - check if new reader is required from iThera
            handles.datainfo.HWDesc.Radius = 0.0405;
            handles.datainfo.HWDesc.StartAngle = -0.7762;
            handles.datainfo.HWDesc.EndAngle = 3.9178;
            handles.datainfo.HWDesc.NumDetectors = 256;
            handles.datainfo.HWDesc.StepAngle = (handles.datainfo.HWDesc.EndAngle - handles.datainfo.HWDesc.StartAngle)/(handles.datainfo.HWDesc.NumDetectors-1);
            handles.angle_sensor = handles.datainfo.HWDesc.StartAngle : handles.datainfo.HWDesc.StepAngle : handles.datainfo.HWDesc.EndAngle;
        end
        
        handles.ts = 0:1/handles.datainfo.HWDesc.SamplingFrequency:(handles.datainfo.MeasurementDesc.RecordLength-1)/handles.datainfo.HWDesc.SamplingFrequency; % sampling instants
        
        if strcmp(handles.recon_method,'MB_Tik')||strcmp(handles.recon_method,'MB_TVL1')||strcmp(handles.recon_method,'WaveletPacket')
            
            handles.n_angles = 2*handles.n;                                     % number of points for discretizing the curve
            limits(1) = handles.datainfo.HWDesc.Radius-(handles.image_width)*sqrt(2)/2;       % limits for the signal
            limits(2) = handles.datainfo.HWDesc.Radius+(handles.image_width)*sqrt(2)/2;       % limits for the signal
            handles.time_res = 2;                                       % time resolution for model-based
            dx = handles.image_width/handles.n;                                 % increment in x
            dt = dx/(handles.time_res*handles.c);                              % increment in t employed to make the model-based reconstruction
            fac = handles.datainfo.HWDesc.SamplingFrequency/handles.c;
            pos_start = max(1,int32((limits(1))*fac));
            pos_end = min(handles.datainfo.MeasurementDesc.RecordLength,int32((limits(2))*fac));
            handles.t = handles.ts(pos_start):dt:handles.ts(pos_end);           % downsampled (& cut) time vector  (less than ts)
            handles.sizeT = length(handles.t);
            
            % check for model matrix if saved, otherwise build it
            if exist([folder_name '\' 'A_mat_t_res_',num2str(handles.time_res),'_',num2str(handles.n),'x',num2str(handles.n),'_width_',num2str(handles.image_width*1e3),'_c_',num2str(handles.c),'.mat'], 'file')
                load([folder_name '\' 'A_mat_t_res_',num2str(handles.time_res),'_',num2str(handles.n),'x',num2str(handles.n),'_width_',num2str(handles.image_width*1e3),'_c_',num2str(handles.c),'.mat'])
            else
                A_mat = Calculate_MatrixMB_Luis(handles.c,handles.n,handles.image_width,handles.t,handles.datainfo.HWDesc.Radius,handles.angle_sensor,handles.n_angles);
                if ~exist(folder_name, 'dir')
                    mkdir(folder_name);
                end
                save(([folder_name '\' 'A_mat_t_res_',num2str(handles.time_res),'_',num2str(handles.n),'x',num2str(handles.n),'_width_',num2str(handles.image_width*1e3),'_c_',num2str(handles.c),'.mat']) , 'A_mat','-v7.3')
            end
        end
        
        if strcmp(handles.recon_method,'WaveletPacket')
            wl_name='db6';
            depth_proj=2;
            depth_im=2;
            Sig=zeros(handles.n,handles.n);
            over_fact=ceil(1.4*size(A_mat,1)/size(A_mat,2));
            [temp,L4]=wavedec2(Sig,depth_im,wl_name);
            max_s=0;
            
            if exist([folder_name '\' 'Ainv_A_GT_',num2str(handles.WP_global_threshold),'_',num2str(handles.truncated),'_t_res_',num2str(handles.time_res),'_',num2str(handles.n),'x',num2str(handles.n),'_width_',num2str(handles.image_width*1e3),'_c_',num2str(handles.c),'.mat']) &&exist (['places_',num2str(handles.n),'x',num2str(handles.n),'.mat'])&& exist (['BB_part_temp_',num2str(handles.n),'x',num2str(handles.n),'.mat'])
                
                load([folder_name '\' 'Ainv_A_GT_',num2str(handles.WP_global_threshold),'_',num2str(handles.truncated),'_t_res_',num2str(handles.time_res),'_',num2str(handles.n),'x',num2str(handles.n),'_width_',num2str(handles.image_width*1e3),'_c_',num2str(handles.c),'.mat'])
                load ([folder_name '\' 'places_',num2str(handles.n),'x',num2str(handles.n),'.mat'])
                load ([folder_name '\' 'BB_part_temp_',num2str(handles.n),'x',num2str(handles.n),'.mat'])
            else
                for kk=1:4^(depth_im)
                    kk
                    base=return_base_arb_lev_even_db(kk,depth_im,wl_name) ;
                    BB_part=calc_B_part_arb_lev(base,handles.n,handles.n,L4,depth_im);
                    places_2_check=find_important_places2(BB_part,0.2);
                    BB_part_temp{kk}=sparse(BB_part(:,places_2_check));
                    AA_temp2=Mat_wlet_dec_arb_lev_even_db(A_mat*BB_part_temp{kk},handles.sizeT,handles.datainfo.HWDesc.NumDetectors,depth_proj,wl_name);
                    places{kk}=choose_time_proj_low_mem2(AA_temp2,over_fact,size(BB_part_temp{kk},2));
                    A_part=AA_temp2(places{kk},:);
                    [u,s,v]=svd(full(A_part),'econ');
                    
                    if handles.WP_global_threshold==1
                        switch kk
                            case 1
                                max_s=max(s(:));
                        end
                    elseif handles.WP_global_threshold==0
                        max_s=max(s(:));
                    end
                    
                    Ainv_A{kk}=TSVD_YYH(u,s,v,handles.truncated,max_s);
                    clear u s v
                end
                
                save (([folder_name '\' 'Ainv_A_GT_',num2str(handles.WP_global_threshold),'_',num2str(handles.truncated),'_t_res_',num2str(handles.time_res),'_',num2str(handles.n),'x',num2str(handles.n),'_width_',num2str(handles.image_width*1e3),'_c_',num2str(handles.c),'.mat']) ,'Ainv_A','-v7.3')
                save (([folder_name '\' 'places_',num2str(handles.n),'x',num2str(handles.n),'.mat']) ,'places')
                save (([folder_name '\' 'BB_part_temp_',num2str(handles.n),'x',num2str(handles.n),'.mat']) ,'BB_part_temp')
                
                clear  AA_temp2 A_part BB_part places_2_check
                
            end
        end
        
        avgs = handles.datainfo.MeasurementDesc.Averages;
        runs = handles.datainfo.RunNum;
        repetitions = handles.datainfo.RepNum;
        slices = numel(handles.datainfo.ZPositions);
        lambdas = numel(handles.datainfo.Wavelengths);
        handles.Recon=zeros(handles.n,handles.n,runs,slices,repetitions,lambdas);
        
        %% extract a slice to test recon with
        for run_idx=handles.CurrRun:handles.CurrRun
            for slc = handles.CurrSlice:handles.CurrSlice
                for wvl = handles.CurrWav:handles.CurrWav
                    for rep = 1:repetitions
%                     for rep = 1:min(2,repetitions)    % Mar - 2-12-2016
                        %% select signals depending on the number of wavelengths
                        if(lambdas == 1 && avgs == 1)
                            sigMat = squeeze(handles.sigMat(:, :, run_idx, slc )); % only 5 dimensions if single wavelength
                        end
                        if (lambdas == 1 && avgs ~= 1)
                            sigMat = squeeze(handles.sigMat(:, :, run_idx, slc, :)); % 6 dimensions if multi-wavelength
                        end
                        if (lambdas ~= 1)
                            sigMat = squeeze(handles.sigMat(:, :, run_idx, slc, rep, wvl)); % 6 dimensions if multi-wavelength
                        end
                        
                        filter_f = [handles.filter_min handles.filter_max];%[0.04e6 8e6];      % a band-pass filter (originlay between 100kHz and 8MHz), or maybe try the low pass a little lower
                        sigMat = filter_function(sigMat, filter_f, handles.datainfo.HWDesc.SamplingFrequency);
                        
                        switch handles.recon_method
                            case 'BackProjection'
                                Recon_tmp = backproject_luis(sigMat,handles.n,handles.datainfo.HWDesc.Radius,handles.angle_sensor,handles.c,'full',handles.ts,handles.datainfo.HWDesc.SamplingFrequency,handles.image_width,0,0);
                                if handles.noneg
                                    Recon_tmp=max(Recon_tmp,0);
                                end
                        end
                        
                        if strcmp(handles.recon_method,'MB_Tik')||strcmp(handles.recon_method,'MB_TVL1')||strcmp(handles.recon_method,'WaveletPacket')
                            for j=1:handles.datainfo.HWDesc.NumDetectors
                                sigMat2(:,j) = interp1(handles.ts,sigMat(:,j),handles.t);       % interpolate raw data to smaller resolution (downsampling)
                            end
                            clear sigMat;
                            b_vec = reshape(sigMat2, handles.sizeT*handles.datainfo.HWDesc.NumDetectors,1);           % reshape raw data matrix to column-major format
                            clear sigMat2
                        end
                        
                        if strcmp(handles.recon_method,'MB_Tik')||strcmp(handles.recon_method,'MB_TVL1')
                            Recon_tmp  = reconstruction( A_mat, b_vec,handles.n,handles.recon_method,handles.MB_regu,handles.noneg);
                        end
                        if strcmp(handles.recon_method,'WaveletPacket')
                            Recon_tmp  =reconstruction_WP(A_mat,Ainv_A, b_vec,places,BB_part_temp,handles.n,handles.t,handles.datainfo.HWDesc.NumDetectors,handles.noneg);
                        end
                        
                        handles.Recon(:, :, run_idx, slc, rep, wvl) = Recon_tmp;
                        clear Recon_tmp;
                        
                    end
                end
            end
        end
        %handles.Recon= Recon;
        % shows only the first rep
        R = squeeze(handles.Recon(:, :, handles.CurrRun, handles.CurrSlice, 1, handles.CurrWav));
        R = addtext(R,handles.n,handles.CurrRun,handles.CurrSlice,handles.datainfo.Wavelengths,handles.CurrWav,handles.image_width,handles.c,handles.recon_method,handles.MB_regu,handles.truncated,handles.noneg);
        axes(handles.axes1);
        imagesc(R(:,:)); colormap('gray'); axis off; axis equal;colorbar ;
        guidata(hObject, handles);
    end
else
    return;
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonCalSpec.
function pushbuttonCalSpec_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCalSpec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%here copy code from eigenspectra.

if sum(strcmp(fieldnames(handles), 'datainfo')) >0
    if sum(strcmp(fieldnames(handles.datainfo), 'MeasurementDesc'))>0
        if sum(strcmp(fieldnames(handles.datainfo.MeasurementDesc), 'WaterAbsorptionCoeff'))>0
            
            abs_energy = handles.datainfo.MeasurementDesc.WaterAbsorptionCoeff;
            Tissue_Fl = ones(1,size(handles.Recon,6)).*exp(-abs_energy.*0.04)';
            Tissue_Fl = Tissue_Fl./norm(Tissue_Fl,2);
            for i=1:size(handles.Recon,6)
                handles.Recon(:,:,:,:,:,i) = handles.Recon(:,:,:,:,:,i)./Tissue_Fl(i);
            end
            MSOT_cal = linspace(1,0.85,length(handles.datainfo.Wavelengths));
            for i=1:size(handles.Recon,6)
                handles.Recon(:,:,:,:,:,i) = handles.Recon(:,:,:,:,:,i)./MSOT_cal(i);
            end
        end
    end
end
set(handles.pushbuttonCalSpec,'enable','off');
guidata(hObject, handles);




% --- Executes on button press in pushbuttonSpec.
function pushbuttonSpec_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSpec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileNameSpec,PathNameSpec,FilterIndexSpec] = uigetfile({'*.xls'});
if FileNameSpec
    target_spec = xlsread([PathNameSpec '\' FileNameSpec]);
    target_spec(target_spec<0) = 0;
    save('spectra\SpectralSpecifications_New\agent_target_spec','target_spec');
    folder_name = 'spectra\SpectralSpecifications_New';
    handles.spectra = LoadSpectra(folder_name,handles.datainfo.Wavelengths);
    guidata(hObject, handles);
    set(handles.pushbuttonAnalyze,'enable','on');
    
    axes(handles.axes2);
    plot(handles.datainfo.Wavelengths,handles.spectra(1,:));
    grid on;
    axis tight;
end

% --- Executes on button press in pushbuttonAnalyze.
function pushbuttonAnalyze_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to pushbuttonAnalyze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
algorithm = get(handles.popupmenu1,'value')
thresh_sel = get(handles.popupmenu4,'value')
vis = get(handles.popupmenu3,'value')
spec = get(handles.popupmenu2,'value')

switch algorithm
    case 2
        method = 'Simple';
        switch thresh_sel
            case 1
                thresh = 0;
            case 2
                %PFA = 0.0025 / 100 pixels per 200x200 image
                thresh = 4.3662e-04;
            case 3
                %PFA = 5e-04 / 20 pixels per 200x200 image
                thresh = 9.8935e-04;
            case 4
                %PFA = 1.25e-04 / 5 pixels per 200x200 image
                thresh = 0.0019;
        end
    case 1
        method = 'QL_shrinkage_adaptive';
        switch thresh_sel
            case 1
                thresh = 0;
            case 2
                %PFA = 0.0025 / 100 pixels per 200x200 image
                thresh = 3.6073e-04;
            case 3
                %PFA = 5e-04 / 20 pixels per 200x200 image
                thresh = 8.7122e-04;
            case 4
                %PFA = 1.25e-04 / 5 pixels per 200x200 image
                thresh = 0.0018;
        end
    case 4
        method = 'OSP';
        thresh = 0;
    case 3
        method = 'RSDF';
        switch thresh_sel
            case 1
                thresh = 0;
            case 2
                %PFA = 0.0025 / 100 pixels per 200x200 image
                thresh = 4.3662e-04;
            case 3
                %PFA = 5e-04 / 20 pixels per 200x200 image
                thresh = 9.8935e-04;
            case 4
                %PFA = 1.25e-04 / 5 pixels per 200x200 image
                thresh = 0.0019;
        end
end

thresh
algorithm
%Check if the required wavelengths are available for RSDF
if algorithm == 1 | algorithm == 3
    req_wav  =700:10:900;
    flag_wav=1;
    for i=1:21
        if min(abs(handles.datainfo.Wavelengths - req_wav(i))) > 0
            warning('RSDF is not applicable with current wavelengths, AMF is used instead');
            method = 'Simple';
            flag_wav = 0;
            break;
        else
            [dummy,tmp] = min(abs(handles.datainfo.Wavelengths - req_wav(i)));
            keep_idx(i) = tmp;
        end
    end
    
    if flag_wav
        for i=1:21
            wavelengths(i) = handles.datainfo.Wavelengths(keep_idx(i));
            ttmp = handles.Recon(:,:,:,:,1,keep_idx(i));
            Recon_MB(:,:,:,:,1,i) = ttmp;
            spectra(:,i) = handles.spectra(:,keep_idx(i));
        end
        handles.datainfo.Wavelengths = wavelengths;
        handles.Recon = Recon_MB;
        handles.spectra = spectra;
        %         warndlg('RSDF is is defined for wavelength sampling of 700:10:900 nm. Extra wavelengths will be rejected.');
        handles.CurrWav = length(wavelengths);
        set(handles.listbox3,'string',handles.datainfo.Wavelengths);
        load('RSDF_data');
    else
        warndlg('RSDF is not available for the current wavelength sampling. 700:10:900 nm is required');
        set(handles.popupmenu1,'value',1);
        return;
    end
else
    cov_gl = [];
    covAll = [];
end
a=1
set(handles.text10,'String','Processing');
for run=1:size(handles.Recon,3)
    for slice=1:size(handles.Recon,4)
        R = squeeze(handles.Recon(:,:,run,slice,1,:));
        HM(:,:,1,:) = R;
        mixed = permute(HM,[4,3,1,2]);
        [umx, A, flag, message] = unmix(mixed, handles.datainfo.Wavelengths, handles.spectra, method, cov_gl, covAll, thresh, size(mixed,3),size(mixed,4));
        if flag == 0
            warndlg(message);
            return;
        end
        
        switch vis
            case 1
                umx = umx;
            case 2
                umx = sqrt(umx);
            case 3
                umx = sqrt(sqrt(umx));
        end
        if max(umx(:)>0)
            umx = umx./max(umx(:));
        end
        %         umx=addtext(umx,handles.n,handles.CurrRun,handles.CurrSlice,handles.datainfo.Wavelengths,handles.CurrWav,handles.image_width,handles.c,handles.recon_method,handles.MB_regu,handles.truncated,handles.noneg);
        handles.Unmixed(:,:,run,slice) = squeeze(umx);
        handles.CurrWav
        im = R(:,:,handles.CurrWav);
        im=im/max(im(:));
        im=addtext(im,handles.n,handles.CurrRun,slice,handles.datainfo.Wavelengths,handles.CurrWav,handles.image_width,handles.c,handles.recon_method,handles.MB_regu,handles.truncated,handles.noneg);
        Analysis(:,:,:,run,slice) = ImOverlayR(im,squeeze(umx),1,'test'); close;
    end
end
set(handles.text10,'String','Ready');
% [umx_rsdf_ecglrt, A] = unmix_AMF_all(mixed, handles.datainfo.Wavelengths, handles.spectra, method, cov_gl,1,[],[],[]);
handles.Analysis = Analysis;
curr_slice = get(handles.listbox_sliceSelect,'value');
axes(handles.axes1);
image(uint8(255*squeeze(handles.Analysis(:,:,:,handles.CurrRun,handles.CurrSlice)))); axis off; axis equal;colorbar ;
guidata(hObject, handles);
set(handles.pushbutton6,'enable','on');
set(handles.pushbuttonSelectROI,'enable','on');
set(handles.pushbuttonExpSpec,'enable','on');



% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

algorithm = get(handles.popupmenu1,'value')
thresh = get(handles.popupmenu4,'value')
vis = get(handles.popupmenu3,'value')
spec = get(handles.popupmenu2,'value')
algorithm
switch algorithm
    case 2
        method = 'Simple';
    case 1
        method = 'AMF_QLShr';
    case 3
        method = 'RSDF'
    case 4
        method = 'OSP'
end


folder_name = [ handles.PathNamesigMat handles.name];

mkdir(folder_name)
mkdir([folder_name '\overlays']);
mkdir([folder_name '\anatomical']);
mkdir([folder_name '\unmixed']);

set(handles.text10,'String','Processing');

for i=1:size(handles.Analysis,4)
    for j=1:size(handles.Analysis,5)
        save_file = [folder_name '\overlays' '\' handles.datainfo.Name '_Overlay_Slice' int2str(j) 'Run' int2str(i) 'Alg' method 'Thresh' num2str(thresh) 'Vis' num2str(vis) '.tif'];
        imwrite(squeeze(handles.Analysis(:,:,:,i,j)),save_file);
    end
end

for i=1:size(handles.Recon,3)
    for j=1:size(handles.Recon,4)
        save_file = [folder_name '\anatomical' '\' handles.datainfo.Name '_Anatomical_Slice' int2str(j) 'Run' int2str(i) 'Alg' method '.tif'];
        R =squeeze(handles.Recon1(:,:,i,j,1,handles.CurrWav));
        R = R-min(R(:)); R = R./max(R(:));
        imwrite(R,save_file);
    end
end

for i=1:size(handles.Unmixed,3)
    for j=1:size(handles.Unmixed,4)
        save_file = [folder_name '\unmixed' '\' handles.datainfo.Name '_Unmixed_Slice' int2str(j) 'Run' int2str(i) 'Alg' method 'Thresh' num2str(thresh) 'Vis' num2str(vis) '.tif'];
        imwrite(squeeze(handles.Unmixed(:,:,i,j)),save_file);
    end
end
if size(handles.datainfo.Wavelengths,2)>1
    xlswrite([folder_name '\unmixed\spectrum.xls'],[handles.datainfo.Wavelengths' handles.spectra(1,:)']);
else
    xlswrite([folder_name '\unmixed\spectrum.xls'],[handles.datainfo.Wavelengths handles.spectra(1,:)']);
end
set(handles.text10,'String','Ready');

% --- Executes on button press in pushbuttonSelectROI.
function pushbuttonSelectROI_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSelectROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if sum(strcmp(fieldnames(handles), 'Analysis')) == 0
    R = squeeze(handles.Recon(:,:,handles.CurrRun,handles.CurrSlice,1,handles.CurrWav));
    figure,imagesc(R(:,:)); colormap('gray'); axis off; axis equal;colorbar ;
else
    figure,image(uint8(255*squeeze(handles.Analysis(:,:,:,handles.CurrRun,handles.CurrSlice)))); axis off; axis equal;colorbar ;
end
roi = roipoly();
R = squeeze(handles.Recon(:,:,handles.CurrRun,handles.CurrSlice,1,:));
for i=1:size(R,3)
    tmp = R(:,:,i);
    spec_int(i) = mean(tmp(roi));
end

close;
axes(handles.axes3);
plot(handles.datainfo.Wavelengths, spec_int);
grid on;
axis tight;
handles.spec_int = spec_int;
handles.x_int = -1;
handles.y_int = -1;
guidata(hObject, handles);

% --- Executes on button press in pushbuttonExpSpec.
function pushbuttonExpSpec_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonExpSpec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder_name = [ handles.PathNamesigMat handles.name];
if handles.x_int>0
    file = [folder_name '\' handles.datainfo.Name '_Run' int2str(handles.CurrRun) '_Slice' int2str(handles.CurrSlice) '_x' int2str(handles.x_int) '_y' int2str(handles.y_int) '.xls'];
else
    file = [folder_name '\' handles.datainfo.Name '_Run' int2str(handles.CurrRun) '_Slice' int2str(handles.CurrSlice) '_ROI.xls'];
end
if size(handles.datainfo.Wavelengths,2)>1
    xlswrite(file,[handles.datainfo.Wavelengths'  handles.spec_int']);
else
    xlswrite(file,[handles.datainfo.Wavelengths  handles.spec_int']);
end


function getMousePositionOnImage(src, event)
handles = guidata(src);

cursorPoint = get(handles.axes1, 'CurrentPoint');
curX = cursorPoint(1,1);
curY = cursorPoint(1,2);

xLimits = get(handles.axes1, 'xlim');
yLimits = get(handles.axes1, 'ylim');

if (curX > min(xLimits) && curX < max(xLimits) && curY > min(yLimits) && curY < max(yLimits))
    if sum(strcmp(fieldnames(handles), 'Recon')) == 1
        R = squeeze(handles.Recon(:,:,handles.CurrRun,handles.CurrSlice,1,:));
        for i=1:size(R,3)
            tmp = R(round(curY),round(curX),i);
            spec_int(i) = tmp;
        end
        axes(handles.axes3);
        plot(handles.datainfo.Wavelengths, spec_int);
        grid on;
        axis tight;
        handles.spec_int = spec_int;
        handles.x_int = round(curX);
        handles.y_int = round(curY);
    end
else
    disp('Cursor is outside bounds of image.');
end
guidata(src, handles);



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5
algoRecon = get(handles.popupmenu5, 'value');
algoRecon

switch algoRecon
    case 1
        handles.recon_method = 'MB_Tik';
    case 2
        handles.recon_method = 'MB_TVL1';
    case 3
        handles.recon_method = 'WaveletPacket';
    case 4
        handles.recon_method = 'BackProjection';
        %handles.recon_method_image_select='full';
end
guidata(hObject, handles);

if strcmp(handles.recon_method,'MB_Tik')||strcmp(handles.recon_method,'MB_TVL1')
    set(handles.text16,'String',handles.MB_regu);
elseif strcmp(handles.recon_method,'WaveletPacket')
    set(handles.text16,'String',handles.truncated)
end

% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object deletion, before destroying properties.
function popupmenu1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over popupmenu1.
function popupmenu1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on key press with focus on popupmenu1 and none of its controls.
function popupmenu1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttonRecon.
function pushbuttonRecon_Callback(hObject, ~, handles)
% hObject    handle to pushbuttonRecon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


init_folder_name = pwd;
sub_folder_name = 'MB matrices';
folder_name = [init_folder_name '\' sub_folder_name];
datainfo=handles.datainfo;
handles.Recon=zeros(handles.n,handles.n);
h = waitbar(0,'1','Name',[' Time left...'],...
    'CreateCancelBtn',...
    'setappdata(gcbf,''canceling'',1)');

if strcmp(handles.recon_method,'MB_Tik')||strcmp(handles.recon_method,'MB_TVL1')||strcmp(handles.recon_method,'WaveletPacket')
    
    if exist([folder_name '\' 'A_mat_t_res_',num2str(handles.time_res),'_',num2str(handles.n),'x',num2str(handles.n),'_width_',num2str(handles.image_width*1e3),'_c_',num2str(handles.c),'.mat'])
        load ([folder_name '\' 'A_mat_t_res_',num2str(handles.time_res),'_',num2str(handles.n),'x',num2str(handles.n),'_width_',num2str(handles.image_width*1e3),'_c_',num2str(handles.c),'.mat'])
    else
        A_mat = Calculate_MatrixMB_Luis(handles.c,handles.n,handles.image_width,handles.t,handles.datainfo.HWDesc.Radius,handles.angle_sensor,handles.n_angles);
        save (([folder_name '\' 'A_mat_t_res_',num2str(handles.time_res),'_',num2str(handles.n),'x',num2str(handles.n),'_width_',num2str(handles.image_width*1e3),'_c_',num2str(handles.c),'.mat']) , 'A_mat','-v7.3')
        
    end
end

if strcmp(handles.recon_method,'WaveletPacket')
    wl_name='db6';
    depth_proj=2;
    depth_im=2;
    Sig=zeros(handles.n,handles.n);
    over_fact=ceil(1.4*size(A_mat,1)/size(A_mat,2));
    [temp,L4]=wavedec2(Sig,depth_im,wl_name);
    
    if exist(['Ainv_A_GT_',num2str(handles.WP_global_threshold),'_',num2str(handles.truncated),'_t_res_',num2str(handles.time_res),'_',num2str(handles.n),'x',num2str(handles.n),'_width_',num2str(handles.image_width*1e3),'_c_',num2str(handles.c),'.mat']) &&exist (['places_',num2str(handles.n),'x',num2str(handles.n),'.mat'])&& exist (['BB_part_temp_',num2str(handles.n),'x',num2str(handles.n),'.mat'])
        
        load (['Ainv_A_GT_',num2str(handles.WP_global_threshold),'_',num2str(handles.truncated),'_t_res_',num2str(handles.time_res),'_',num2str(handles.n),'x',num2str(handles.n),'_width_',num2str(handles.image_width*1e3),'_c_',num2str(handles.c),'.mat'])
        load (['places_',num2str(handles.n),'x',num2str(handles.n),'.mat'])
        load (['BB_part_temp_',num2str(handles.n),'x',num2str(handles.n),'.mat'])
    else
        
        max_s=0;
        
        for kk=1:4^(depth_im)
            %kk
            base=return_base_arb_lev_even_db(kk,depth_im,wl_name) ;
            BB_part=calc_B_part_arb_lev(base,handles.n,handles.n,L4,depth_im);
            places_2_check=find_important_places2(BB_part,0.2);
            BB_part_temp{kk}=sparse(BB_part(:,places_2_check));
            AA_temp2=Mat_wlet_dec_arb_lev_even_db(A_mat*BB_part_temp{kk},handles.sizeT,handles.datainfo.HWDesc.NumDetectors,depth_proj,wl_name);
            places{kk}=choose_time_proj_low_mem2(AA_temp2,over_fact,size(BB_part_temp{kk},2));
            A_part=AA_temp2(places{kk},:);
            [u,s,v]=svd(full(A_part),'econ');
            
            if handles.WP_global_threshold==1
                switch kk
                    case 1
                        max_s=max(s(:));
                end
            elseif handles.WP_global_threshold==0
                max_s=max(s(:));
            end
            
            Ainv_A{kk}=TSVD_YYH(u,s,v,handles.truncated,max_s);
            clear u s v
            
        end
        
        save (([folder_name '\' 'Ainv_A_GT_',num2str(handles.WP_global_threshold),'_',num2str(handles.truncated),'_t_res_',num2str(handles.time_res),'_',num2str(handles.n),'x',num2str(handles.n),'_width_',num2str(handles.image_width*1e3),'_c_',num2str(handles.c),'.mat']) ,'Ainv_A','-v7.3')
        save (([folder_name '\' 'places_',num2str(handles.n),'x',num2str(handles.n),'.mat']) ,'places')
        save (([folder_name '\' 'BB_part_temp_',num2str(handles.n),'x',num2str(handles.n),'.mat']) ,'BB_part_temp')
        
        clear  AA_temp2 A_part BB_part places_2_check
        
    end
end

set(handles.text10,'String','Reconstructing');

avgs = handles.datainfo.MeasurementDesc.Averages;
runs = handles.datainfo.RunNum;
repetitions = handles.datainfo.RepNum;
slices = numel(handles.datainfo.ZPositions);
lambdas = numel(handles.datainfo.Wavelengths);
total_num = runs*slices*lambdas*repetitions;
cnt = 1;

%% extract a slice to test recon with
for run_idx=1:runs
    for slc = 1:slices
        for wvl = 1:lambdas
            for rep = 1:repetitions
                tic
                
                %% select signals depending on the number of wavelengths
                if(lambdas == 1 && avgs == 1)
                    sigMat = squeeze(handles.sigMat(:, :, run_idx, slc )); % only 5 dimensions if single wavelength
                end
                if (lambdas == 1 && avgs ~= 1)
                    sigMat = squeeze(handles.sigMat(:, :, run_idx, slc, :)); % 6 dimensions if multi-wavelength
                end
                if (lambdas ~= 1)
                    sigMat = squeeze(handles.sigMat(:, :, run_idx, slc, rep, wvl)); % 6 dimensions if multi-wavelength
                end
                
                filter_f = [handles.filter_min handles.filter_max];%[0.04e6 8e6];      % a band-pass filter (originlay between 100kHz and 8MHz), or maybe try the low pass a little lower
                sigMat = filter_function(sigMat, filter_f, handles.datainfo.HWDesc.SamplingFrequency);
                
                switch handles.recon_method
                    case 'BackProjection'
                        Recon_tmp = backproject_luis(sigMat,handles.n,handles.datainfo.HWDesc.Radius,handles.angle_sensor,handles.c,'full',handles.ts,handles.datainfo.HWDesc.SamplingFrequency,handles.image_width,0,0);
                        if handles.noneg
                            Recon_tmp=max(Recon_tmp,0);
                        end
                end
                
                if strcmp(handles.recon_method,'MB_Tik')||strcmp(handles.recon_method,'MB_TVL1')||strcmp(handles.recon_method,'WaveletPacket')
                    for j=1:handles.datainfo.HWDesc.NumDetectors
                        sigMat2(:,j) = interp1(handles.ts,sigMat(:,j),handles.t);       % interpolate raw data to smaller resolution (downsampling)
                    end
                    clear sigMat;
                    b_vec = reshape(sigMat2, handles.sizeT*handles.datainfo.HWDesc.NumDetectors,1);           % reshape raw data matrix to column-major format
                    clear sigMat2
                end
                if strcmp(handles.recon_method,'MB_Tik')||strcmp(handles.recon_method,'MB_TVL1')
                    Recon_tmp  = reconstruction( A_mat, b_vec,handles.n,handles.recon_method,handles.MB_regu,handles.noneg);
                end
                if strcmp(handles.recon_method,'WaveletPacket')
                    Recon_tmp=reconstruction_WP(A_mat,Ainv_A, b_vec,places,BB_part_temp,handles.n, handles.t,handles.datainfo.HWDesc.NumDetectors,handles.noneg);
                end
                
                handles.Recon(:, :, run_idx, slc, rep, wvl) = Recon_tmp;
                Recon_tmp=addtext(Recon_tmp,handles.n,run_idx,slc,handles.datainfo.Wavelengths,wvl,handles.image_width,handles.c,handles.recon_method,handles.MB_regu,handles.truncated,handles.noneg);
                handles.Recon1(:, :, run_idx, slc, rep, wvl) = Recon_tmp;
                
                clear Recon_tmp;
                time(cnt) = toc;
                frames_left = total_num - cnt;
                time_left = mean(time).*frames_left./(60);
                disp(['Time to completion: ', num2str(time_left), ' mins...']);
                drawnow();
                if getappdata(h,'canceling')
                    delete(h)
                    return
                end
                waitbar(cnt/total_num,h,sprintf(['Time to completion: ', num2str(time_left), ' mins...']))
                cnt = cnt+1;
                
            end
        end
    end
end

delete(h)


if sum(strcmp(fieldnames(datainfo), 'MeasurementDesc'))>0
    if sum(strcmp(fieldnames(datainfo.MeasurementDesc), 'WaterAbsorptionCoeff'))>0
        
        abs_energy = handles.datainfo.MeasurementDesc.WaterAbsorptionCoeff;
        Tissue_Fl = ones(1,size(handles.Recon,6)).*exp(-abs_energy.*0.04)';
        Tissue_Fl = Tissue_Fl./norm(Tissue_Fl,2);
        for i=1:size(handles.Recon,6)
            handles.Recon(:,:,:,:,:,i) = handles.Recon(:,:,:,:,:,i)./Tissue_Fl(i);
        end
        MSOT_cal = linspace(1,0.85,length(datainfo.Wavelengths));
        for i=1:size(handles.Recon,6)
            handles.Recon(:,:,:,:,:,i) = handles.Recon(:,:,:,:,:,i)./MSOT_cal(i);
        end
    end
end

if sum(strcmp(fieldnames(datainfo), 'MeasurementDesc'))>0
    if sum(strcmp(fieldnames(datainfo.MeasurementDesc), 'WaterAbsorptionCoeff'))>0
        
        abs_energy = handles.datainfo.MeasurementDesc.WaterAbsorptionCoeff;
        Tissue_Fl = ones(1,size(handles.Recon1,6)).*exp(-abs_energy.*0.04)';
        Tissue_Fl = Tissue_Fl./norm(Tissue_Fl,2);
        for i=1:size(handles.Recon1,6)
            handles.Recon1(:,:,:,:,:,i) = handles.Recon1(:,:,:,:,:,i)./Tissue_Fl(i);
        end
        MSOT_cal = linspace(1,0.85,length(datainfo.Wavelengths));
        for i=1:size(handles.Recon1,6)
            handles.Recon1(:,:,:,:,:,i) = handles.Recon1(:,:,:,:,:,i)./MSOT_cal(i);
        end
    end
end

recon_method=handles.recon_method ;
image_width=handles.image_width;
switch handles.recon_method
    case 'BackProjection'
        
        name=['\Recon_' handles.datainfo.Name '_BackProjection_full_' num2str(handles.n) 'x' num2str(handles.n) '_width_' num2str(handles.image_width*1e3),'_c_',num2str(handles.c) '_NN_' num2str(handles.noneg) '_filter_' num2str(handles.filter_min/1e6) '_' num2str(handles.filter_max/1e6)];
        path = [handles.PathNamesigMat name];
        mkdir(path);
        save_file=[handles.PathNamesigMat name name '.mat'];
        Recon =squeeze(handles.Recon1(:,:,:,:,:,:));
        save(save_file,'Recon','datainfo','recon_method','image_width');
        
    case 'MB_Tik'
        
        name=['\Recon_' handles.datainfo.Name '_MB_Tik_regu_',num2str(handles.MB_regu),'_t_res_',num2str(handles.time_res),'_' num2str(handles.n) 'x' num2str(handles.n) '_width_' num2str(handles.image_width*1e3),'_c_',num2str(handles.c) '_NN_' num2str(handles.noneg) '_filter_' num2str(handles.filter_min/1e6) '_' num2str(handles.filter_max/1e6)];
        path = [handles.PathNamesigMat name];
        mkdir(path);
        save_file=[handles.PathNamesigMat name name '.mat'];
        Recon =squeeze(handles.Recon1(:,:,:,:,:,:));
        save(save_file,'Recon','datainfo','recon_method','image_width');
    case 'MB_TVL1'
        
        name=['\Recon_' handles.datainfo.Name '_MB_TVL1_regu_',num2str(handles.MB_regu),'_t_res_',num2str(handles.time_res),'_' num2str(handles.n) 'x' num2str(handles.n) '_width_' num2str(handles.image_width*1e3),'_c_',num2str(handles.c) '_NN_' num2str(handles.noneg) '_filter_' num2str(handles.filter_min/1e6) '_' num2str(handles.filter_max/1e6)];
        path = [handles.PathNamesigMat name];
        mkdir(path);
        save_file=[handles.PathNamesigMat name name '.mat'];
        Recon =squeeze(handles.Recon1(:,:,:,:,:,:));
        save(save_file,'Recon','datainfo','recon_method','image_width');
    case 'WaveletPacket'
        
        name=['\Recon_' handles.datainfo.Name '_WaveletPacket_',num2str(handles.truncated),'_t_res_',num2str(handles.time_res),'_' num2str(handles.n) 'x' num2str(handles.n) '_width_' num2str(handles.image_width*1e3),'_c_',num2str(handles.c) '_NN_' num2str(handles.noneg) '_filter_' num2str(handles.filter_min/1e6) '_' num2str(handles.filter_max/1e6)];
        path = [handles.PathNamesigMat name];
        mkdir(path);
        save_file=[handles.PathNamesigMat name name '.mat'];
        Recon =squeeze(handles.Recon1(:,:,:,:,:,:));
        save(save_file,'Recon','datainfo','recon_method','image_width');
end

handles.name=name;

set(handles.text10,'String','Recon saved');

R = squeeze(handles.Recon(:,:,handles.CurrRun,handles.CurrSlice,1,handles.CurrWav));
R=addtext(R,handles.n,handles.CurrRun,handles.CurrSlice,handles.datainfo.Wavelengths,handles.CurrWav,handles.image_width,handles.c,handles.recon_method,handles.MB_regu,handles.truncated,handles.noneg);
axes(handles.axes1);
imagesc(R(:,:)); colormap('gray'); axis off; axis equal;colorbar ;

guidata(hObject, handles);


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttonReconCurrSlice.
function pushbuttonReconCurrSlice_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonReconCurrSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
init_folder_name = pwd;
sub_folder_name = 'MB matrices';
folder_name = [init_folder_name '\' sub_folder_name];

if strcmp(handles.recon_method,'MB_Tik')||strcmp(handles.recon_method,'MB_TVL1')||strcmp(handles.recon_method,'WaveletPacket')
    
    % check for model matrix if saved, otherwise build it
    if exist([folder_name '\' 'A_mat_t_res_',num2str(handles.time_res),'_',num2str(handles.n),'x',num2str(handles.n),'_width_',num2str(handles.image_width*1e3),'_c_',num2str(handles.c),'.mat'], 'file')
        load([folder_name '\' 'A_mat_t_res_',num2str(handles.time_res),'_',num2str(handles.n),'x',num2str(handles.n),'_width_',num2str(handles.image_width*1e3),'_c_',num2str(handles.c),'.mat'])
    else
        A_mat = Calculate_MatrixMB_Luis(handles.c,handles.n,handles.image_width,handles.t,handles.datainfo.HWDesc.Radius,handles.angle_sensor,handles.n_angles);
        if ~exist(folder_name, 'dir')
            mkdir(folder_name);
        end
        save(([folder_name '\' 'A_mat_t_res_',num2str(handles.time_res),'_',num2str(handles.n),'x',num2str(handles.n),'_width_',num2str(handles.image_width*1e3),'_c_',num2str(handles.c),'.mat']) , 'A_mat','-v7.3')
    end
end

if strcmp(handles.recon_method,'WaveletPacket')
    wl_name='db6';
    depth_proj=2;
    depth_im=2;
    Sig=zeros(handles.n,handles.n);
    over_fact=ceil(1.4*size(A_mat,1)/size(A_mat,2));
    [temp,L4]=wavedec2(Sig,depth_im,wl_name);
    
    if exist([folder_name '\' 'Ainv_A_GT_',num2str(handles.WP_global_threshold),'_',num2str(handles.truncated),'_t_res_',num2str(handles.time_res),'_',num2str(handles.n),'x',num2str(handles.n),'_width_',num2str(handles.image_width*1e3),'_c_',num2str(handles.c),'.mat']) &&exist (['places_',num2str(handles.n),'x',num2str(handles.n),'.mat'])&& exist (['BB_part_temp_',num2str(handles.n),'x',num2str(handles.n),'.mat'])
        
        load ([folder_name '\' 'Ainv_A_GT_',num2str(handles.WP_global_threshold),'_',num2str(handles.truncated),'_t_res_',num2str(handles.time_res),'_',num2str(handles.n),'x',num2str(handles.n),'_width_',num2str(handles.image_width*1e3),'_c_',num2str(handles.c),'.mat'])
        load ([folder_name '\' 'places_',num2str(handles.n),'x',num2str(handles.n),'.mat'])
        load ([folder_name '\' 'BB_part_temp_',num2str(handles.n),'x',num2str(handles.n),'.mat'])
    else
        
        max_s=0;
        
        for kk=1:4^(depth_im)
            kk
            base=return_base_arb_lev_even_db(kk,depth_im,wl_name) ;
            BB_part=calc_B_part_arb_lev(base,handles.n,handles.n,L4,depth_im);
            places_2_check=find_important_places2(BB_part,0.2);
            BB_part_temp{kk}=sparse(BB_part(:,places_2_check));
            AA_temp2=Mat_wlet_dec_arb_lev_even_db(A_mat*BB_part_temp{kk},handles.sizeT,handles.datainfo.HWDesc.NumDetectors,depth_proj,wl_name);
            places{kk}=choose_time_proj_low_mem2(AA_temp2,over_fact,size(BB_part_temp{kk},2));
            A_part=AA_temp2(places{kk},:);
            [u,s,v]=svd(full(A_part),'econ');
            
            if handles.WP_global_threshold==1
                switch kk
                    case 1
                        max_s=max(s(:));
                end
            elseif handles.WP_global_threshold==0
                max_s=max(s(:));
            end
            
            Ainv_A{kk}=TSVD_YYH(u,s,v,handles.truncated,max_s);
            clear u s v
            
        end
        
        save (([folder_name '\' 'Ainv_A_GT_',num2str(handles.WP_global_threshold),'_',num2str(handles.truncated),'_t_res_',num2str(handles.time_res),'_',num2str(handles.n),'x',num2str(handles.n),'_width_',num2str(handles.image_width*1e3),'_c_',num2str(handles.c),'.mat']) ,'Ainv_A','-v7.3')
        save (([folder_name '\' 'places_',num2str(handles.n),'x',num2str(handles.n),'.mat']) ,'places')
        save (([folder_name '\' 'BB_part_temp_',num2str(handles.n),'x',num2str(handles.n),'.mat']) ,'BB_part_temp')
        
        clear  AA_temp2 A_part BB_part places_2_check
        
    end
end

avgs = handles.datainfo.MeasurementDesc.Averages;
runs = handles.datainfo.RunNum;
repetitions = handles.datainfo.RepNum;
slices = numel(handles.datainfo.ZPositions);
lambdas = numel(handles.datainfo.Wavelengths);
total_num = runs*slices*lambdas*repetitions;
cnt = 1;

%% extract a slice to test recon with
for run_idx=handles.CurrRun:handles.CurrRun
    for slc = handles.CurrSlice:handles.CurrSlice
        for wvl = handles.CurrWav:handles.CurrWav
            for rep = 1:repetitions
                %% select signals depending on the number of wavelengths
                if(lambdas == 1 && avgs == 1)
                    sigMat = squeeze(handles.sigMat(:, :, run_idx, slc )); % only 5 dimensions if single wavelength
                end
                if (lambdas == 1 && avgs ~= 1)
                    sigMat = squeeze(handles.sigMat(:, :, run_idx, slc, :)); % 6 dimensions if multi-wavelength
                end
                if (lambdas ~= 1)
                    sigMat = squeeze(handles.sigMat(:, :, run_idx, slc, rep, wvl)); % 6 dimensions if multi-wavelength
                end
                
                filter_f = [handles.filter_min handles.filter_max];%[0.04e6 8e6];      % a band-pass filter (originlay between 100kHz and 8MHz), or maybe try the low pass a little lower
                sigMat = filter_function(sigMat, filter_f, handles.datainfo.HWDesc.SamplingFrequency);
                tic
                switch handles.recon_method
                    case 'BackProjection'
                        Recon_tmp = backproject_luis(sigMat,handles.n,handles.datainfo.HWDesc.Radius,handles.angle_sensor,handles.c,'full',handles.ts,handles.datainfo.HWDesc.SamplingFrequency,handles.image_width,0,0);
                        if handles.noneg
                            Recon_tmp=max(Recon_tmp,0);
                        end
                end
                
                if strcmp(handles.recon_method,'MB_Tik')||strcmp(handles.recon_method,'MB_TVL1')||strcmp(handles.recon_method,'WaveletPacket')
                    for j=1:handles.datainfo.HWDesc.NumDetectors
                        sigMat2(:,j) = interp1(handles.ts,sigMat(:,j),handles.t);       % interpolate raw data to smaller resolution (downsampling)
                    end
                    clear sigMat;
                    b_vec = reshape(sigMat2, handles.sizeT*handles.datainfo.HWDesc.NumDetectors,1);           % reshape raw data matrix to column-major format
                    clear sigMat2
                end
                if strcmp(handles.recon_method,'MB_Tik')||strcmp(handles.recon_method,'MB_TVL1')
                    Recon_tmp  = reconstruction( A_mat, b_vec,handles.n,handles.recon_method,handles.MB_regu,handles.noneg);
                end
                if strcmp(handles.recon_method,'WaveletPacket')
                    Recon_tmp  =reconstruction_WP(A_mat,Ainv_A, b_vec,places,BB_part_temp,handles.n, handles.t,handles.datainfo.HWDesc.NumDetectors,handles.noneg);
                end
                
                Recon_tmp=addtext(Recon_tmp,handles.n,run_idx,slc,handles.datainfo.Wavelengths,wvl,handles.image_width,handles.c,handles.recon_method,handles.MB_regu,handles.truncated,handles.noneg);
                
                % clear handles.Recon handles.axes1
                handles.Recon=zeros(handles.n,handles.n);
                %                  handles.n
                %                  size(handles.Recon)
                %                  size(Recon_tmp)
                handles.Recon(:, :, run_idx, slc, rep, wvl) = Recon_tmp;
                clear Recon_tmp;
                toc
            end
        end
    end
end

R = squeeze(handles.Recon(:,:,handles.CurrRun,handles.CurrSlice,1,handles.CurrWav));
R=addtext(R,handles.n,handles.CurrRun,handles.CurrSlice,handles.datainfo.Wavelengths,handles.CurrWav,handles.image_width,handles.c,handles.recon_method,handles.MB_regu,handles.truncated,handles.noneg);
axes(handles.axes1);
imagesc(R(:,:)); colormap('gray'); axis off; axis equal;colorbar ;
% figure, imagesc(R(:,:)); colormap('gray'); axis off; axis equal;colorbar ;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pushbutton_LoadsigMat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_LoadsigMat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in popupmenu7.
function popupmenu7_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu7 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu7
image_width = get(handles.popupmenu7, 'value');
image_width

switch image_width
    case 1
        handles.image_width = 20e-3;
        handles.n=200;
        set(handles.edit6,'enable','off');
        set(handles.edit7,'enable','off');
    case 2
        handles.image_width = 20e-3;
        handles.n=267;
        set(handles.edit6,'enable','off');
        set(handles.edit7,'enable','off');
    case 3
        handles.image_width = 25e-3;
        handles.n=250;
        set(handles.edit6,'enable','off');
        set(handles.edit7,'enable','off');
    case 4
        handles.image_width = 25e-3;
        handles.n=334;
        set(handles.edit6,'enable','off');
        set(handles.edit7,'enable','off');
    case 5
        set(handles.edit6,'enable','on');
        set(handles.edit7,'enable','on');
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes during object creation, after setting all properties.
function text9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
handles.image_width=1e-3*str2double(get(handles.edit5,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
resolution=1e-6*str2double(get(handles.edit7,'String'))
handles.n=round(handles.image_width/resolution);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
noneg=get(hObject,'Value');
noneg
handles.noneg=noneg;
guidata(hObject, handles);

function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
handles.filter_min=1e6*str2double(get(handles.edit7,'String'));

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
handles.filter_max=1e6*str2double(get(handles.edit8,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider7_Callback(hObject, ~, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.c=get(handles.slider7,'Value');
handles.c=round(handles.c);
set(handles.slider7,'Value',handles.c);
guidata(hObject, handles);
set(handles.text18,'String',handles.c);


% --- Executes during object creation, after setting all properties.
function slider7_CreateFcn(hObject, ~, ~)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider8_Callback(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.regu=get(handles.slider8,'Value');
handles.regu=round(handles.regu);
set(handles.slider8,'Value',handles.regu);

handles.truncated = 0.1*(1.1^(handles.regu-20));
handles.MB_regu=1e6*(2^((handles.regu-20)));

if strcmp(handles.recon_method,'MB_Tik')||strcmp(handles.recon_method,'MB_TVL1')
    set(handles.text16,'String',handles.MB_regu);
elseif strcmp(handles.recon_method,'WaveletPacket')
    set(handles.text16,'String',handles.truncated)
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function slider8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function text18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function uipanel3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
