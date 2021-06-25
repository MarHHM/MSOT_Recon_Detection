function varargout = testGui(varargin)
% TESTGUI MATLAB code for testGui.fig
%      TESTGUI, by itself, creates a new TESTGUI or raises the existing
%      singleton*.
%
%      H = TESTGUI returns the handle to a new TESTGUI or the handle to
%      the existing singleton*.
%
%      TESTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TESTGUI.M with the given input arguments.
%
%      TESTGUI('Property','Value',...) creates a new TESTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before testGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to testGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help testGui

% Last Modified by GUIDE v2.5 01-Jul-2016 16:18:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @testGui_OpeningFcn, ...
                   'gui_OutputFcn',  @testGui_OutputFcn, ...
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


% --- Executes just before testGui is made visible.
function testGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to testGui (see VARARGIN)

% re-name
clc;
javaaddpath MSOTBeans\msotbeans.jar
javaaddpath java_class\recon.jar                   
javaaddpath java_class\patbeans.jar                
javaaddpath java_class\xmlbeans-2.5.0\lib\xbean.jar
    
set(handles.figure1, 'Name', 'MSOT Molecular Target Detection');

% Choose default command line output for testGui
handles.output = hObject;
addpath('unmixing code');
axis off; axis equal;

handles = guidata(hObject);
handles.output = 1;
guidata(hObject, handles);

handles = guidata(hObject);
guidata(hObject, handles);

set(handles.pushbuttonSpec,'enable','off');
set(handles.pushbuttonAnalyze,'enable','off');
set(handles.pushbutton6,'enable','off');
set(handles.pushbuttonCalSpec,'enable','off');

set(gcf, 'WindowButtonDownFcn', @getMousePositionOnImage);
pan off % Panning will interfere with this code

guidata(hObject, handles);

% UIWAIT makes testGui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = testGui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
slice = get(handles.listbox1,'value');
handles.CurrSlice = slice;
guidata(hObject, handles);
if sum(strcmp(fieldnames(handles), 'Analysis')) == 0
%     R = squeeze(handles.Recon(:,:,handles.CurrRep,handles.CurrSlice,1,handles.CurrWav));
    R = squeeze(handles.Recon(:,:,handles.CurrSlice,handles.CurrWav));          % edited to work with my data format
    axes(handles.axes1);
    imagesc(R(:,:)); colormap('gray'); axis off; axis equal;
else
    axes(handles.axes1);
    image(uint8(255*squeeze(handles.Analysis(:,:,:,handles.CurrRep,handles.CurrSlice)))); axis off; axis equal;
end

% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
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
handles.CurrRep = run;
guidata(hObject, handles);
if sum(strcmp(fieldnames(handles), 'Analysis')) == 0
    R = squeeze(handles.Recon(:,:,handles.CurrRep,handles.CurrSlice,1,handles.CurrWav));
    axes(handles.axes1);
    imagesc(R(:,:)); colormap('gray'); axis off; axis equal;
else
    axes(handles.axes1);
    image(squeeze(uint8(255*handles.Analysis(:,:,:,handles.CurrRep,handles.CurrSlice)))); axis off; axis equal;
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
%     R = squeeze(handles.Recon(:,:,handles.CurrRep,handles.CurrSlice,1,handles.CurrWav));
    R = squeeze(handles.Recon(:,:,handles.CurrSlice,handles.CurrWav));
    axes(handles.axes1);
    imagesc(R(:,:)); colormap('gray'); axis off; axis equal;
else
    axes(handles.axes1);
    image(uint8(255*squeeze(handles.Analysis(:,:,:,handles.CurrRep,handles.CurrSlice)))); axis off; axis equal;
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
            handles.spectra = LoadSpectra(folder_name,handles.datainfo.Wavelengths)
            guidata(hObject, handles);
         %Alexa 750
        case 3
            folder_name = 'spectra\SpectralSpecifications_AF750';
            handles.spectra = LoadSpectra(folder_name,handles.datainfo.Wavelengths)
            guidata(hObject, handles);
         %Alexa 790
        case 4
            folder_name = 'spectra\SpectralSpecifications_AF790';
            handles.spectra = LoadSpectra(folder_name,handles.datainfo.Wavelengths)
            guidata(hObject, handles);
         %IntSense
        case 5
            folder_name = 'spectra\SpectralSpecifications_IntSense';
            handles.spectra = LoadSpectra(folder_name,handles.datainfo.Wavelengths)
            guidata(hObject, handles);
        case 6
            folder_name = 'spectra\Spectral_Specifications_DiR';
            handles.spectra = LoadSpectra(folder_name,handles.datainfo.Wavelengths)
            guidata(hObject, handles);
        case 7
            folder_name = 'spectra\Spectral specifications Pigment';
            handles.spectra = LoadSpectra(folder_name,handles.datainfo.Wavelengths)
            guidata(hObject, handles);
        case 8
            folder_name = 'spectra\SpectralSpecifications_iRFP';
            handles.spectra = LoadSpectra(folder_name,handles.datainfo.Wavelengths)
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
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on button press in pushbutton_LoadRecon.
function pushbutton_LoadRecon_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_LoadRecon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileNameRecon,PathNameRecon,FilterIndexRecon] = uigetfile({'*.mat';'*.msot'});     %.mat->load recon || .msot->load raw acquisition data

if FileNameRecon

    javaaddpath MSOTBeans\msotbeans.jar
    javaaddpath java_class\recon.jar                   
    javaaddpath java_class\patbeans.jar                
    javaaddpath java_class\xmlbeans-2.5.0\lib\xbean.jar

    if FileNameRecon(end-3:end) == 'msot'           % loading raw acquisition data
        
        [datainfo]  = loadMSOT( [PathNameRecon '\' FileNameRecon] );
        
        c=0;
        for ri = 1:numel(datainfo.ReconNode)
            c = c+1;
            reconlist{c} = [datainfo.Name ' -- ' datainfo.ReconNode(ri).Name];
            handles.reconref{c} = [num2str(1) '-' num2str(ri)];
        end
        handles.reconNodes = c;
        set(handles.lb_recon,'String',reconlist');
        set(handles.lb_recon,'Value',1);
        
        %to be copied in listbox
        selInd = get(handles.lb_recon,'Value');
        contents = get(handles.lb_recon,'String');
        FileNameRecon = contents{selInd};   

        [Recon_MB_L1 wls zpos ts datainfo] = loadMSOTRecon(datainfo,selInd);

        % Need to make sure on how the data are loaded
        if length(size(Recon_MB_L1)) == 6
            Nx = size(Recon_MB_L1,1);
            Ny = size(Recon_MB_L1,2);
            Nr = size(Recon_MB_L1,3);
            Ns = size(Recon_MB_L1,4);
            Nt = size(Recon_MB_L1,5);
            Nw = size(Recon_MB_L1,6);
        elseif length(size(Recon_MB_L1)) == 7
            Recon_MB_L1 = mean(Recon_MB_L1,7);
            Nx = size(Recon_MB_L1,1);
            Ny = size(Recon_MB_L1,2);
            Nr = size(Recon_MB_L1,3);
            Ns = size(Recon_MB_L1,4);
            Nt = size(Recon_MB_L1,5);
            Nw = size(Recon_MB_L1,6);
            Na = size(Recon_MB_L1,7);
            Recon_MB_L1= mean(Recon_MB_L1,7);
        else
            warning('inconsistent data');
        end
        handles.Recon = Recon_MB_L1;

        handles.datainfo = datainfo;
        handles.datainfoKeep = handles.datainfo;
        set(handles.pushbuttonSpec,'enable','on');
        set(handles.pushbuttonCalSpec,'enable','on');

        % set(handles.listbox1)
        set(handles.listbox1,'string',1:length(datainfo.ZPositions));
        set(handles.listbox2,'string',1:length(datainfo.RunNum));
        set(handles.listbox3,'string',handles.datainfo.Wavelengths);

        handles.CurrSlice = 1;
        handles.CurrRep = 1;
        handles.CurrWav = length(handles.datainfo.Wavelengths);
    else            % load .mat (reconstructed data)
        % Find any variable in the file that starts with Recon
        % Depending on the dimensionality annotate
        handles.reconNodes = 1;
        disp("Loading Recon structure..");
        handles.ReconLoad = load([PathNameRecon '\' FileNameRecon]);
        load_vars = fieldnames(handles.ReconLoad);
        %iterate through the fileds until you find a 'Recon' field
        for i=1:length(load_vars)
            tmp_var = load_vars{i};
            if length(tmp_var)>4
                %if a field is found update the field index.
                if strcmp(tmp_var(1:5),'Recon')
                    idx_field = i;
                end
            end
        end
        %get the reconstructed images
        tmpRecon = getfield(handles.ReconLoad, load_vars{idx_field});

        switch length(size(tmpRecon))
            case 3
                % x-y-wavelengths
                tmp(:,:,1,1,1,:) = tmpRecon(:,:,:);
                handles.Recon = tmp;
            case 4
                % x-y-slices-wavelengths
                tmp(:,:,1,:,1,:) = tmpRecon(:,:,:,:);
                handles.Recon = tmp;
            case 5
                % x-y-runs-slices-wavelengths
                tmp(:,:,:,:,1,:) = tmpRecon(:,:,:,:,:);
                handles.Recon = tmp;
            case 6
                % x-y-runs-slices-reps-wavelengths
                tmp(:,:,:,:,1,:) = mean(tmpRecon(:,:,:,:,:,:),5);       % average across reps
                handles.Recon = tmp;
        end

        if any(strcmp('datainfo',fieldnames(handles.ReconLoad)))
            handles.datainfo = handles.ReconLoad.datainfo;
            set(handles.pushbuttonSpec,'enable','on');
            set(handles.pushbuttonCalSpec,'enable','on');

            % set(handles.listbox1)
            clear handles.listbox1;
            clear handles.listbox2;
            clear handles.listbox3;
            set(handles.listbox1,'string',1:length(handles.datainfo.ZPositions));
            set(handles.listbox2,'string',1:length(handles.datainfo.RunNum));
            set(handles.listbox3,'string',handles.datainfo.Wavelengths);
        else
            warndlg('datainfo is not provided, please provide datainfo','Warning!');
        end

        handles.CurrSlice = 1;
        handles.CurrRep = 1;
%         handles.CurrWav = length(handles.datainfo.Wavelengths);
        handles.CurrWav = 1;

    end
    
    set(handles.pushbuttonAnalyze,'enable','off');
%     clear handles.Analysis 
%     clear handles.Unmixed 
%     clear handles.spectra
    if any(strcmp('Analysis',fieldnames(handles)))
        handles = rmfield(handles,'Analysis');
%     delete(handles.Analysis);
        disp('remove handles')
    end
    
    if any(strcmp('Unmixed',fieldnames(handles)))
        handles = rmfield(handles,'Unmixed');
        disp('remove handles')
    end
    
    if any(strcmp('spectra',fieldnames(handles)))
        handles = rmfield(handles,'spectra');
        disp('remove handles')
    end

    set(handles.text9,'String',handles.datainfo.Name);
    
    R = squeeze(handles.Recon(:,:,handles.CurrRep,handles.CurrSlice,1,end));
    axes(handles.axes1);
    imagesc(R(:,:)); colormap('gray'); axis off; axis equal;
    guidata(hObject, handles);

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
function pushbuttonAnalyze_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAnalyze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
algorithm = get(handles.popupmenu1,'value');
thresh_sel = get(handles.popupmenu4,'value');
vis = get(handles.popupmenu3,'value');
% spec = get(handles.popupmenu2,'value');

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
        method = 'QL_shrinkage_adaptive_interp';
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
        method = 'RSDF_interp';
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
% thresh
% algorithm
%Check if the required wavelengths are available for RSDF
if algorithm == 1 || algorithm == 3
    req_wav  =700:1:900;
%     flag_wav=1;
    w_subset = intersect(handles.datainfo.Wavelengths,req_wav);

    if length(w_subset)<length(handles.datainfo.Wavelengths)        
        uiwait( warndlg('This algorithm was developed for wavelengths within [700-900] nm and multiple of 10. Some wavelengths of the current dataset will be discarded','modal'));
    end
    
    for i=1:length(w_subset) 
       [~,tmp] = min(abs(handles.datainfo.Wavelengths - w_subset(i)));
       keep_idx(i) = tmp;
    end

    for i=1:length(w_subset) 
        ttmp = handles.Recon(:,:,:,:,1,keep_idx(i));
        Recon_MB(:,:,:,:,1,i) = ttmp;
        spectra(:,i) = handles.spectra(:,keep_idx(i));
    end
    handles.datainfo.Wavelengths = w_subset;
    handles.Recon = Recon_MB;
    handles.spectra = spectra;
%         warndlg('RSDF is is defined for wavelength sampling of 700:10:900 nm. Extra wavelengths will be rejected.');
    if handles.CurrWav > length(w_subset) 
        handles.CurrWav = length(w_subset);
    end
    set(handles.listbox3,'string',handles.datainfo.Wavelengths);
    load('RSDF_data');  


else
    cov_gl = [];
    covAll = [];
end
set(handles.text10,'String','Processing...');
for run=1:size(handles.Recon,3)
    for slice=1:size(handles.Recon,4)
        R = squeeze(handles.Recon(:,:,run,slice,1,:));
        HM(:,:,1,:) = R;
        mixed = permute(HM,[4,3,1,2]);
        [umx, ~, flag, message] = unmix(mixed, handles.datainfo.Wavelengths, handles.spectra, method, cov_gl, covAll, thresh, size(mixed,3),size(mixed,4));
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
        handles.Unmixed(:,:,run,slice) = squeeze(umx);
        im = R(:,:,handles.CurrWav);
        Analysis(:,:,:,run,slice) = ImOverlayR(im,squeeze(umx),1,'test'); close;
    end
end
set(handles.text10,'String','Ready');
    % [umx_rsdf_ecglrt, A] = unmix_AMF_all(mixed, handles.datainfo.Wavelengths, handles.spectra, method, cov_gl,1,[],[],[]);
handles.Analysis = Analysis;
curr_slice = get(handles.listbox1,'value');
axes(handles.axes1);
image(uint8(255*squeeze(handles.Analysis(:,:,:,handles.CurrRep,handles.CurrSlice)))); axis off; axis equal;
guidata(hObject, handles);
set(handles.pushbutton6,'enable','on');


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

init_folder_name = uigetdir();

sub_folder_name = handles.datainfo.Name;
folder_name = [init_folder_name '\' sub_folder_name];

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
         R =squeeze(handles.Recon(:,:,i,j,1,handles.CurrWav));
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

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Message = {'This software implements spectral detection algorithms for MSOT molecular target detection.' ,
           'For a thorough description of the algorithms refer to ',
           '',
           'S. Tzoumas, et. al., "Unmixing molecular agents from absorbing tissue in multispectral optoacoustic tomography, Medical Imaging, IEEE Transactions on, 2014. ' ,
           '',
           'S. Tzoumas, et. al., "Statistical molecular target detection framework for multispectral optoacoustic tomography," Medical Imaging, IEEE Transactions on, 2016.' ,
           '',
           'For information on how to use the GUI please read the readme.doc file' ,
           'For reporting bugs or for any other issues please contact strtzoumas@gmail.com'}
       
Title = 'Info';
h = msgbox(Message,Title);

% --- Executes on button press in pushbuttonSelectROI.
function pushbuttonSelectROI_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSelectROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if sum(strcmp(fieldnames(handles), 'Analysis')) == 0
    R = squeeze(handles.Recon(:,:,handles.CurrRep,handles.CurrSlice,1,handles.CurrWav));
    figure,imagesc(R(:,:)); colormap('gray'); axis off; axis equal;
else
    figure,image(uint8(255*squeeze(handles.Analysis(:,:,:,handles.CurrRep,handles.CurrSlice)))); axis off; axis equal;
end
roi = roipoly();
R = squeeze(handles.Recon(:,:,handles.CurrRep,handles.CurrSlice,1,:));
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
folder_name = uigetdir();
if handles.x_int>0
    file = [folder_name '\' handles.datainfo.Name '_Run' int2str(handles.CurrRep) '_Slice' int2str(handles.CurrSlice) '_x' int2str(handles.x_int) '_y' int2str(handles.y_int) '.xls'];
else
    file = [folder_name '\' handles.datainfo.Name '_Run' int2str(handles.CurrRep) '_Slice' int2str(handles.CurrSlice) '_ROI.xls'];
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
        R = squeeze(handles.Recon(:,:,handles.CurrRep,handles.CurrSlice,1,:));
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


%Version 2

% --- Executes on selection change in lb_recon.
function lb_recon_Callback(hObject, eventdata, handles)
% hObject    handle to lb_recon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lb_recon contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lb_recon

if handles.reconNodes >1;

   if any(strcmp('datainfoKeep',fieldnames(handles)))
        handles.datainfo = handles.datainfoKeep;
   end
        
        
    recon_node = get(handles.lb_recon,'value');

    %to be copied in listbox
    selInd = get(handles.lb_recon,'Value');
    contents = get(handles.lb_recon,'String');
    FileNameRecon = contents{selInd};   

    [Recon_MB_L1 wls zpos ts datainfo] = loadMSOTRecon(handles.datainfo,selInd);

    % Need to make sure on how the data are loaded
    if length(size(Recon_MB_L1)) == 6
        Nx = size(Recon_MB_L1,1);
        Ny = size(Recon_MB_L1,2);
        Nr = size(Recon_MB_L1,3);
        Ns = size(Recon_MB_L1,4);
        Nt = size(Recon_MB_L1,5);
        Nw = size(Recon_MB_L1,6);
    elseif length(size(Recon_MB_L1)) == 7
        Recon_MB_L1 = mean(Recon_MB_L1,7);
        Nx = size(Recon_MB_L1,1);
        Ny = size(Recon_MB_L1,2);
        Nr = size(Recon_MB_L1,3);
        Ns = size(Recon_MB_L1,4);
        Nt = size(Recon_MB_L1,5);
        Nw = size(Recon_MB_L1,6);
        Na = size(Recon_MB_L1,7);
        Recon_MB_L1= mean(Recon_MB_L1,7);
    else
        warning('inconsistent data');
    end
    handles.Recon = Recon_MB_L1;

    set(handles.pushbuttonSpec,'enable','on');
    set(handles.pushbuttonCalSpec,'enable','on');
    handles.datainfo.Wavelengths

    % set(handles.listbox1)
    clear handles.listbox1;
    clear handles.listbox2;
    clear handles.listbox3;
    set(handles.listbox1,'string',1:length(datainfo.ZPositions));
    set(handles.listbox2,'string',1:length(datainfo.RunNum));
    set(handles.listbox3,'string',handles.datainfo.Wavelengths);

    handles.CurrSlice = 1;
    handles.CurrRep = 1;
    handles.CurrWav = length(handles.datainfo.Wavelengths);

     set(handles.pushbuttonAnalyze,'enable','off');
    %     clear handles.Analysis 
    %     clear handles.Unmixed 
    %     clear handles.spectra
    if any(strcmp('Analysis',fieldnames(handles)))
    handles = rmfield(handles,'Analysis');
    %     delete(handles.Analysis);
    disp('remove handles')
    end

    if any(strcmp('Unmixed',fieldnames(handles)))
    handles = rmfield(handles,'Unmixed');
    disp('remove handles')
    end

    if any(strcmp('spectra',fieldnames(handles)))
    handles = rmfield(handles,'spectra');
    disp('remove handles')
    end

    set(handles.text9,'String',handles.datainfo.Name);

    R = squeeze(handles.Recon(:,:,handles.CurrRep,handles.CurrSlice,1,end));
    axes(handles.axes1);
    imagesc(R(:,:)); colormap('gray'); axis off; axis equal;
    guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function lb_recon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lb_recon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% loading Study
