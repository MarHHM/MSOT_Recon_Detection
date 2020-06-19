classdef msotProbe < handle
    % msotProbe Collector Class of all Probe Properties
    %   Can be used to generate sensor positions file and also to plot the
    %   probe geometry
    
    properties
      % meta properties
      ProbeID = 'PR010';
      DesignID = 'cdc12358-1';
      
      % geometry
      Type@char = 'Parametric';                % Mode of Definition, 'Parametric' (requires Radius, N, Coverage) or 'Coordinates' (requiers Sensors, N)
      Radius@double = 40.5e-3;                 % radius of detector (only Parametric)
      Curvature@double = 40e-3;                % radius of curvature (in elevation)
      Coverage@double = 125;                   % coverage in degree (only Parametric)
      N@uint16 = uint16(256);                  % number of elements
      Spacing@double = 0.1e-3;                 % inter element spacing
      Height@double = 13e-3;                   % height of elements
      
      % casing related properties
      AxialOffset = 0.07;               % Distance of Membrane from COR
      LightSpotSize = 0.02;             % Width of LightSpot (FWHM, in plane)
      LightSpotHeight = 0.006;           % Height of LightSpot (FWHM, out of plane)
      
      % Couplant Properties
      CouplantSOS = 0;                  % Couplant Speed Of Sound
      CouplantCOR = 0;                  % Couplant Center of Rotation
      CouplantCurvature = 0;            % Couplant Radius of Curvature
        
      % acoustic properties
      CenterFrequency@double = 3e6;            % Center Frequency in Hz (default: 5MHz)
      Bandwidth@double = 0.6;           % Relative T/R Bandwidth (default: 0.6)
      PhaseShift@double = 0.0;         % Phase Shift (default: 0)
      Slope@double = 0;                 % Phase Slope (default: 0)
    end
    
    properties (Dependent = true)
      NumDetectors@uint16;              % Number of Detectors, read only, derived from N
    end
    
    methods
       % *****************************************************************
       %% CONSTRUCTOR
        function obj = msotProbe(varargin)
            % Create MSOT Probe Object
            %
            % msotProbe() - create with default properties
%             obj@handle;         % superclass constructor
            if nargin == 0, return; end
        end
        
        
        
        
        
        
        
        % *****************************************************************
        %% ACCESSORS
        function obj = set.Type(obj,Type)
            if strcmp(Type,'Parametric') || strcmp(Type,'Coordinates')
                obj.Type = Type;
                return
            end
            error('Invalid Type - valid types are: Parametric, Coordinates');
        end
        
        function N = get.NumDetectors(obj)
            N = obj.N;
        end
        
        function obj = set.N(obj)
            obj.N = N;
            obj.Sensors = [];
        end
        
        
        function Angles = get.Angles(obj)
            if strcmp(obj.Type,'Coordinates')
                error('Cannot output Angles for Coordinate Definition');
                return;
            end
            
        end
        
        function obj = set.Angles(obj,Angles)
            if strcmp(obj.Type,'Coordinates')
                error('Cannot set Angles for Coordinate Definition');
                return;
            end
            
            if diff(abs(Angles(obj.N/2+uint16([0 1]))+pi/2)) > eps
                fprintf('WARNING: Angles are not symmetric\n');
            end
            obj.N = uint16(numel(Angles));
            obj.Angles = reshape(Angles,1,obj.N);
            obj.Coverage = Angles(end)-Angles(1) / pi * 180;
            calculateSensors(obj);
        end
            
        function Sensors = get.Sensors(obj)
            if isempty(obj.Sensors)
                obj = calculateSensors(obj);
            end
            Sensors = obj.Sensors;
        end

        function obj = set.Sensors(obj,Sensors)
            if strcmp(obj.Type,'Parametric')
                error('Cannot set Coordinates for Parametric Definition');
                return;
            end
            if isempty(Sensors), obj.Sensors = []; return; end;
            if size(Sensors,2) ~= 3
                error('Sensor Position Array does not have three dimensions');
            end
            obj.N = uint16(size(Sensors,1));
            obj.Sensors = Sensors;
            
        end
        
        function str = toString(obj)
            str = sprintf('%.1fMHz @ %i%%, %.2fpi',...
                obj.f_center*1e-6,...
                obj.bw * 100,...
                obj.phase_shift);
        end
        
        
        
        
       % *****************************************************************
       %% CALCULATIONS
        function obj = calculateSensors(obj)
            if strcmp(obj.Type,'Coordinates')
                error('Cannot Calculate Sensors for Coordinate Definition');
                return;
            end
            
            fprintf('Calculating Cartesian Coordinates from Spherical Coordinates\n');
            pos_elements=zeros(obj.N,3);
            [pos_elements(:,1),pos_elements(:,3),pos_elements(:,2)]=sph2cart(zeros(obj.N,1),obj.Angles',ones(obj.N,1)*obj.Radius);
            obj.Sensors = pos_elements;
        end

        
        % *****************************************************************
       %% CONFIGURATION
       function configure(obj)
           ret = inputdlg({'Center Frequency (MHz)','Bandwidth (% T/R)','Phase Offset'},'Impulse Response',1,...
                { num2str(obj.f_center*1e-6,'%.1f'),...
                  num2str(round(obj.bw*1e2),'%i'),...
                  num2str(obj.phase_shift,'%.2f')   });
           if isempty(ret), return; end;
           fc = str2double(ret{1}); if isnan(fc) || ~(fc > 1), errordlg('Cannot convert Frequency %s to a number (%.1f), aborting',ret{1},fc); return; end;
           baw = str2double(ret{2}); if isnan(baw) ||~(baw > 1), errordlg('Cannot convert bandwidth %s to a number (%.1f), aborting',ret{2},baw); return; end;
           ph = str2double(ret{3}); if isnan(ph), errordlg('Cannot convert phase %s to a number (%.2f), aborting',ret{3},ph); return; end;
           obj.f_center = fc*1e6;
           obj.bw = baw*1e-2;
           obj.phase_shift = ph;
       end
        
        
        
        % *****************************************************************
       %% FILE ACCESS
        
        function filename = write(obj,filename)
            % write   Write Time Domain IR to File
            % Writes Impulse Response to file compatible with ViewMSOT and
            % ViewMSOTc. If no parameters are given, a generic file name
            % based on the properties is used and returned
            %
            % Parameters:
            % (1): Filename [optional]

            if nargin == 1 || isempty(filename),
                filename = ['IR_' strrep(sprintf('%.1f',obj.f_center*1e-6),'.','c') 'MHz_'...
                            num2str(round(obj.bw*100)) 'perc' '.irf'];
            end

            FID = fopen(filename,'w');
            fwrite(FID,obj.tdom,'double');
            fclose(FID);
            fprintf('Impulse response written to %s\n',filename);
        end
        
    end
end