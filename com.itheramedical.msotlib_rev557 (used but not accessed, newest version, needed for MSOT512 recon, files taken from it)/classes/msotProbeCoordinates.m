classdef msotProbeCoordinates < msotProbe
    properties (Constant = true)
        Type@char vector = 'Coordinates';
    end

    properties (SetObservable)
        Sensors@double = zeros(0,3);
    end
    
    properties
        StartAngle@double = nan;                % Angle of first Element, read only
        StepAngle@double = nan;                 % Step Angle, read only 
        EndAngle@double = nan;                  % End Angle, read only 
        Pitch@double = nan;                    % Width of an Element, in Plane (only Parametric)
    end
    
    properties (Dependent = true)
        Radius@double;
        Coverage@double;
        NumDetectors@uint16;
        Angles@double;
        Spherical@double;
    end
  
methods
    % *****************************************************************
    %% CONSTRUCTOR
    function obj = msotProbeCoordinates(ID,DID,S,varargin)
        % Create MSOT Coordinate Probe Object
        %
        % msotProbeParametric(ProbeID,DesignID,Radius,Coverage,N)
        obj@msotProbe(ID,DID,varargin{:});
        obj.Sensors = S;
    end
    
    
    % *****************************************************************
    %% ACCESSORS
    function R = get.Radius(obj)
        if ~obj.muted, warning('Radius may be inaccurate for Probes with Coordinate Definition'); end
        R=mean(obj.Spherical(:,3));
    end

    function C = get.Coverage(obj)
        if isempty(obj.Spherical), C = 0; return; end
        if ~obj.muted, warning('Coverage may be inaccurate for Probes with Coordinate Definition'); end
        ang_=obj.Spherical(:,2);
        C = abs(diff([min(ang_),max(ang_)])) / pi * 180;
        C = C*(obj.N+1)/obj.N;
    end

    function N = get.NumDetectors(obj)
        N = uint16(size(obj.Sensors,1));
    end

    function A = get.Angles(obj) 
        if ~obj.muted, warning('Angles may be inaccurate for Probes with Coordinate Definition'); end        
        A = obj.Spherical(:,2);
    end
          
    function obj = set.Sensors(obj,Sensors)
        if isempty(Sensors), obj.Sensors = []; return; end;
        if size(Sensors,2) ~= 3
            error('Sensor Position Array does not have three dimensions');
        end
        obj.Sensors = Sensors;
        if ~obj.isSymmetric
            warning('Provided Geometry is not symmetric');
        end
    end
       
    function S = get.Spherical(obj)
        Sens = obj.Sensors;
        [S(:,1), S(:,2), S(:,3)] = cart2sph(Sens(:,1),Sens(:,3),Sens(:,2));
        loc = S(:,1) == pi;
        S(loc,2) = -(pi+S(loc,2));
        S(loc,1) = 0;
    end
    
    
    
    % *****************************************************************
    %% FUNCTIONALITY
    function str = writeXML_(~)
        str = [];
    end
    
    function obj = readXML_(obj,~)
    end
  end
end