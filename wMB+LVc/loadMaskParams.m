function [xc, yc, Rc, logicalMask] = loadMaskParams(targetPath, mask)

% switch targetPath
%   case 'D:\DATASETS\MSOT 256\Marwan\(2017-03-03) testing phantom ink+straw\both horizontal\'
%     %%% dims for my straw phantom in MSOT 256 - both horizontal (true spectrum)
%     xc = -0.4e-3;   % centre of the phantom in x direction (w.r.t to whole image centers)
%     yc = -1.8e-3;  % centre of the phantom in y direction
%     Rc = 10.3e-3;    % radius of the phantom (m) (phantom was done inside the 20ml syringe whose diam = 20mm)
%   case 'D:\DATASETS\MSOT 256\Marwan\(2017-03-03) testing phantom ink+straw\straw up\'
%     if mask == 'A'
%       xc = -0.000305259384889830;
%       yc = -0.000203950338600455;
%       Rc = 0.0102063046023447;
%     else          % 'B' --> reflector ROI
%       xc = 0.00213280451237456;
%       yc = -0.00638530997304583;
%       Rc = 0.00212130339074157;
%     end
%       
%   case 'D:\DATASETS\MSOT 256\Scan_96 (Ivan - cleanest data with tumor)\'
%     xc = -0.5e-3;
%     yc = 0e-3;
%     Rc = 9.3e-3;    % radius of the target (m) (nominally 8e-3)
%   case 'D:\DATASETS\MSOT 256\Evangelos (high reflections)\Scan_81\'
%     xc = 0.000;
%     yc = -.002;
%     Rc = 0.01;
%   case 'D:\DATASETS\MSOT 256\Evangelos (high reflections)\Scan_90\'
%     xc = 2.106481481481524e-05;
%     yc = -7.686397984886611e-04;
%     Rc = 0.011500000000000;
%   case 'Evangelos - scan 94'
%     xc = 8.891203703703697e-04;
%     yc = -5.797229219143560e-04;
%     Rc = 0.011500000000000;
%   case 'D:\DATASETS\MSOT 256\Brain (Ivan)\study 2\Scan_65\'
%     xc = 2.685185185185333e-05;
%     yc = -0.002267380352645;
%     Rc = 0.0085;
%   case 'D:\DATASETS\MSOT 256\Scan_26 (Josefine - reflections from black spots)\Scan_26\'
%     xc = 0.000340772038885246;
%     yc = -0.00246650558914710;
%     Rc = 0.0095;
%   case 'D:\DATASETS\MSOT 256\Qutaiba blood phantom\Scan_6\'
%     xc = -0.000388005390835576;
%     yc = 0.000252156334231803;
%     Rc = 0.012;
% end

%%% dims for my straw phantom in MSOT 256 (false spectrum)
% xc = -0.5e-3;   % centre of the phantom in x direction (w.r.t to whole image centers)
% yc = 0;  % centre of the phantom in y direction
% Rc = 10.3e-3;    % radius of the phantom (m) (phantom was done inside the 20ml syringe whose diam = 20mm)
%%% dims for my straw phantom in MSOT 256 - ink up (true spectrum)
% xc = -0.7e-3;   % centre of the phantom in x direction (w.r.t to whole image centers)
% yc = -0.4e-3;  % centre of the phantom in y direction
% Rc = 10.3e-3;    % radius of the phantom (m) (phantom was done inside the 20ml syringe whose diam = 20mm)
%%% dims for my straw phantom in MSOT 256 - straw up (true spectrum)
% xc = -0.3e-3;   % centre of the phantom in x direction (w.r.t to whole image centers)
% yc = -0.1e-3;  % centre of the phantom in y direction
% Rc = 10.3e-3;    % radius of the phantom (m) (phantom was done inside the 20ml syringe whose diam = 20mm)

%%% dims for Steven mouse data from MSOT 512
% xc = 1e-3;   % centre of the phantom in x direction (w.r.t to whole image centers)
% yc = 1.3e-3;  % centre of the phantom in y direction
% Rc = 8e-3;    % radius of the target (m) (nominally 8e-3)

