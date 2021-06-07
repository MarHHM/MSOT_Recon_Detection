% executed when Matlab starts (rename to startup.m & put in C:\Users\ibn_a\Documents\MATLAB)

%% PATHS 
onedrive__path = "C:\Users\ibn_a\OneDrive\";
MB_MAT__path = "C:\PA_local\MB_matrices";
datasets__path =  'C:\PA_local\DATASETS';
MSOT_Recon_Detection__path = onedrive__path+"\PA_imaging\wMB+LVc\MSOT_Recon_Detection";

addpath(genpath(onedrive__path+"\Others\My_Matlab_toolBox"));
addpath(genpath(onedrive__path+"\PA_imaging\wMB+LVc\MSOT_Recon_Detection\"));


% addpath(genpath(MSOT_Recon_Detection__path));
% addpath(genpath('D:\OneDrive\PA_imaging\msotrecon_handheld'));
% addpath(genpath('D:\OneDrive\PA_imaging\Luis code\'));
% addpath(genpath('D:\OneDrive\PA_imaging\MSOT_Analysis_GUI (IBMI gui for recon & analysis)'));
% addpath(genpath('D:\OneDrive\PA_imaging\acousticSimWrapper-master\'));

javaaddpath(MSOT_Recon_Detection__path+"\MSOTBeans\msotbeans.jar");
javaaddpath(MSOT_Recon_Detection__path+"\java_class\recon.jar");
javaaddpath(MSOT_Recon_Detection__path+"\java_class\patbeans.jar");
javaaddpath(MSOT_Recon_Detection__path+"\java_class\xmlbeans-2.5.0\lib\xbean.jar");

%% cosmetics
% set(groot,'defaultLineLineWidth', 3);
% set(groot,'defaultLineMarkerSize', 6);
% % set(groot,'defaultAxesFontSize', 18);            % Jaber: 16 to 28 are good for publications
% % set(groot,'defaultTextColor', [0.5 0.5 0.5]);           % to be apropriate for both white & black bgds (used in ppt)
% set(groot,'defaultTextColor', 'k');
% set(groot,'defaultAxesFontWeight','bold');
% set(groot,'defaultColorbarTickLength', 0.03);
% set(groot,'defaultColorbarFontSize', 15);         % NOT WORKING!!
% set(groot,'defaultColorbarLineWidth', 2);

%% for the 'Matlab-Editor-Plugin' at 'https://github.com/GavriYashar/Matlab-Editor-Plugin/wiki/Setup'
% at.mep.Start.start('C:\Users\marwan.muhammad\Google Drive\My_Matlab_toolBox\Matlab-Editor-Plugin\CustomProps.properties',...
%   'C:\Users\marwan.muhammad\Google Drive\My_Matlab_toolBox\Matlab-Editor-Plugin\DefaultProps.properties');
% import at.mep.util.*;
% import at.mep.editor.*;
% 
% ctrl = true;
% shift = true;
% alt = true;
% 
% ks_MyKeyReleaseFunction = KeyStrokeUtil.getKeyStroke(java.awt.event.KeyEvent.VK_INSERT, ~ctrl, ~shift, alt, false);
% EditorApp.addMatlabCallback('MyKeyReleaseFunction', ks_MyKeyReleaseFunction , 'MyKeyReleaseFunction');
