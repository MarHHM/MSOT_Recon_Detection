%% reconstruct single scan (conv MB)  
fprintf("\n------------ Script 'reconstruct single scan (conv MB)' ------------\r");

scan_path = [datasets__path '\Marwan\(2017-03-03) testing phantom ink+straw\Scan_16'];    %  || Marwan\(2017-03-03) testing phantom ink+straw\Scan_16 || \Qutaiba\phantom_3\Scan_2 || '\Qutaiba\Study_4\Scan_4'
% scan_path_parts = strsplit(scan_path, '\');
% datainfo = loadMSOT((scan_path+"\"+scan_path_parts(end)+".msot"));
nonNeg = [
%           false;
          true;
                ];
RECON_ALL = false;       % !! if true, 'recon_lims' is ignored!!
% recon_lims.zpos= [1 length(datainfo.ZPositions)];
% recon_lims.wls = [1 length(datainfo.Wavelengths)];
recon_lims.zpos = [15 15];
recon_lims.wls = [1 1];
addtnl_note = "";       % additional note to be added at the end of the recon file name (e.g. --noReg || --f_min_0 || --f_max_5)

for i = 1:length(nonNeg)
  recon_single_msot_scan(scan_path, nonNeg(i), recon_lims, RECON_ALL, addtnl_note);
end

%% showing single recon (conv MB)
disp('------------ Script "showing single recon (conv MB)" ------------');
study = "Marwan\(2017-03-03) testing phantom ink+straw";     %"Marwan\(2017-03-03) testing phantom ink+straw" || "Qutaiba\phantom_2" ||
scan = "Scan_16";
recon = "MB_Tik-nonNeg_1-zposs_1-wls_1--thetaCvrg_180";
recon_path = datasets__path+"\"+study+"\"+scan+"\recons\"+recon+"\"+recon+".mat";
% zpos = 59.6200;
% wl = 850;
zpos_idx = 1;
wl_idx = 1;


disp("Loading recon structure..");
load(recon_path);
disp_last_line("done.");
% wl_idx = find(datainfo.Wavelengths==wl);
% zpos_idx = find(datainfo.ZPositions==zpos);
figure, imagesc(Recon(:,:,1,zpos_idx,1,wl_idx));
    title(("conv MB"+" (recon: "+study+"\"+scan+"\"+recon+") @ zpos="+zpos+" & wl="+wl), 'interpreter', 'none'), colormap(bone), colorbar, axis image off
    
%% wMB --> reconstruct single scan
fprintf("\n------------ Script 'wMB --> reconstruct single scan' ------------\r");

scan_path = datasets__path+"\Marwan\(2017-03-03) testing phantom ink+straw\Scan_16";    % \Qutaiba\phantom_3\Scan_2 || Marwan\(2017-03-03) testing phantom ink+straw\Scan_16 ||
base_convMB_recon__nonNeg = [
%                              false;
                             true
                                    ];
LVc = [
%        0;       % original wMB (Luis 2013)
       1        % my method
         ];
RECON_ALL = false;      % !! if true, 'recon_lims' is ignored!!
% recon_lims.zpos= [1 length(datainfo.ZPositions)];
% recon_lims.wls = [1 length(datainfo.Wavelengths)];
recon_lims.zpos = [15 15];
recon_lims.wls = [1 1];
addtnl_note = "";       % additional note to be added at the end of the recon file name (e.g. --noReg || --f_min_0 || --f_max_5)

for i = 1:length(base_convMB_recon__nonNeg)
    for j = 1:length(LVc)
      base_convMB_recon__fName = ("MB_Tik-nonNeg_"+ int2str(base_convMB_recon__nonNeg(i)) +"-zposs_20-wls_8"+addtnl_note);
      recon_wMB(scan_path, base_convMB_recon__nonNeg(i), LVc(j), RECON_ALL, recon_lims, base_convMB_recon__fName, addtnl_note);
    end
end

%% showing single recon (wMB)
disp('------------ Script "showing single recon (wMB)" ------------');
study = "Qutaiba\phantom_3";     %"Marwan\(2017-03-03) testing phantom ink+straw" || "Qutaiba\phantom_2" ||
scan = "Scan_2";
recon = "wMB - nonNeg_0-zpos_1-reps_1-wls_1-reconRes_7.5e-05-imW_0.03-w_1-LVc_1shift0";
recon_path = datasets__path+"\"+study+"\"+scan+"\RECONs\"+recon+".mat";

load(recon_path)
figure, imagesc(ReconW(:,:,1,1,1,1)),...
    title("wMB" + " (scan: " + study+"\"+scan + ")"), colormap(bone), colorbar, axis image off
%% Reconstructing several scans
disp('------------ Script "Reconstructing several scans" ------------');
dataset_path = datasets__path+"\Qutaiba phantom 2";
scan_names = ["Scan_1"; "Scan_2";];
nonNeg = [false; true];
RECON_ALL = true;

for i = 1:size(scan_names, 1)
  for j = 1:length(nonNeg)
%                 disp([[dataset_path "\" scan_names(i,:)] int2str(nonNeg(j))])
    recon_single_msot_scan([dataset_path "\" scan_names(i,:)], nonNeg(j), RECON_ALL);       % \Hong\phantom\Scan_69  ||  \Marwan\(2017-03-03) testing phantom ink+straw\both horizontal || \Qutaiba phantom\Scan_1
  end
end

%% showing all wls @ certain zpos
disp('------------ Script "showing all wls @ certain zpos" ------------');
% close all;
% scan_path = [datasets__path '\Marwan\(2017-03-03) testing phantom ink+straw\Scan_16'];    % \Qutaiba phantom\Scan_5 || Marwan\(2017-03-03) testing phantom ink+straw\Scan_16
zpos = 1;
study = "Qutaiba\Study_4\";     %"Marwan\(2017-03-03) testing phantom ink+straw" || "Qutaiba phantom 2" ||
scan = "Scan_4";
recon_fName = "MB_Tik-nonNeg_0-zposs_8-wls_31";
recon_path = datasets__path+"\"+study+"\"+scan+"\recons\"+recon_fName+"\"+recon_fName+".mat";

load(recon_path)
for wl_idx = 1:length(datainfo.Wavelengths)
  figure, imagesc(Recon(:,:,1,zpos,1,wl_idx)),...
    title("conv MB - wl = "+ int2str(datainfo.Wavelengths(wl_idx)) + " (scan: " + study+scan + ")"), colormap(bone), colorbar, axis image off
end

%% adapt hist eq to enahnce contrast of Qutaib phantom_2 & show borders
% disp('------------ Script "adapt hist eq to enahnce contrast of Qutaiba phantom_2 & show borders" ------------');
% study = "Qutaiba\phantom_2";
% scan = "Scan_5";
% recon_to_show = "MB_Tik - nonNeg_0-zpos_1-reps_1-wls_11-reconRes_7.5e-05-imW_0.03";
% recon_path = datasets__path+"\"+study+"\"+scan+"\recons\"+recon_to_show+".mat";
% 
% load(recon_path)
% test_im = mat2gray(Recon(:,:,1,1,1,2));    % normalize image before histogram eqlztn
% test_im_enhncd = adapthisteq(test_im, "NumTiles", [4,4]);
% figure, montage({test_im, test_im_enhncd});
% % figure, imagesc(test_im_enhncd), colormap(bone), colorbar, axis image off

%% for wl_idx -> show all zpos (cross-sections)
disp('------------ Script "for wl_idx -> show all zpos (cross-sections)" ------------');
wl_idx = 1;
for zpos_idx = 1:length(datainfo.ZPositions)
  figure, imagesc(Recon(:,:,1,zpos_idx,1,wl_idx)), title(["conv MB - wl = " int2str(datainfo.Wavelengths(wl_idx)) " - zpos = " int2str(datainfo.ZPositions(zpos_idx))]), colormap(bone), colorbar, axis image off
end

%% Reconstructing several scans (wMB)
disp('------------ Script "Reconstructing several scans (wMB)" ------------');
dataset_path = datasets__path+"\Qutaiba phantom 2";
scan_names = ["Scan_2"; "Scan_5"];
RECON_ALL = 1;
nonNeg = [0; 1];
LVc = [0; 1];
shift = [0];
% shift = [5:5:50];

for j = 1:size(scan_names, 1)
  for j = 1:length(nonNeg)
    for k = 1:length(LVc)
      for l = 1:length(shift)
        %             disp([[dataset_path "\" scan_names(i,:)] int2str(nonNeg(j))  int2str(LVc(k))])
        recon_wMB((dataset_path + "\" + scan_names(j,:)), nonNeg(j), LVc(k), RECON_ALL,  shift(l));       % \Hong\phantom\Scan_69  ||  \Marwan\(2017-03-03) testing phantom ink+straw\both horizontal || \Qutaiba phantom\Scan_1
        disp("-------------------------------------------------------------------------------");
      end
    end
  end
end


%% Unmix several scans & recons
disp('------------ Script "Unmix several scans & recons" ------------');
dataset_path = datasets__path+"\Marwan\(2017-03-03) testing phantom ink+straw\";    % Marwan\(2017-03-03) testing phantom ink+straw\Scan_16 || Qutaiba\phantom_2
scan_names = ["Scan_16";];
agent = "ink";     % ICG || iRFP || AF750 || ink
z_pos__idx = 10;
NORMALIZE_INPUT_SPECTRA = true;     % true-> binary detection || false-> concentration detection (better, specially if unmixing for only one agent)
% recons_to_unmix = ["Scan_8";];
% nonNeg = [0;];
% LVc = [1];
% for i = 1:size(scan_names, 1)
%   for j = 1:length(nonNeg)
%     for k = 1:length(LVc)
%       for l = 1:length(shift)
        %             disp([[dataset_path "\" scan_names(i,:)] int2str(nonNeg(j))  int2str(LVc(k))])
recon_fName = "MB_Tik-nonNeg_0-zposs_20-wls_8";
spectralAnalysis(agent, (dataset_path+"\"+scan_names+"\recons\"+recon_fName+"\"+recon_fName+".mat"),...
                  "MB_Tik", z_pos__idx, MSOT_Recon_Detection__path, NORMALIZE_INPUT_SPECTRA);
% spectralAnalysis([dataset_path "\" scan_names "\recons\wMB - nonNeg_0-zpos_8-reps_1-wls_13-reconRes_7.5e-05-imW_0.03-w_1-LVc_0.mat"],...
%                   "wMB", z_pos__idx);
% spectralAnalysis([dataset_path "\" scan_names "\recons\wMB - nonNeg_0-zpos_8-reps_1-wls_13-reconRes_7.5e-05-imW_0.03-w_1-LVc_1.mat"],...
%                   "wMB", z_pos__idx);
%       end
%     end
%   end
% end

%%