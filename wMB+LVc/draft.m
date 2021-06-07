%% reconstruct single scan (conv MB)
disp('------------ Script "reconstruct single scan (conv MB)" ------------');
scan_path = [datasets__path '\MSOT 256\Qutaiba\phantom_1\Scan_8'];    %  || Marwan\(2017-03-03) testing phantom ink+straw\Scan_16
nonNeg = false;
RECON_ALL = false;
recon_single_msot_scan(scan_path, nonNeg, RECON_ALL);

%% showing single recon (conv MB)
disp('------------ Script "showing single recon (conv MB)" ------------');
study = "Qutaiba\phantom_1";     %"Marwan\(2017-03-03) testing phantom ink+straw" || "Qutaiba\phantom_2" ||
scan = "Scan_2";
recon_to_show = "MB_Tik - nonNeg_0-zpos_6-reps_1-wls_55-reconRes_7.5e-05-imW_0.03";
recon_path = datasets__path+"\MSOT 256\"+study+"\"+scan+"\RECONs\"+recon_to_show+".mat";

load(recon_path)
figure, imagesc(Recon_MB(:,:,1,1,1,1)),...
    title("conv MB" + " (scan: " + study+"\"+scan + ")"), colormap(bone), colorbar, axis image off

%% showing single recon (wMB)
disp('------------ Script "showing single recon (wMB)" ------------');
study = "Qutaiba\phantom_2";     %"Marwan\(2017-03-03) testing phantom ink+straw" || "Qutaiba\phantom_2" ||
scan = "Scan_3";
recon_to_show = "wMB - nonNeg_0-zpos_1-reps_1-wls_11-reconRes_7.5e-05-imW_0.03-w_1-LVc_0shift0";
recon_path = datasets__path+"\MSOT 256\"+study+"\"+scan+"\RECONs\"+recon_to_show+".mat";

load(recon_path)
figure, imagesc(ReconW(:,:,1,1,1,1)),...
    title("wMB" + " (scan: " + study+"\"+scan + ")"), colormap(bone), colorbar, axis image off
%% Reconstructing several scans
disp('------------ Script "Reconstructing several scans" ------------');
dataset_path = datasets__path+"\MSOT 256\Qutaiba phantom 2";
scan_names = ["Scan_1"; "Scan_2";];
nonNeg = [false; true];
RECON_ALL = true;

for i = 1:size(scan_names, 1)
  for j = 1:length(nonNeg)
%                 disp([[dataset_path "\" scan_names(i,:)] int2str(nonNeg(j))])
    recon_single_msot_scan([dataset_path "\" scan_names(i,:)], nonNeg(j), RECON_ALL);       % \Hong\phantom\Scan_69  ||  \Marwan\(2017-03-03) testing phantom ink+straw\both horizontal || \Qutaiba phantom\Scan_1
  end
end

%% reconstruct single scan (wMB)
disp('------------ Script "reconstruct single scan (wMB)" ------------');
scan_path = datasets__path+"\MSOT 256\Marwan\(2017-03-03) testing phantom ink+straw\Scan_16";    % \Qutaiba phantom\Scan_5 || Marwan\(2017-03-03) testing phantom ink+straw\Scan_16
nonNeg = false;
LVc = 0;
RECON_ALL = false;
recon_wMB(scan_path, nonNeg, LVc, RECON_ALL);

% LVc = 1;
% recon_wMB(scan_path, nonNeg, LVc, RECON_ALL); 
%% showing all wls @ certain zpos
disp('------------ Script "showing all wls @ certain zpos" ------------');
% close all;
% scan_path = [datasets__path '\MSOT 256\Marwan\(2017-03-03) testing phantom ink+straw\Scan_16'];    % \Qutaiba phantom\Scan_5 || Marwan\(2017-03-03) testing phantom ink+straw\Scan_16
zpos = 1;
study = "Marwan\(2017-03-03) testing phantom ink+straw";     %"Marwan\(2017-03-03) testing phantom ink+straw" || "Qutaiba phantom 2" ||
scan = "Scan_16";
recon_path = datasets__path+"\MSOT 256\"+study+"\"+scan+"\RECONs\MB_Tik - nonNeg_0-zpos_1-reps_1-wls_1-reconRes_7.5e-05-imW_0.03__hhhhhhh.mat";

load(recon_path)
for wl_idx = 1:length(datainfo.Wavelengths)
  figure, imagesc(Recon_MB(:,:,1,zpos,1,wl_idx)),...
    title("conv MB - wl = "+ int2str(datainfo.Wavelengths(wl_idx)) + " (scan: " + study+"\"+scan + ")"), colormap(bone), colorbar, axis image off
end

%% adapt hist eq to enahnce contrast of Qutaib phantom_2 & show borders
disp('------------ Script "adapt hist eq to enahnce contrast of Qutaiba phantom_2 & show borders" ------------');
study = "Qutaiba\phantom_2";
scan = "Scan_5";
recon_to_show = "MB_Tik - nonNeg_0-zpos_1-reps_1-wls_11-reconRes_7.5e-05-imW_0.03";
recon_path = datasets__path+"\MSOT 256\"+study+"\"+scan+"\recons\"+recon_to_show+".mat";

load(recon_path)
test_im = mat2gray(Recon_MB(:,:,1,1,1,2));    % normalize image before histogram eqlztn
test_im_enhncd = adapthisteq(test_im, "NumTiles", [4,4]);
figure, montage({test_im, test_im_enhncd});
% figure, imagesc(test_im_enhncd), colormap(bone), colorbar, axis image off

%% for wl_idx -> show all zpos (cross-sections)
disp('------------ Script "for wl_idx -> show all zpos (cross-sections)" ------------');
wl_idx = 1;
for zpos_idx = 1:length(datainfo.ZPositions)
  figure, imagesc(Recon_MB(:,:,1,zpos_idx,1,wl_idx)), title(["conv MB - wl = " int2str(datainfo.Wavelengths(wl_idx)) " - zpos = " int2str(datainfo.ZPositions(zpos_idx))]), colormap(bone), colorbar, axis image off
end

%% Reconstructing several scans (wMB)
disp('------------ Script "Reconstructing several scans (wMB)" ------------');
dataset_path = datasets__path+"\MSOT 256\Qutaiba phantom 2";
scan_names = ["Scan_2"; "Scan_5"];
RECON_ALL = 1;
nonNeg = [0; 1];
LVc = [0; 1];
shift = [0];
% shift = [5:5:50];

for i = 1:size(scan_names, 1)
  for j = 1:length(nonNeg)
    for k = 1:length(LVc)
      for l = 1:length(shift)
        %             disp([[dataset_path "\" scan_names(i,:)] int2str(nonNeg(j))  int2str(LVc(k))])
        recon_wMB((dataset_path + "\" + scan_names(i,:)), nonNeg(j), LVc(k), RECON_ALL,  shift(l));       % \Hong\phantom\Scan_69  ||  \Marwan\(2017-03-03) testing phantom ink+straw\both horizontal || \Qutaiba phantom\Scan_1
        disp("-------------------------------------------------------------------------------");
      end
    end
  end
end


%% Unmix several scans & recons
disp('------------ Script "Unmix several scans & recons" ------------');
dataset_path = datasets__path+"\MSOT 256\Qutaiba phantom";
scan_names = ["Scan_8";];
z_pos__idx = 8;
% recons_to_unmix = ["Scan_8";];
% nonNeg = [0;];
% LVc = [1];
% for i = 1:size(scan_names, 1)
%   for j = 1:length(nonNeg)
%     for k = 1:length(LVc)
%       for l = 1:length(shift)
        %             disp([[dataset_path "\" scan_names(i,:)] int2str(nonNeg(j))  int2str(LVc(k))])
spectralAnalysis([dataset_path "\" scan_names "\recons\MB_Tik - nonNeg_0-zpos_8-reps_1-wls_13-reconRes_7.5e-05-imW_0.03.mat"],...
                  "MB_Tik", z_pos__idx);
spectralAnalysis([dataset_path "\" scan_names "\recons\wMB - nonNeg_0-zpos_8-reps_1-wls_13-reconRes_7.5e-05-imW_0.03-w_1-LVc_0.mat"],...
                  "wMB", z_pos__idx);
spectralAnalysis([dataset_path "\" scan_names "\recons\wMB - nonNeg_0-zpos_8-reps_1-wls_13-reconRes_7.5e-05-imW_0.03-w_1-LVc_1.mat"],...
                  "wMB", z_pos__idx);
%       end
%     end
%   end
% end

%%