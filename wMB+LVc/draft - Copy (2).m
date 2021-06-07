%% reconstruct single scan (conv MB)
scan_path = 'S:\PA_DATASETS\MSOT 256\Marwan\(2017-03-03) testing phantom ink+straw\Scan_16';    % \Qutaiba phantom\Scan_5 || Marwan\(2017-03-03) testing phantom ink+straw\Scan_16
nonNeg = false;
RECON_ALL = false;
recon_single_msot_scan(scan_path, nonNeg, RECON_ALL);

%% Reconstructing several scans
dataset_path = 'S:\PA_DATASETS\MSOT 256\Qutaiba phantom 2';
scan_names = ['Scan_5'; 'Scan_6';];
nonNeg = [false; true];
RECON_ALL = true;

for i = 1:size(scan_names, 1)
  for j = 1:length(nonNeg)
%                 disp([[dataset_path '\' scan_names(i,:)] int2str(nonNeg(j))])
    recon_single_msot_scan([dataset_path '\' scan_names(i,:)], nonNeg(j), RECON_ALL);       % \Hong\phantom\Scan_69  ||  \Marwan\(2017-03-03) testing phantom ink+straw\both horizontal || \Qutaiba phantom\Scan_1
  end
end

%% reconstruct single scan (wMB)
scan_path = 'S:\PA_DATASETS\MSOT 256\Marwan\(2017-03-03) testing phantom ink+straw\Scan_16';    % \Qutaiba phantom\Scan_5 || Marwan\(2017-03-03) testing phantom ink+straw\Scan_16
nonNeg = false;
LVc = 0;
RECON_ALL = false;
recon_wMB(scan_path, nonNeg, LVc, RECON_ALL);

% LVc = 1;
% recon_wMB(scan_path, nonNeg, LVc, RECON_ALL);
%% showing all wls at the first zpos
for wl_idx = 1:length(datainfo.Wavelengths)
  figure, imagesc(Recon_MB(:,:,1,1,1,wl_idx)), title(['conv MB - wl = ' int2str(datainfo.Wavelengths(wl_idx))]), colormap(bone), colorbar, axis image off
end

%% for wl_idx -> show all zpos (cross-sections)
wl_idx = 1;
for zpos_idx = 1:length(datainfo.ZPositions)
  figure, imagesc(Recon_MB(:,:,1,zpos_idx,1,wl_idx)), title(['conv MB - wl = ' int2str(datainfo.Wavelengths(wl_idx)) ' - zpos = ' int2str(datainfo.ZPositions(zpos_idx))]), colormap(bone), colorbar, axis image off
end

%% Reconstructing several scans (wMB)
dataset_path = 'S:\PA_DATASETS\MSOT 256\Qutaiba phantom';
scan_names = ['Scan_8';];
RECON_ALL = 0;
nonNeg = [0;];
LVc = [1];
shift = [5:5:50];

for i = 1:size(scan_names, 1)
  for j = 1:length(nonNeg)
    for k = 1:length(LVc)
      for l = 1:length(shift)
        %             disp([[dataset_path '\' scan_names(i,:)] int2str(nonNeg(j))  int2str(LVc(k))])
        recon_wMB([dataset_path '\' scan_names(i,:)], nonNeg(j), LVc(k), RECON_ALL,  shift(l));       % \Hong\phantom\Scan_69  ||  \Marwan\(2017-03-03) testing phantom ink+straw\both horizontal || \Qutaiba phantom\Scan_1
      end
    end
  end
end


%% Unmix several scans & recons
dataset_path = 'S:\PA_DATASETS\MSOT 256\Qutaiba phantom';
scan_names = ['Scan_8';];
z_pos__idx = 8;
% recons_to_unmix = ['Scan_8';];
% nonNeg = [0;];
% LVc = [1];
% for i = 1:size(scan_names, 1)
%   for j = 1:length(nonNeg)
%     for k = 1:length(LVc)
%       for l = 1:length(shift)
        %             disp([[dataset_path '\' scan_names(i,:)] int2str(nonNeg(j))  int2str(LVc(k))])
spectralAnalysis([dataset_path '\' scan_names '\recons\MB_Tik - nonNeg_0-zpos_8-reps_1-wls_13-reconRes_7.5e-05-imW_0.03.mat'],...
                  'MB_Tik', z_pos__idx);
spectralAnalysis([dataset_path '\' scan_names '\recons\wMB - nonNeg_0-zpos_8-reps_1-wls_13-reconRes_7.5e-05-imW_0.03-w_1-LVc_0.mat'],...
                  'wMB', z_pos__idx);
spectralAnalysis([dataset_path '\' scan_names '\recons\wMB - nonNeg_0-zpos_8-reps_1-wls_13-reconRes_7.5e-05-imW_0.03-w_1-LVc_1.mat'],...
                  'wMB', z_pos__idx);
%       end
%     end
%   end
% end

%%