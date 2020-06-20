%% Reconstructing several scans
dataset_path = 'S:\PA_DATASETS\MSOT 256\Qutaiba phantom';
scan_names = ['Scan_8';];
for i = 1:size(scan_names, 1)
    recon_single_msot_scan([dataset_path '\' scan_names(i,:)]);       % \Hong\phantom\Scan_69  ||  \Marwan\(2017-03-03) testing phantom ink+straw\both horizontal || \Qutaiba phantom\Scan_1
end

%% reconstruct single scan (conv MB)
scan_path = 'S:\PA_DATASETS\MSOT 256\Marwan\(2017-03-03) testing phantom ink+straw\Scan_16';    % \Qutaiba phantom\Scan_5 || Marwan\(2017-03-03) testing phantom ink+straw\Scan_16
nonNeg = false;
RECON_ALL = false;
recon_single_msot_scan(scan_path, nonNeg, RECON_ALL);

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
nonNeg = [1;];
LVc = [0; 1];
         
for i = 1:size(scan_names, 1)
    for j = 1:length(nonNeg)
        for k = 1:length(LVc)
%             disp([[dataset_path '\' scan_names(i,:)] int2str(nonNeg(j))  int2str(LVc(k))])
            recon_wMB([dataset_path '\' scan_names(i,:)], nonNeg(j), LVc(k));       % \Hong\phantom\Scan_69  ||  \Marwan\(2017-03-03) testing phantom ink+straw\both horizontal || \Qutaiba phantom\Scan_1
        end
    end
end