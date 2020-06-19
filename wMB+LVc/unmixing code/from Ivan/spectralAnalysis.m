%% NOTES
% - choose ur datacube for which u will do the unmixing
% - the datacube DATA should be arranged as a 3D array of (y,x,wls)
% - ask the source of the dataset if any agents were used (for the A_est choice)

%% PATHS & PARAMS
SHIFT_RECONS = false;       % to avoid neg values (make min = 0) - i think it should be true if negatives were allowed during the reconstruction

%% %%%%%%%%%%%%%%%%%%%%%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% choose input recon (between convMB or ReconW)
DATA = squeeze(Recon_MB);
% DATA = squeeze(ReconW);

%% load the spectra
spectraPath = 'C:\Users\bzfmuham\OneDrive\PA_imaging\wMB+LVc\MSOT_Recon_Detection\spectra\SpectralSpecifications_AF750';  % \SpectralSpecifications_AF750 || \SpectralSpecifications_iCG
wavelengths = datainfo.Wavelengths;
spectra = LoadSpectra(spectraPath, wavelengths);
water = datainfo.MeasurementDesc.WaterAbsorptionCoeff;
a = 18.2; b = 0.6;
melanin = a.*(wavelengths./500).^(-b);
fat_spectra = load('fat_spectra.mat');
lipid = fat_spectra.fat( ismember(fat_spectra.fat(:,1), wavelengths), 2);

%% read data & spectra for wavelengths<=900 (to avoid dominance of water & fat signals in NIR)
Recon_MB_so2 = DATA(:,:,wavelengths<=900);
wavelengths_so2 = wavelengths(wavelengths<=900);
water_so2 = water(wavelengths<=900);
spectra_so2 = spectra(:, wavelengths<=900);
melanin_so2 = melanin(wavelengths<=900);
lipid_so2 = lipid(wavelengths<=900);

if SHIFT_RECONS
    Recon_MB_so2(:,:,1:end) = Recon_MB_so2(:,:,1:end) - min(min(Recon_MB_so2(:,:,1:end)));      % vectorized - instead of using a 'for' loop
end

%% set up the absorbers matrix (choose spectra according to what chromophores are significant in ur wavelength range)
A_est = spectra_so2'; % only hemoglobin (oxy & deoxy)
% A_est = [spectra_so2', melanin_so2]; % only hemoglobin + melanin
% A_est = [spectra_so2', water_so2, lipid_so2];  % everything but melanin
% A_est = [spectra_so2', water_so2, lipid_so2, melanin_so2];  % everything

oxy_map = zeros(size(Recon_MB_so2,1), size(Recon_MB_so2,2));
deoxy_map = zeros(size(Recon_MB_so2,1), size(Recon_MB_so2,2));
so2_map = zeros(size(Recon_MB_so2,1), size(Recon_MB_so2,2));
res_map = zeros(size(Recon_MB_so2,1), size(Recon_MB_so2,2));
norm_im = zeros(size(Recon_MB_so2,1), size(Recon_MB_so2,2));

%% NNLS fitting of spectrum to loaded spectra (for each pixel)
for i=1:size(Recon_MB_so2,1)
    for j=1:size(Recon_MB_so2,2)
        
        pixSpect = squeeze(Recon_MB_so2(i,j,:));         % single pixel spectrum
        if sum(pixSpect<0) == 0                     % if no negative vals, o.w. return NaN
            norm_im(i,j) = norm(pixSpect);          % calc Euclidean norm of pixel spectrum
            
            if norm(pixSpect)~= 0                   % if non-zero spectrum
                sp_norm = pixSpect/norm_im(i,j);
                [coeffs, resnorm] = lsqnonneg(A_est, sp_norm);
                oxy_map(i, j) = coeffs(2);
                deoxy_map(i, j) = coeffs(1);
                so2_map(i, j) = coeffs(2)./( coeffs(2) + coeffs(1) );
                res_map(i,j) = resnorm;                     % residual
            else
                oxy_map(i, j) = 0;
                deoxy_map(i, j) = 0;
                so2_map(i, j) = 0;
                res_map(i,j) = NaN;
            end
        else
            oxy_map(i, j) = NaN;
            deoxy_map(i, j) = NaN;
            so2_map(i, j) = NaN;
            res_map(i,j) = NaN;
        end
    end
end

%%% plot unmixed spectral maps
% figure; imagesc(oxy_map); axis image off; colorbar; title('oxy');
% figure; imagesc(deoxy_map); axis image off; colorbar; title('deoxy');
figure; imagesc(so2_map); axis image off; colorbar; colormap jet; title('sO2');


%% overlay spectral maps over recon (wl 800nm) (i.e. anatomy)
% VIP: here, colorbars are not meaningful yet!!

anatomy = squeeze(DATA(:,:,wavelengths==800));
mask_oxy = zeros(size(oxy_map));
mask_oxy(oxy_map > 0) =1;
oxy_map_to_overl = ( oxy_map - min(oxy_map(:)) )./max( oxy_map(:) - min(oxy_map(:)) );
% map = [ [0.01:(1/6):1]', zeros(6,1), zeros(6,1)];
% map = [ [0.0, 0.05, 0.1, 0.15, 0.2, 1.0]', zeros(6,1), zeros(6,1)];
colorMap_oxy = autumn(64);
colorMap_oxy = colorMap_oxy(end:-1:1,:);

oxy_ovr = overlay_multipurpose(oxy_map_to_overl, anatomy, anatomy, mask_oxy, colorMap_oxy);
% figure; imagesc(oxy_ovr); axis image off; title('oxy overlayed'), colormap(colorMap_oxy), colorbar;


mask_deoxy = zeros(size(deoxy_map));
mask_deoxy(deoxy_map > 0) =1;
deoxy_map_to_overl = ( deoxy_map - min(deoxy_map(:)) )./max( deoxy_map(:) - min(deoxy_map(:)) );
colorMap_deoxy = winter(64);
colorMap_deoxy = colorMap_deoxy(end:-1:1,:);

deoxy_ovr = overlay_multipurpose(deoxy_map_to_overl, anatomy, anatomy, mask_deoxy, colorMap_deoxy);
% figure; imagesc(deoxy_ovr); axis image off; title('deoxy overlayed'), colormap(colorMap_deoxy),colorbar;


%%% Overlay sO2
mask_so2 = zeros(size(so2_map));
mask_so2(so2_map >= 0) = 1;
so2_map_to_overl = so2_map;
colorMap_sO2 = jet(64);
% show half of the 'colorMap_sO2' colorbar
colorBar_sO2 = colorMap_sO2;
colorBar_sO2 = colorBar_sO2(size(colorBar_sO2,1)/2:end,:);

so2_ovr = overlay_multipurpose(so2_map_to_overl, anatomy, anatomy, mask_so2, colorMap_sO2);
%%% BAD for the DR of the colormap (the ROI is better!!)
figure; imagesc(so2_ovr); axis image off; title('sO2 overlayed'), colormap(colorBar_sO2), colorbar;

%%% try showing only small ROI to have better DR of the colormap
% figure; imagesc(so2_ovr(115:193,47:101)); axis image off; title('sO2 overlayed'), colormap(colorBar_sO2), colorbar;   %% reflection ROI

%%% TRIAL : show only sO2 for an ROI (NOT WORKING TILL NOW!!)
% mask = roipoly;
% masked = so2_map.*mask;
% so2_ovr_masked = overlay_multipurpose(masked, anatomy, anatomy, mask_so2, colorMap_sO2);
% figure; imagesc(so2_ovr_masked); axis image off; title('sO2 overlayed'), colormap(colorBar_sO2), colorbar;


% so2_ovr = overlay(so2_map_to_overl, anatomy, anatomy, mask_so2);
% figure; imagesc(so2_ovr); axis image off; title('sO2 overlayed');


% figure; imagesc(res_map); colormap jet; axis image off; colorbar

%% plot spectrum of single pixel
% figure, imagesc(squeeze(Recon_MB_so2(:,:,end))); colormap gray; axis image off;
% % imagesc(res_map.*rejection_mask); colormap jet; colorbar; axis equal tight;
% 
% while(1)
%     [y, x] = ginput(1);
%     x = round(x); y = round(y);
%     
%     pixSpect = squeeze( Recon_MB_so2( x, y, : ) ); pixSpect = pixSpect./norm(pixSpect,2);
%     %     [ weight_rev_dummy, spec_dummy, measured_s_norm, proj, res ] = get_dummy_fit( sp', M_dummy, Model_dummy, 0 );
%     [coeffs, resnorm] = lsqnonneg(A_est, pixSpect);
%     sp_rec = A_est*coeffs
%     g = figure;
%     plot(wavelengths_so2, [pixSpect, sp_rec], 'LineWidth', 2); legend('orig', 'fitted')
%     % 	plot(wavelengths_so2, [measured_s_norm; spec_dummy], 'LineWidth', 2); legend('orig', 'fitted')
%     title(['Res: ', num2str(res), ', Ref: ', num2str(res_map_dummy(x,y))]);
%     pause;
%     close(g);
%     figure(f);
%     
% end
