%% Spectral unmixing (inversion) (linear non-neg fitting)
% VIP --> don't use with negative recon data (data should be reconstructed using 'nnls') so that the spectrum at any pixel doesn't include negative values

%%% load reconStack
full_path = [sigMat_pathName 'recons\reconMB_Tik_imSz_200_slcs_1_wls_8_reps_1.mat'];
reconStack = load(full_path);
reconStack = reconStack.reconStack;
wls = datainfo.Wavelengths;

slicNum = size(reconStack, 4);
wlNum = size(reconStack, 6);
n = size(reconStack, 1);

deoxy = zeros(n,n,slicNum);
oxy = zeros(n,n,slicNum);
SO2 = zeros(n,n,slicNum);
OverlapOxyDeoxy = zeros(n,n,slicNum);
agent_ws = zeros(n,n,slicNum);

for slice_idx=1:slicNum
    
    spectraPath = 'C:\Users\marwan.muhammad\Dropbox\PA_imaging\Unmixing code (linear - Jaber)\SpectralSpecifications1'; % path to oxy- and deoxy-hemoglobin spectra
    spectra = LoadSpectra(spectraPath,  wls);      % loads the spectra (deoxy - oxy) for the selected wavelengths (normalized at 800nm)
    %%%% with agent: spectra(1,:),spectra(2,:) and spectra(3,:) are related to agent, Hb, and Hbo2 respetively.
    %%%% without agent: spectra(2,:) is related to Hbo2 and spectra(1,:) is related to Hb
    %              spectra(3,:) = 1.7e12*(wavelengths).^-3.48;
    
    % aa=ones(1,11);
    % spectra(1,:)=100*spectra(1,:);
    
    % test with nico recon
    
    %             for ind_lam = 1:lambdas
    %                 maximum1(ind_lam) = max(max(Recon_MB(:,:,slice_idx,ind_lam)));
    %             end
    %             maximum2 = max(maximum1);
    %             Recon_MB(:,:,slice_idx,:) = Recon_MB(:,:,slice_idx,:)/maximum2;
    datacube = squeeze(reconStack(:,:, slice_idx,:));
    %             datacube = squeeze(fliplr(Recon_MB(:,:,pos, slice_idx,:)));
    
    
    %%% trim neg values (as mu_a shouldn't be neg) or (better instead to re-recon using 'nnls' inversion)
    for i=1:n
        for j=1:n
            for wl_idx=1:length(wls)
                if datacube(j,i,wl_idx) < 0
                    datacube(j,i,wl_idx) = 0;
                end
            end
        end
    end
    
    %%% unmix deoxy & oxy conc at each pixel
    oxyConcIm = zeros(n,n);       deoxyConcIm = zeros(n,n);
    sO2Im = zeros(n,n);
    SO2Im_raw = zeros(n,n);
    Overlap_Oxy_deoxy = zeros(n,n);
    agent_wsIm = zeros(n,n);
    
    % linearly unmix each pixel
    for i=1:n
        for j=1:n
            [umx, ~] = lsqnonneg(spectra', squeeze(datacube(j,i,:)));
            %     [umx, f0, t] = nnls_conjgrad_armijo(spectra',squeeze(datacube(y,x,:)),zeros(2,1),0.01,100,2);
            %       agent_wsIm(y,x) = umx(1);       deoxyConcIm(y,x) = umx(2);       oxyConcIm(y,x) = umx(3);
            
            deoxyConcIm(j,i) = umx(1);     oxyConcIm(j,i) = umx(2);
        end
    end
    
    figure, imagesc(deoxyConcIm);
    figure, imagesc(oxyConcIm);
    
    
    %% calculate so2 and Overlap of oxy and deoxy without artifact
    
    max_oxy = max(max(oxyConcIm));
    max_deoxy = max(max(deoxyConcIm));
    norm_oxy = oxyConcIm/max_oxy;
    norm_deoxy = deoxyConcIm/max_deoxy;
    
    %%% calc sO2 for each pixel
    for i = 1:n
        for j = 1:n
            Overlap_Oxy_deoxy(i,j) = norm_oxy(i,j).*norm_deoxy(i,j);
            sO2Im(i,j) = ((norm_oxy(i,j)).^2/(norm_deoxy(i,j)+norm_oxy(i,j)));
            if  norm_oxy(i,j)==0
                sO2Im(i,j) = 0;
            end
        end
    end
    
    figure, imagesc(sO2Im), colormap(jet), colorbar;
    
    
    for i = 1:size(oxy,1)
        for j = 1:size(oxy,2)
            if oxyConcIm(i,j)==0 && deoxyConcIm(i,j)==0
                SO2Im_raw(i,j) = 1e-3;
            else
                SO2Im_raw(i,j)=(oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j)));
            end
        end
    end
    
    
    oxy(:,:,slice_idx) = oxyConcIm;
    deoxy(:,:,slice_idx) = deoxyConcIm;
    agent_ws(:,:,slice_idx) = agent_wsIm;
    SO2(:,:,slice_idx) = sO2Im;
    SO2_raw(:,:,slice_idx) = SO2Im_raw;
    OverlapOxyDeoxy(:,:,slice_idx) = Overlap_Oxy_deoxy;
    
    
end
%
%         file_name_deoxy = [datainfo.Name, ' deoxy', '.mat'];
%         file_name_oxy = [datainfo.Name, ' oxy', '.mat'];
%         file_name_agent_ws = [datainfo.Name, ' agent_ws', '.mat'];
%         file_name_SO2 = [datainfo.Name, ' SO2', '.mat'];
%         file_name_SO2_raw = [datainfo.Name, ' SO2_raw', '.mat'];
%         file_name_OverlapOxyDeoxy = [datainfo.Name, ' OverlapOxyDeoxy', '.mat'];
%
%         save_path_deoxy = [current_path, file_name_deoxy];
%         save_path_oxy = [current_path, file_name_oxy];
%         save_path_agent_ws = [current_path, file_name_agent_ws];
%         save_path_SO2 = [current_path, file_name_SO2];
%         save_path_SO2_raw = [current_path, file_name_SO2_raw];
%         save_path_OverlapOxyDeoxy = [current_path, file_name_OverlapOxyDeoxy];
%
%         save(save_path_deoxy, 'deoxy');
%         save(save_path_oxy, 'oxy');
%         save(save_path_agent_ws, 'agent_ws');
%         save(save_path_SO2, 'SO2');
%         save(save_path_SO2_raw, 'SO2_raw');
%         save(save_path_OverlapOxyDeoxy, 'OverlapOxyDeoxy');





%% SO2

%             if (1 >= norm_oxy(i,j) && norm_oxy(i,j) > 0.9) || (1 >= norm_deoxy(i,j) && norm_deoxy(i,j) > 0.9)
% %                SO2Im(i,j)=(oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j))).*0.95;
%                SO2Im(i,j)=(oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j)));
%             elseif (0.9 >= norm_oxy(i,j) && norm_oxy(i,j) > 0.8) || (0.9 >= norm_deoxy(i,j) && norm_deoxy(i,j) > 0.8)
% %                 SO2Im(i,j) = (oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j))).*0.85;
%                 SO2Im(i,j) = (oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j)));
%             elseif (0.8 >= norm_oxy(i,j) && norm_oxy(i,j) > 0.7) || (0.8 >= norm_deoxy(i,j) && norm_deoxy(i,j) > 0.7)
% %                 SO2Im(i,j) = (oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j))).*0.75;
%                 SO2Im(i,j) = (oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j)));
%             elseif (0.7 >= norm_oxy(i,j) && norm_oxy(i,j) > 0.6) || (0.7 >= norm_deoxy(i,j) && norm_deoxy(i,j) > 0.6)
% %                 SO2Im(i,j) = (oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j))).*0.65;
%                 SO2Im(i,j) = (oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j)));
%             elseif (0.6 >= norm_oxy(i,j) && norm_oxy(i,j) > 0.5) || (0.6 >= norm_deoxy(i,j) && norm_deoxy(i,j) > 0.5)
% %                 SO2Im(i,j) = (oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j))).*0.55;
%                 SO2Im(i,j) = (oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j)));
%             elseif (0.5 >= norm_oxy(i,j) && norm_oxy(i,j) > 0.4) || (0.5 >= norm_deoxy(i,j) && norm_deoxy(i,j) > 0.4)
% %                 SO2Im(i,j) = (oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j))).*0.45;
%                 SO2Im(i,j) = (oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j)));
%             elseif (0.4 >= norm_oxy(i,j) && norm_oxy(i,j) > 0.3) || (0.4 >= norm_deoxy(i,j) && norm_deoxy(i,j) > 0.3)
% %                 SO2Im(i,j) = (oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j))).*0.35;
%                 SO2Im(i,j) = (oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j)));
%             elseif (0.3 >= norm_oxy(i,j) && norm_oxy(i,j) > 0.2) || (0.3 >= norm_deoxy(i,j) && norm_deoxy(i,j) > 0.2)
% %                 SO2Im(i,j) = (oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j))).*0.25;
%                 SO2Im(i,j) = (oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j)));
%             elseif (0.2 >= norm_oxy(i,j) && norm_oxy(i,j) > 0.1) || (0.2 >= norm_deoxy(i,j) && norm_deoxy(i,j) > 0.1)
% %                 SO2Im(i,j) = (oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j))).*0.15;
%                 SO2Im(i,j) = (oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j)));
%             elseif (0.1 >= norm_oxy(i,j) && norm_oxy(i,j) > 0) || (0.1 >= norm_deoxy(i,j) && norm_deoxy(i,j) > 0.0)
%                 SO2Im(i,j) = (oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j))).*0.05;
% %                 SO2Im(i,j) = (oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j)));
%             else
%                 SO2Im(i,j) = 0;
%             end


%             if (0.2 >= norm_oxy(i,j)) && (0.2 >= norm_deoxy(i,j))
%                 SO2Im(i,j) = (oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j)));
%             elseif (0.15 >= norm_oxy(i,j)) && (0.15 >= norm_deoxy(i,j))
%                 SO2Im(i,j) = (oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j))).*0.15;
%             elseif (0.1 >= norm_oxy(i,j)) && (0.1 >= norm_deoxy(i,j))
%                 SO2Im(i,j) = (oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j))).*0.1;
%             elseif (0.05 >= norm_oxy(i,j)) && (0.05 >= norm_deoxy(i,j))
%                 SO2Im(i,j) = (oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j))).*0.05;
%             elseif  norm_oxy(i,j)==0
%                 SO2Im(i,j) = 0;
%             else
%                 SO2Im(i,j) = (oxyConcIm(i,j)/(deoxyConcIm(i,j)+oxyConcIm(i,j)));
%             end










