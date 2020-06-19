% function [ mu_p0, sigma_p0, reconStack_p0ColMajor] = form_p0_priors(reconStack)
% calculates mean vector & cov matrix of recons reps
%   Detailed explanation goes here

%% PARAMS
% mu_p0 = 0;
% sigma_p0 = 0;
% load any recon just to read the dims from it
reconStack_currSOS = load('D:\Marwan\Datasets_offline\Marwan_1\Scan_6\last recons 20 sos 100 n MB\reconMB_Tik_imSz100_sos1480_rep1to500');
reconStack_currSOS = reconStack_currSOS.reconStack;
Zpos_begin = 1;
Zpos_end = size(reconStack_currSOS, 4);
slcNum = Zpos_end-Zpos_begin+1;
wl_begin = 1;
wl_end = size(reconStack_currSOS, 3);
wlNum = wl_end-wl_begin+1;
rep_begin = 1;
rep_end = size(reconStack_currSOS, 5);
repNum = rep_end-rep_begin+1;
im_ht = size(reconStack_currSOS,1);
im_wd = size(reconStack_currSOS,2);
im_sz = im_ht*im_wd;

%% concatenate reps for each pixel from diffenret SOS recons (to increase the effective #obs to be = #RVs for the Cov Mat to be pos def)
sos_vec = linspace(1480,1540,20);
sos_vec = [sos_vec 1543.1579];
reconStack_imColMajor = zeros(slcNum, wlNum, repNum*length(sos_vec), im_sz);
pxlVsRep_currSOS = zeros(repNum,im_sz);
for sos_idx = 1:length(sos_vec)
    sos_current = sos_vec(sos_idx);
    reconStack_currSOS = load(['D:\Marwan\Datasets_offline\Marwan_1\Scan_6\last recons 20 sos 100 n MB\reconMB_Tik_imSz100_sos' num2str(sos_current) '_rep1to500.mat' ]);
    reconStack_currSOS = reconStack_currSOS.reconStack;
    
    for slc_idx = 1 : slcNum
        for wl_idx = 1 : wlNum
            for rep_idx = 1 : repNum
                % convert recon images to column-major vecs
                pxlVsRep_currSOS(rep_idx,:) = reshape(reconStack_currSOS(:,:,wl_idx,slc_idx,rep_idx), im_sz, 1);
            end
            % rearrange reconstack dims to a logical order (descending --> Zpos, wl, rep, im)
            reconStack_imColMajor(slc_idx,wl_idx,((sos_idx-1)*rep_idx)+1:sos_idx*rep_idx,:) = pxlVsRep_currSOS;
        end
    end
end

%% form priors (per wl per Zpos)
p0mean = zeros(slcNum, wlNum, im_sz);

% take the mean over the "reps" dim
for slc_idx = 1 : slcNum
    for wl_idx = 1 : wlNum
        p0mean(slc_idx,wl_idx,:) = mean(reconStack_imColMajor(slc_idx, wl_idx, :, :), 3);
    end
end

% form the Cov matrix (per Zpos, per wl)
p0cov = zeros(slcNum, wlNum, im_sz, im_sz);
% p0covInv = zeros(slcNum, wlNum, im_sz, im_sz);
% p0covInvCholsky = zeros(slcNum, wlNum, im_sz, im_sz);
for slc_idx = 1 : slcNum
    for wl_idx = 1 : wlNum
        tic;
        p0cov(slc_idx,wl_idx,:,:) = cov(squeeze(reconStack_imColMajor(slc_idx, wl_idx, :, :)));
        toc;
    end
end

%% calc cholesky decomp & its inv of every cov mat (per slc per wl)
for slc_idx = 1 : slcNum
    for wl_idx = 1 : wlNum
        tic;
        temp = chol(squeeze(p0cov(slc_idx,wl_idx,:,:)), 'lower');
        p0covInvCholsky(slc_idx,wl_idx,:,:) = inv(temp);
        toc;
    end
end
% end

