% function [ mu_p0, sigma_p0, reconStack_p0ColMajor] = form_p0_priors(reconStack)
% calculates mean vector & cov matrix of recons reps
%   Detailed explanation goes here

%% PARAMS
% mu_p0 = 0;
% sigma_p0 = 0;
Zpos_begin = 1;
Zpos_end = size(reconStack, 4);
wl_begin = 1;
wl_end = size(reconStack, 3);
rep_begin = 1;
rep_end = size(reconStack, 5);
im_ht = size(reconStack,1);
im_wd = size(reconStack,2);
im_sz = im_ht*im_wd;

%% convert recon images to column-major vecs 
% & rearrange reconstack dims to a logical order (descending --> Zpos, wl,
% rep, im) & normalize each rep (Q: is normalization reqd?)
reconStack_p0ColMajor = zeros(Zpos_end-Zpos_begin+1, wl_end-wl_begin+1, rep_end-rep_begin+1, im_sz);
for slc_idx = 1 : Zpos_end-Zpos_begin+1
    for wl_idx = 1 : wl_end-wl_begin+1
        for rep_idx = 1 : rep_end-rep_begin+1
            reconStack_p0ColMajor(slc_idx,wl_idx,rep_idx,:) = reshape(reconStack(:,:,wl_idx,slc_idx,rep_idx), im_sz, 1);
%             reconStack_p0ColMajor(slc_idx,wl_idx,rep_idx,:) = mat2gray(reconStack_p0ColMajor(slc_idx,wl_idx,rep_idx,:));
        end
    end
end

%% form priors (per wl per Zpos)
p0mean = zeros(Zpos_end-Zpos_begin+1, wl_end-wl_begin+1, im_sz);

% take the mean over the "reps" dim
for slc_idx = 1 : Zpos_end-Zpos_begin+1
    for wl_idx = 1 : wl_end-wl_begin+1
        p0mean(slc_idx,wl_idx,:) = mean(reconStack_p0ColMajor(slc_idx, wl_idx, :, :), 3);
    end
end

% form the Cov matrix (per Zpos, per wl)
p0cov = zeros(Zpos_end-Zpos_begin+1, wl_end-wl_begin+1, im_sz, im_sz);
% p0covInv = zeros(Zpos_end-Zpos_begin+1, wl_end-wl_begin+1, im_sz, im_sz);
% p0covInvCholsky = zeros(Zpos_end-Zpos_begin+1, wl_end-wl_begin+1, im_sz, im_sz);
for slc_idx = 1 : Zpos_end-Zpos_begin+1
    for wl_idx = 1 : wl_end-wl_begin+1
        tic;
        p0cov(slc_idx,wl_idx,:,:) = cov(squeeze(reconStack_p0ColMajor(slc_idx, wl_idx, :, :)));
        toc;
    end
end

% calc cholesky decomp & its inv of every cov mat (per slc per wl)
for slc_idx = 1 : Zpos_end-Zpos_begin+1
    for wl_idx = 1 : wl_end-wl_begin+1
        tic;
        temp = chol(squeeze(p0cov(slc_idx,wl_idx,:,:)), 'lower');
        p0covInvCholsky(slc_idx,wl_idx,:,:) = inv(temp);
        toc;
    end
end
% end

