function [umx A flag message] = unmix(mixed, w, spectra, method, K, covAll,Thresh, n1, n2)

message = [];
flag = 1;
nm_sp = size(spectra,1);
A_est = spectra';
%Normalize Spectra
for i=1:size(A_est,2)
    A_est(:,i) = A_est(:,i)./norm(A_est(:,i),2);
end
S = A_est(:,1);

switch method
    case 'OSP'
        tmp = squeeze(mixed);
        nm_sp = size(spectra,1);
        
        if nm_sp<3
            message = 'No hemoglobin spectra loaded';
            flag = 0;
            return;            
        end
        %         spectra(nm_sp+1,:) = ones(size(w));
        A_est = spectra';     
        pixels = tmp(:,:);

        %Normalize Spectra
        for i=1:size(A_est,2)
           A_est(:,i) = A_est(:,i)./norm(A_est(:,i));  
        end
      
        p = size(A_est,2);
        for i=1:1
           d = A_est(:,i);
           U = [A_est(:,1:i-1)  A_est(:,i+1:end) ];
           P = (eye([length(d) length(d)]) - U*pinv(U));
           operator = d'*P;
           unmixed = operator*pixels ;
           umx(i,1,:,:) = reshape(unmixed,n1,n2)./(d'*P*d); 
        end
        A = A_est;          
        umx(umx<0) = 0;
        umx = squeeze(umx);
        flag = 1;
    case 'Simple'  % in case of contaminated background or simple diagonal loading
        
        tmp = squeeze(mixed);
        x = tmp(:,:);
        [x, mixedmean] = remmean(x);
        psi1 = n1*n2; psi2 = 0;
        
        G = cov(x');
        G_inv = inv(G);
        
        %Adaptive Matched Filter
        unmixed = (S'*G_inv*x);
        unmixed(unmixed<0) = 0;
        den = psi1*(S'*G_inv*S);
        unmixed = (unmixed.^2)./den;

        umx = reshape(unmixed,n1,n2);
        umx(umx<Thresh) = 0;
        p = size(A_est,2);
        A = S;
        flag = 1;
    case 'EC_GLRT'
        
        tmp = squeeze(mixed);
        x = tmp(:,:);
        [x, mixedmean] = remmean(x);
        psi1 = n1*n2; psi2 = 0;

        addpath('D:\Users\stratis.tzoumas\Documents\MATLAB\Robust_AMF_27_10_2015\StatisticalCharacterization\');
        v = find_t_DOF(x);
        G = cov(x');
        
        G_inv = inv(G);
        %EC GLRT
        unmixed = zeros(1,size(x,2));
        for j=1:size(x,2)
            unmixed(j) = sqrt((v-1)/((v-2)+x(:,j)'*G_inv*x(:,j)))*(S'*G_inv*x(:,j))/sqrt((S'*G_inv*S));
        end 
        unmixed(unmixed<0) = 0;
        umx = reshape(unmixed,n1,n2);
        umx = umx.^2;
        p = size(A_est,2);
        A = S;
     case 'QL_shrinkage_adaptive'
        
        tmp = squeeze(mixed);
        x = tmp(:,:);
        [x, mixedmean] = remmean(x);
        psi1 = n1*n2; psi2 = 0;

        size(covAll)
        %Quasi local estimator of covariance matrix
        G_init = cov(x');
        [V D] = eig(K); %K=V*D*V'
        C_QL = V*diag(diag(V'*G_init*V))*V';
        
        G_init = G_init./norm(G_init,'fro');
        C_QL = C_QL./norm(C_QL,'fro');
  
        for kk=1:length(covAll)
            G_dist_All(kk) = norm(G_init./norm(G_init,'fro')-covAll{kk}./norm(covAll{kk},'fro'),'fro');
        end
%                 a1 =  sort(tmp1);
        a2 =  sort(G_dist_All);
%         figure,plot(sort(G_dist_All));
%                 metric1(mouse_idx,i,j,ii) = a1(10);
        metric = a2(1);
        metric

        betha = 4*sqrt(metric); if betha>1 betha=1; end
        betha
        
        G_init = cov(x');
        [V D] = eig(K); %K=V*D*V'
        C_QL = V*diag(diag(V'*G_init*V))*V';
        G = (1-betha)*G_init+betha*C_QL;
        G_inv = inv(G);    
        %Adaptive Matched Filter
        unmixed = (S'*G_inv*x);
        unmixed(unmixed<0) = 0;
        den = psi1*(S'*G_inv*S);
        unmixed = (unmixed.^2)./den;
        umx = reshape(unmixed,n1,n2);
        umx(umx<Thresh) = 0;
%         umx = sqrt(umx);
        p = size(A_est,2);
        A = betha;
     case 'RSDF'
         
        tmp = squeeze(mixed);
        x = tmp(:,:);
        [x, mixedmean] = remmean(x);
        psi1 = n1*n2; psi2 = 0;

        addpath('D:\Users\stratis.tzoumas\Documents\MATLAB\Robust_AMF_27_10_2015\Covariance_metrics\');
        addpath('D:\Users\stratis.tzoumas\Documents\MATLAB\Robust_AMF_27_10_2015\StatisticalCharacterization\');       
     
        %Quasi local estimator of covariance matrix
        G_init = cov(x');
        [V D] = eig(K); %K=V*D*V'
        C_QL = V*diag(diag(V'*G_init*V))*V';
        
        G_init = G_init./norm(G_init,'fro');
        C_QL = C_QL./norm(C_QL,'fro');
               
        for kk=1:length(covAll)
            G_dist_All(kk) = norm(G_init./norm(G_init,'fro')-covAll{kk}./norm(covAll{kk},'fro'),'fro');
        end
        a =  sort(G_dist_All);
        metric = a(1);
        metric

        betha = 4*sqrt(metric); if betha>1 betha=1; end
        betha
        %Quasi local estimator of covariance matrix
        G_init = cov(x');
        [V D] = eig(K); %K=V*D*V'
        C_QL = V*diag(diag(V'*G_init*V))*V';
        
        G = (1-betha)*G_init+betha*C_QL;
        G_inv = inv(G);    
        
        unmixed = (S'*G_inv*x);
        unmixed(unmixed<0) = 0;
        den = psi1*(S'*G_inv*S);
        unmixed = (unmixed.^2)./den;
        tt = sort(unmixed);
        
        test = mean(tt(end-10:end));
        
        if test<4e-3 & betha<1
            v = find_t_DOF(x);
            v
            unmixed = zeros(1,size(x,2));
            for j=1:size(x,2)
                unmixed(j) = sqrt((v-1)/((v-2)+x(:,j)'*G_inv*x(:,j)))*(S'*G_inv*x(:,j))/sqrt((S'*G_inv*S));
            end 
            unmixed(unmixed<0) = 0;
            umx = reshape(unmixed,n1,n2);
            umx = umx.^2;
            umx(umx<Thresh) = 0;
            umx = umx.^2;
            umx = umx.^2;

        else
            unmixed = (S'*G_inv*x);
            unmixed(unmixed<0) = 0;
            den = psi1*(S'*G_inv*S);
            unmixed = (unmixed.^2)./den;
%             unmixed = sqrt(unmixed);
%             unmixed = sqrt(unmixed);
            umx = reshape(unmixed,n1,n2);
            umx(umx<Thresh) = 0;
        end

        p = size(A_est,2);
        A = S;
    case 'QL_shrinkage_adaptive_interp'
        
        tmp = squeeze(mixed);
        x = tmp(:,:);
        [x, mixedmean] = remmean(x);
        psi1 = n1*n2; psi2 = 0;
         
        %% based on the current  wavelength sampling, recover the index of the prior
        if min(w)<700 || max(w)>900
            disp('No prior modeling available for current wavelengths');
            return;
        end
        w_prior = 700:10:900;
        for w_idx=1:length(w)
           [dummy idx_sampling(w_idx)] = min(abs(w_prior - w(w_idx))) ;
        end
     
        %Quasi local estimator of covariance matrix
        G_init = cov(x');
        [V D] = eig(K(idx_sampling,idx_sampling)); %K=V*D*V'
        C_QL = V*diag(diag(V'*G_init*V))*V';
        
        G_init = G_init./norm(G_init,'fro');
        C_QL = C_QL./norm(C_QL,'fro');
               
        for kk=1:length(covAll)
            G_prior_i = covAll{kk}; G_prior_i_sampled = G_prior_i(idx_sampling,idx_sampling);
            G_prior_i_sampled = G_prior_i_sampled./norm(G_prior_i_sampled,'fro');
            G_dist_All(kk) = norm(G_init./norm(G_init,'fro')-G_prior_i_sampled,'fro');
        end
        
        a2 =  sort(G_dist_All);
        metric = a2(1);
        metric

        betha = 4*sqrt(metric); if betha>1 betha=1; end
        betha
         
        %Quasi local estimator of covariance matrix
        G_init = cov(x');
        [V D] = eig(K(idx_sampling,idx_sampling)); %K=V*D*V'
        C_QL = V*diag(diag(V'*G_init*V))*V';
        
        G = (1-betha)*G_init+betha*C_QL;
        G_inv = inv(G);    

        unmixed = (S'*G_inv*x);
        unmixed(unmixed<0) = 0;
        den = psi1*(S'*G_inv*S);
        unmixed = (unmixed.^2)./den;
        umx = reshape(unmixed,n1,n2);
        umx(umx<Thresh) = 0;
%         umx = sqrt(umx);
        p = size(A_est,2);
        A = betha;
    
    case 'RSDF_interp'
         
        tmp = squeeze(mixed);
        x = tmp(:,:);
        [x, mixedmean] = remmean(x);
        psi1 = n1*n2; psi2 = 0;
        
        %% based on the current  wavelength sampling, recover the index of the prior
        if min(w)<700 || max(w)>900
            disp('No prior modeling available for current wavelengths');
            return;
        end
        w_prior = 700:10:900;
        for w_idx=1:length(w)
           [dummy idx_sampling(w_idx)] = min(abs(w_prior - w(w_idx))) ;
        end
     
        %Quasi local estimator of covariance matrix
        G_init = cov(x');
        [V D] = eig(K(idx_sampling,idx_sampling)); %K=V*D*V'
        C_QL = V*diag(diag(V'*G_init*V))*V';
        
        G_init = G_init./norm(G_init,'fro');
        C_QL = C_QL./norm(C_QL,'fro');
               
        for kk=1:length(covAll)
            G_prior_i = covAll{kk}; G_prior_i_sampled = G_prior_i(idx_sampling,idx_sampling);
            G_prior_i_sampled = G_prior_i_sampled./norm(G_prior_i_sampled,'fro');
            G_dist_All(kk) = norm(G_init./norm(G_init,'fro')-G_prior_i_sampled,'fro');
        end
        a =  sort(G_dist_All);
        metric = a(1);
        metric

        betha = 4*sqrt(metric); if betha>1 betha=1; end
 
        
        %Quasi local estimator of covariance matrix
        G_init = cov(x');
        [V D] = eig(K(idx_sampling,idx_sampling)); %K=V*D*V'
        C_QL = V*diag(diag(V'*G_init*V))*V';
        
        G = (1-betha)*G_init+betha*C_QL;
        G_inv = inv(G);    

        unmixed = (S'*G_inv*x);
        unmixed(unmixed<0) = 0;
        den = psi1*(S'*G_inv*S);
        unmixed = (unmixed.^2)./den;
        tt = sort(unmixed);
        
        test = mean(tt(end-10:end));
        
        if test<4e-3 & betha<1
            v = find_t_DOF(x);
            v
            unmixed = zeros(1,size(x,2));
            for j=1:size(x,2)
                unmixed(j) = sqrt((v-1)/((v-2)+x(:,j)'*G_inv*x(:,j)))*(S'*G_inv*x(:,j))/sqrt((S'*G_inv*S));
            end 
            unmixed(unmixed<0) = 0;
            umx = reshape(unmixed,n1,n2);
            umx = umx.^2;
            umx(umx<Thresh) = 0;
            umx = umx.^2;
            umx = umx.^2;

        else
            unmixed = (S'*G_inv*x);
            unmixed(unmixed<0) = 0;
            den = psi1*(S'*G_inv*S);
            unmixed = (unmixed.^2)./den;
%             unmixed = sqrt(unmixed);
%             unmixed = sqrt(unmixed);
            umx = reshape(unmixed,n1,n2);
            umx(umx<Thresh) = 0;
        end

        p = size(A_est,2);
        A = S;
end

end