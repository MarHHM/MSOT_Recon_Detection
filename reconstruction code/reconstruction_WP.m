function [ Recon_sum ] = reconstruction_WP(A_mat,Ainv_A, b_vec,places,BB_part_temp,n,t,n_proj,noneg)

wl_name='db6';
depth_proj=2;
depth_im=2;
sizeT = length(t);
PPa_wlet= Mat_wlet_dec_arb_lev_even_db( b_vec, sizeT,n_proj,depth_proj,wl_name);        
recon=zeros(n*n,4^(depth_im)); 
Recon_sum=0;

 for kk=1:4^(depth_im)
%kk
     PPa_wlet_part=PPa_wlet(places{kk}).';
     recon_wlet=Ainv_A{kk}*sparse(PPa_wlet_part(:));
     recon_wlet_cell{kk}=recon_wlet(:);
     recon(:,kk)=BB_part_temp{kk}*recon_wlet(:);
      Recon=reshape(recon(:,kk),n,n);
     Recon_sum=Recon_sum+Recon;
     
 end
 
if noneg
Recon_sum=max(Recon_sum,0);
end

% for rr=2:10
% %rr
%     Recon_sum_old=Recon_sum;
%     Recon_sum=0;
%     sig_low_time=A_mat*sum(recon,2);
%     sig_low_temp= Mat_wlet_dec_arb_lev_even_db(sig_low_time,sizeT,n_proj,depth_proj,wl_name);
%     
%     for kk=1:4^(depth_im) 
%     sig_low=sig_low_temp(places{kk});
%     PPa_wlet_part=PPa_wlet(places{kk}).'-1*sig_low(:).';
%     recon_wlet=Ainv_A{kk}*sparse(PPa_wlet_part(:));
%     recon_wlet_cell{kk}=recon_wlet_cell{kk}+0.2*recon_wlet(:);
%     recon(:,kk)=BB_part_temp{kk}*recon_wlet_cell{kk}(:);
%     Recon=reshape(recon(:,kk),n,n);
%     Recon_sum=Recon_sum+Recon;
%     
%     
%     end 
%     if noneg
% Recon_sum=max(Recon_sum,0);
%     end
% 
% end 

 disp('WaveletPacket done!')





