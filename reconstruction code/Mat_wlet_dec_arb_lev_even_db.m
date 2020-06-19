function [ A_out ] = Mat_wlet_dec_arb_lev_even_db( A, n_t,n_proj,depth,wl_name ) %Performs time-proj 2d-wlet directly on the matrix
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%temp=calc_db2_wlet_dec2_pack_arb_lev_even_db(full(reshape(A(:,1),n_t,n_proj)),depth,wl_name).';
%A_out=(zeros(length(temp),size(A,2)));

%if you want to conserve memory use A_out=sparse(zeros(length(temp),size(A,2))); 
%this option is slower though

for kk=1:size(A,2)  %this calculates the multipication of C_doulbe with A_mat, more accurate than the matrix
  A_out(:,kk)=calc_db2_wlet_dec2_pack_arb_lev_even_db(full(reshape(A(:,kk),n_t,n_proj)),depth,wl_name).';
%   if mod(kk,100)==1
%       kk
%   end
end
A_out=sparse(A_out);

end

