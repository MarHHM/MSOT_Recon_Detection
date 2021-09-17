function [ Recon ] = reconstruction(A_mat, b_vec,n,method,regu,noneg)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

switch method
    case 'MB_Tik'
        nn = n*n;
        FILT = regu*calculate_matrix_highpass_1(n);
        L = sparse(FILT);
        cell_mat{1,1} = A_mat;
        cell_mat{2,1} = L;
        A_mat = cell2mat(cell_mat);
        b_vec = [b_vec; zeros(nn,1)];
        if noneg
            R=nnls_conjgrad_armijo(A_mat,b_vec,zeros(n*n,1),0.001,5,3);
        else
            R=lsqr_b(A_mat,b_vec,50);
        end
        Recon=reshape(R(:,end),n,n);
    case 'MB_TVL1'
        TVWeight = 0e5; 	% Weight for TV penalty
        xfmWeight = regu;	% Weight for Transform L1 penalty
        Itnlim = 50;		% Number of iterations
        Recon  = TVL1_yiyong(A_mat,b_vec,TVWeight,xfmWeight,Itnlim,noneg);
end
end