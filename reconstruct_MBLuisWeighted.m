function Recon = reconstruct_MBLuisWeighted(A_mat, b_vec, WF_pixels_diag, n, MB_regu, noneg)

        nn = n*n;
        L = MB_regu*calculate_matrix_highpass_1(n);
        % L = sparse(FILT);
        WF_pixels_diag_sparse = MB_regu*sparse(WF_pixels_diag);
        cell_mat{1,1} = A_mat;
        cell_mat{2,1} = L;
        cell_mat{3,1} = WF_pixels_diag_sparse ;     % Luis idea for 'pixels with less coverage' compensation (still to implement)
       
        A_mat = cell2mat(cell_mat);
        b_vec = [b_vec; zeros(nn,1); zeros(nn,1)];
        if noneg
            R=nnls_conjgrad_armijo(A_mat,b_vec,zeros(n*n,1),0.001,5,3);
        else
            R=lsqr_b(A_mat,b_vec,50);
        end
        Recon = reshape(R(:,end),n,n);
        disp('MB_Tik done!')
       




