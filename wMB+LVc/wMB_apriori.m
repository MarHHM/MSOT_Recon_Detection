% main idea: use expression for 'P_r(t_ij)' when acoustic apriori are known (eq. 7 in APL 2011) in the expression for 
% 'P_d(t_ij)' for the MB case (eq.9 in PMB 2013)


%% This part is from APL 2011 (wBP+apriori)
%%% estimate 'P_r(t_ij)' by Monte Carlo --> probability that the wave detected at (transducer i, time j) is a reflected or scattered wave wit unit amplitude 
% (eq 7. in APL 2011)
P_r.hist = calculate_time_histograms_refs( logicalMask_A, logicalMask_B, t, P_r.nPairs, r_sensor, angle_sensor, c, P_r.nHist, im_w, n );

%%% calc 'P_d(t_ij)' --> probability that the wave detected at (transducer i, time j) is a direct wave (from a point on the circumference del(A_ij,int)) --
% eq 2. in APL 2011
P_d = zeros( 1, nDet*sizeT );
WT = 1;
for i = 1:nDet          % for each transducer
  switch WT
    case 1
      P_d((i-1)*sizeT+1 : i*sizeT) = max( 0, 1-w*(P_r.hist(:,i)/max( P_r.hist(:,i) )) );    % by substituing eq. 8 in eq. 2
      
    case 2
      P_d((i-1)*sizeT+1 : i*sizeT) = max( 0, 1-w*(P_r.hist(:,i)/max( P_r.hist(:,i) )).^2 );
      
  end
end


%% This part is from PMB 2013 (wMB)  
% put 'P_d' in matrix form for use with MB (eq.9 in PMB 2013)
pos = linspace( 1, nDet*sizeT, nDet*sizeT );
W = sparse( pos, pos, P_d, nDet*sizeT, nDet*sizeT );
% weighting A_mat & b_vec
A_matW = W*A_mat;
b_vec_In = prepareMeasuredPressureVec( sigMat_current, ts, t, nDet );       % interploate to lower res & reshape to col-major
b_vecW = W*b_vec_In;

%% do recon
ReconW(:, :, 1, slc_idx, rep_idx, wl_idx) = reconstruction( A_matW, b_vecW, n, RECON_METHOD, MB_regu, IMPOSE_NONNEGATIVITY );

