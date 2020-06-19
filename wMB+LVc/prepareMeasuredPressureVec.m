function b_vec = prepareMeasuredPressureVec( sigMat_current, ts, t, NumDetectors )

% interpolate raw data to lower resolution (downsampling) & truncate
sigMat2 = zeros(length(t), size(sigMat_current, 2));
for j = 1:NumDetectors
  sigMat2(:,j) = interp1(ts, sigMat_current(:,j), t);
end

% reshape interpolated measures to column-major format
b_vec = reshape( sigMat2, length(t)*NumDetectors, 1);