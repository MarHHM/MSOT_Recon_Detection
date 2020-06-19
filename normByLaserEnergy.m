function sigMat = normByLaserEnergy(sigMat, datainfo)

laserEnergy_vec = [datainfo.ScanFrames.LaserEnergy];
laserEnergy_vec = laserEnergy_vec(1:length(datainfo.Wavelengths));

for run_idx = 1:1
    for slc_idx = 1 : size(sigMat, 4)
        for wl_idx = 1 : size(sigMat, 6)
            for rep_idx = 1 : size(sigMat, 5)
                sigMat(:, :, run_idx, slc_idx, rep_idx, wl_idx) = ...
                    sigMat(:, :, run_idx, slc_idx, rep_idx, wl_idx)/laserEnergy_vec(wl_idx);
            end
        end
    end
end