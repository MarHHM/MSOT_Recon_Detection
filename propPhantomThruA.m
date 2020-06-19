% propagate phantom thru model matrix
addReconPaths;

%% PARAMS
n = 100;        % im width (pixels)
image_width = 25e-3;        % im width (m)
time_res = 4;                                       % time resolution for model-based (make it bigger for the FWD propagation)
muaTarget = 1;
MB_MAT_FOLDER = 'E:\Marwan\MB_matrices';       % 'E:\Marwan\MB_matrices'  'D:\MB_matrices'

%% Main code
n_angles = 2*n;                                     % number of points for discretizing the curve
T = datainfo.AverageTemperature;
c = 12 + round(1.402385 * 1e3 + 5.038813 * T - 5.799136 * 1e-2 * T^2 + 3.287156 * 1e-4 * T^3 - 1.398845 * 1e-6 * T^4 + 2.787860 * 1e-9 * T^5 );
angle_sensor = datainfo.HWDesc.StartAngle : datainfo.HWDesc.StepAngle : datainfo.HWDesc.EndAngle;
[t, ts] = formInterpolationVec(datainfo, n, time_res, image_width);
% ts = 0 : 1/datainfo.HWDesc.SamplingFrequency : (datainfo.MeasurementDesc.RecordLength-1)/datainfo.HWDesc.SamplingFrequency;     % sampling instants
muaBgd = muaTarget/6;       % add bgd absorption of 1/6 the target absorption (like Luis 2011)

phntmIm = load('E:\Marwan\Dropbox\PA_imaging\MSOT_Recon_Detection\BV2_Derenzo.mat');
phntmIm = phntmIm.BV2;
phntmIm = imresize(phntmIm, n/size(phntmIm,1));
phntmIm = phntmIm + muaBgd;  % add some bgd absorption
for i = 1:size(phntmIm,1)       % truncate to muaTarget
    for j = 1:size(phntmIm,2)
        if phntmIm(i,j)>muaTarget
            phntmIm(i,j) = muaTarget;
        end
    end
end

%figure, imagesc(phntmIm), colormap(bone);

phntmIm = reshape(phntmIm, [n^2 1]);  % reshape im to col major

%%% load model matrix if saved, otherwise build & save it
A_matPath = [MB_MAT_FOLDER '\A_mat_t_res_' num2str(time_res) '_' num2str(n) 'x' num2str(n) '_width_' num2str(image_width*1e3) ...
            '_c_' num2str(c) '_nDet_' num2str(length(angle_sensor)) '_t_' num2str(length(t)) '.mat'];
A_mat = compOrLoadA_mat(A_matPath, c, n, image_width, t, datainfo.HWDesc.Radius, angle_sensor, n_angles);       % make sure that aprams here matches the nes in the A-matPath

%%% propagate
sigMat_numerical= A_mat*phntmIm; 
save(['numerical phantoms\pDetIn_derenzo_timeRes_' num2str(time_res) '.mat'], 'sigMat_numerical', 'n', 'image_width', 'time_res', 'muaTarget', 'muaBgd', 't');


%%% try recon to check
NONNEG = 0;                     % 0 --> normal lsqr
regu = 20;
MB_regu = 1e6*(2^((regu-20)));
% p0convMB = reconstruction(A_mat, pDetIn, n_repIn, 'MB_Tik', MB_regu, NONNEG);
p0convMB = reconstruction(A_mat, sigMat_numerical, n, 'MB_Tik', MB_regu, NONNEG);
figure, imagesc(p0convMB), title('conventional MB Tik recon'), colormap(bone);

% figure, imagesc(phntmIm, [0 1]), colormap('gray');


