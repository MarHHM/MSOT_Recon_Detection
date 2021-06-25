function [t, ts] = formInterpolationVec(datainfo, n, time_res, image_width, c)
% ts --> sampling instants
% t --> interpolated (& cut) time instants (used usually with MB recons)

% forming truncated interpolation t vec as in code for reconing corresp images (taken AS IS from recon code)
RecordLength = double(datainfo.MeasurementDesc.RecordLength);
ts = 0 : 1/datainfo.HWDesc.SamplingFrequency : (RecordLength-1)/datainfo.HWDesc.SamplingFrequency; % sampling instants

limits(1) = datainfo.HWDesc.Radius-(image_width)*sqrt(2)/2;       % limits for the signal
limits(2) = datainfo.HWDesc.Radius+(image_width)*sqrt(2)/2;       % limits for the signal
dx = image_width/n;                                 % increment in x
dt = dx/(time_res*c);                              % increment in t employed to make the model-based reconstruction
fac = datainfo.HWDesc.SamplingFrequency/c;
pos_start = max(1,int32((limits(1))*fac));
pos_end = min(RecordLength,int32((limits(2))*fac));
t = ts(pos_start):dt:ts(pos_end);           % downsampled (& cut) time vector  (less than ts)