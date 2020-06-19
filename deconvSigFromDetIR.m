struct = load('IR_MSOT256_microsphere100um_based_on_integration.mat');
IR = struct.IR;
FR = struct.FR;

sigBefore = sigMat_truncated(:, 1, 1, 1, 1, 1);

% [sigAfter, remainder] = deconv(sigBefore, IR);

Y = fft(sigBefore);
figure, plot(abs(fftshift(Y)));
figure, plot(abs(FR));
H = Y./FR;
figure, plot(abs(H));
% sigAfter = flipud(abs(ifft(H(:))));
sigAfter = abs(ifft(H(:)));

figure, subplot(2,1, 1), plot(sigBefore);
        subplot(2,1, 2), plot(sigAfter);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test
close all;
Fs = datainfo.HWDesc.SamplingFrequency;
fc = 3.5e6;     % close to central freq of the detector
tVec = 0:1/Fs:(length(sigBefore)-1)/Fs;
mySig = sin(2*pi*fc*tVec);
figure, plot(mySig);

convolved = conv(mySig, IR, 'same');
figure, plot(convolved);

deconvolved = deconv(convolved, IR);
figure, plot(deconvolved);