function [ sigMat ] = filter_function( sigMat, filter_f,fs)
%FILTER_FUNCTION Summary of this function goes here
%   Detailed explanation goes here
% filter parameters
% disp('-> Filtering 1D channels data before reconstruction...');
if( filter_f(2) ~= 0 )
    f_LPF = filter_f(2);
    [b_LPF,a_LPF] = cheby1( 8, .01, 2 * f_LPF/fs * .92 );
end
if( filter_f(1) ~= 0 )
    f_HPF = filter_f(1);
    if( 2 * f_HPF/fs < 0.016 )
        [b_HPF,a_HPF] = cheby1( 2, .01, 2 * f_HPF/fs * 3.3, 'high' );
    else
        [b_HPF,a_HPF] = cheby1( 4, .01, 2 * f_HPF/fs * 1.46, 'high' );
    end
end

% filtering
if( filter_f(2) ~= 0 )
    sigMat = filtfilt(b_LPF, a_LPF, sigMat);
end
if( filter_f(1) ~= 0 )
    sigMat = filtfilt(b_HPF, a_HPF, sigMat);
end

end

