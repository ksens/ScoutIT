function arrEG = edgeGauss1D_v2(arr, mm, s)
% edgeGauss1d apply edge gauss on 1D signal
%   arr -> array on which to apply edge-gauss
%   mm -> half-width of Gaussian (depth till which array is affected)
%   s -> standard deviation of gaussian

arr = edgeGauss1D_left(arr, mm, s);
arrEG = edgeGauss1D_right(arr, mm, s);

% plot(arr), hold on, plot(arr.*Filt, 'r'), hold off
end

