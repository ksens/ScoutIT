function arrEG = edgeGauss1D_right(arr, mm, s)

Filt = ones(size(arr));

GG = gauss_distribution(-mm:mm, 0, s); 
GG = GG / max(GG(:));
Filt(end-mm:end) = GG(mm+1:end);

arrEG = arr.*Filt;
end
