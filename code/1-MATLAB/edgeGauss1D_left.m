function arrEG = edgeGauss1D_left(arr, mm, s)

Filt = ones(size(arr));

GG = gauss_distribution(-mm:mm, 0, s); 
GG = GG / max(GG(:));
Filt(1:mm) = GG(2:mm+1); 

arrEG = arr.*Filt;
end
