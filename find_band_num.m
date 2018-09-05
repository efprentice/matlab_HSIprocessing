function [band_num] = find_band_num(wavelength)
%FIND_BAND_NUM Finds the nearest band number to a desired wavelength;
% wavelength input must be one-dimensional, any length
% HSI_wlens (wlens variable) should be updated for every new camera

load('HSI_wlens.mat')

band_num = zeros(1,length(wavelength));
for i = 1:length(wavelength)
    [~,bn] = min(abs(wlens-wavelength(i)));
    band_num(i) = bn;
end

end

