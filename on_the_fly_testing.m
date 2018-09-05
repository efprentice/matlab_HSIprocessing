% HSI v4 Data Processing: Trondheimsfjorden ground based measurements
% E.F. Prentice
% 16 July 2018
% yuhuuuu

% Description: To view raw sample spectral intensity plots and a spatial
% image of the HSI scan, edit .png folder location (variable 'path_name').
% functions 'spec_lim_matching' and 'find_band_num' required

% Research questions: Can I determine if the images are over- / under-
% exposed?  What is the quality of the spatial scan?  Can I discern objects
% and landmarks? What parts of the image are over/under exposed?

% Results: Raw spectral intensity plots quickly show over/under exposure.
% Yes, with a leveled tripod, image quality is ok. Auto-rotate device 
% coordinated with fps would be better. Spatial image helps with visibly 
% identifying exposure limitations, objects.

% Next steps: Select only the spectra of water pixels to evaluate exposure.
% Show min/max exposure plots, not just 6 evenly spaced samples. Trim left,
% right edges of scans. Make RGB band selection more robust. Characterize 
% ocean pixel spectra and see if algae pixels look different. Write a
% script to create wlens.mat files easily. Speed up code.

close all; clear all; clc

% Assumptions ( + running required input functions)

load('HSI_wlens.mat')                    % scaled wavelength values, called wlens (x-axis vals) THIS SHOULD BE DONE FOR EACH CAMERA
band_num = find_band_num([638,543,460]); % band numbers of a desired wavelengths, I like RGB [638,543,460]
spec_max = 40950;                        % maximum intensity count, from 10*(-1+2^12)


% Load HSI .png images (specify folder location)

%path_name = uigetdir('.','Select HSI folder');
path_name = 'C:\Users\elizabep\Documents\fieldwork_data\180724_Trondfj-T2\smiskaret\scan2';
%path_name = 'C:\Users\elizabep\Desktop\vmshare\test';

file_names = dir(strcat(path_name,[filesep '*.png']));
n_frames = size(file_names,1);
samp_frames = floor(linspace(1,n_frames,6));
test_frame = imread([file_names(samp_frames(1)).folder, filesep, file_names(samp_frames(1)).name]);


% Plot spectral intensity of sample frames

figure('Name','Spectral Signatures of Raw Sample Frames')
for idx = 1:length(samp_frames)
    frame = imread([file_names(samp_frames(idx)).folder, filesep, file_names(samp_frames(idx)).name]);
    subplot(2,3,idx)
    plot(wlens,frame,'k')
    title({['Spectra of Frame ' num2str(samp_frames(idx))],''})
    xlabel('Wavelength [nm]')
    ylabel('Intensity [counts]')
    ylim([0 spec_max+(.1*spec_max)])
    hold on
end


% Get band (RGB) value matrix ... (make faster)

n_bands = length(band_num);
n_vpixels = size(frame,1);
rgb_matrix = zeros(n_vpixels,n_frames,n_bands);

for idx_band = 1:n_bands
    for idx_fr = 1:n_frames
        frame = imread([file_names(idx_fr).folder,filesep,file_names(idx_fr).name]);
        for idx_vpix = 1:n_vpixels
            rgb_matrix(idx_vpix,idx_fr,idx_band) = frame(idx_vpix,band_num(idx_band));
        end
    end
end


% Plot RGB spatial reconstruction with noted sample frames

rgb_matrix_tc = (rgb_matrix/spec_max);

figure('Name','Raw RGB Spatial Composite')
imshow(rgb_matrix_tc)
hold on

for i = 1:length(samp_frames)
    line([samp_frames(i) samp_frames(i)],[0 idx_vpix],'Color','red')
end


% Stretch spatial reconstruction for more realistic look

rgb_matrix_tcstretch = imresize(rgb_matrix_tc, [idx_vpix 1800]);

figure('Name','Stretched RGB Spatial Composite')
imshow(rgb_matrix_tcstretch)
%image(rgb_matrix_tcstretch)
