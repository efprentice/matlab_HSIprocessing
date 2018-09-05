% HSI v4 Data Processing: Trondheimsfjorden ground based measurements
% E.F. Prentice
% 2 Aug 2018

% Description: Choose a scene (folder) by specifying the PathName variable.
% Decide on your tool for selection (imfreehand, impoly, etc) or variable
% h. Select a small portion of the scene and view intensity counts of that
% selection (or class, if chosen).

% Research questions: Can just pixels of a certain class be selected
% instead of viewing all spectrum per frame? Are the spectra of different
% classes discernable?

% Results: Yes, with GUI and spatial reconstruction of scene. Spectra of
% vegetation v sky are very different but water and algae water still look
% quite similar.

% Next steps: Intensity counts aren't very helpful, get reflectance or
% radiance. Compare wavelengths other than RGB for biology. Look at pixels
% near each other.

clear all; clc

% Assumptions
load('HSI_wlens.mat') % scaled wavelength values, called wlens (x-axis vals)
spec_max = 36855;     % max spectral values for scaling
bandnum = [82,54,30];         % default [80,56,33], I like [82,54,30]
ln = 322;                     % line number for trimming fucntinon [1 650]
tr_pk = 10;                   % number of peaks to consider when trimming
wl_red = bandnum(1);          % Red wavelength band number (of 180) = 630 nm default
wl_green = bandnum(2);        % Green wavelength band number (of 180) = 550 nm default
wl_blue = bandnum(3);         % Blue wavelength band number (of 180) = 470 nm default


% Load data
%PathName = uigetdir('.','Select HSI folder');
%imagefiles = dir(strcat(PathName,[ filesep '*.png']));

PathName = 'C:\Users\elizabep\Desktop\vmshare\test';
%PathName = 'C:\Users\elizabep\Documents\fieldwork_data\180726_Trondfj-T4\vattakammen\scan1';
imagefiles = dir(strcat(PathName,[filesep '*.png']));
nFrames = size(imagefiles,1);


% Get RGB values
frame_sz = imread([imagefiles(1).folder,filesep,imagefiles(1).name]);
nHeight = size(frame_sz,1);
rgb_matrix = zeros(nHeight,nFrames,3);
for idx_fr = 1:nFrames
    frame = imread([imagefiles(idx_fr).folder,filesep,imagefiles(idx_fr).name]);
    for idx_ht = 1:nHeight
        rgb_matrix(idx_ht,idx_fr,1) = frame(idx_ht,wl_red);
        rgb_matrix(idx_ht,idx_fr,2) = frame(idx_ht,wl_green);
        rgb_matrix(idx_ht,idx_fr,3) = frame(idx_ht,wl_blue);
    end
end

% Trim repeated-pixel edges of RGB spatial reconstruction
[m,idx_peaks] = maxk(diff(rgb_matrix(ln,:,3)),tr_pk,'ComparisonMethod','abs');
cut1 = min(idx_peaks);
cut2 = max(idx_peaks);

rgb_matrix_trimmed = (rgb_matrix(:,cut1:cut2,:)/spec_max);

% Show spatial rep of scans
mini_sz = size(rgb_matrix_trimmed);
figure('Name','Spectral Signature of Selected Pixels')
subplot(1,2,1)
image(rgb_matrix_trimmed)
title({'click and drag to highlight desired area,','double-click inside shape to finalize'},'FontSize',8)

% Visually select desired area to analyze
h = imfreehand; % or imfreehand or impoly or imrect(?) or imellipse
pts_poly = wait(h);

% Create mesh of possible pixels
pix_vert = 1:mini_sz(1);
pix_horz = 1:mini_sz(2);
[xx,yy] = meshgrid(pix_horz,pix_vert);
pts_mesh = [xx(:),yy(:)];

% Plot points contained in desired area
in = inpolygon(pts_mesh(:,1),pts_mesh(:,2),pts_poly(:,1),pts_poly(:,2));
idx_in = find(in);
pts_inpoly = pts_mesh(idx_in,:);

hold on
plot(pts_inpoly(:,1),pts_inpoly(:,2),'r.')

% Plot spectral signature of selected points (pixels)
lbnd_fr = min(pts_inpoly(:,1))+(cut1-1);
ubnd_fr = max(pts_inpoly(:,1))+(cut1-1);
wl_scale = linspace(min(wlens),max(wlens),180); %%%%

%figure('Name','Spectral Signature of Selected Pixels')
subplot(1,2,2)
for idx_fr2 = lbnd_fr:ubnd_fr
    frame2 = imread([imagefiles(idx_fr2).folder,filesep,imagefiles(idx_fr2).name]);
    idx_pixinfr = find(pts_inpoly(:,1) == idx_fr2-(cut1-1));
    for idx_vpix = min(idx_pixinfr):max(idx_pixinfr)
        plot(wl_scale,frame2(pts_inpoly(idx_vpix,2),:),'k')
        title({'Spectral Signature of Selected Pixels',''},'FontSize',8)
        xlabel('Wavelength [nm]')
        ylabel('Intensity [counts]')
        ylim([0 4e4])
        hold on
    end
end
    




