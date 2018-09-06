% HSI v4 Data Processing: Trondheimsfjorden ground based measurements
% E.F. Prentice
% 9 Aug 2018

% Description: Make a 3D mesh plot of a square of intensity values at peak
% chlorophyll-sensitive wavelength.

% Research questions: How do the intensity counts from water pixels vary 
% spatially (in 3D) in chlorophyll-sensitive bands?

% Results: very difficult to say what we see

% Next steps: get radience!!!

close all; clear all; clc

% Assumptions
load('HSI_wlens.mat')           % scaled wavelength values, called wlens (x-axis vals)
spec_max = 36855;               % max spectral values for scaling
ln = 322;                       % line number for trimming fucntinon [1 650]
tr_pk = 10;                     % number of peaks to consider when trimming
sf = 1800;                      % scale factor for image stretching
band_num = find_band_num([638,543,460,675,683]); % find band number of a desired wavelength


% Load data
%PathName = uigetdir('.','Select HSI folder');
%imagefiles = dir(strcat(PathName,[ filesep '*.png']));

PathName = 'C:\Users\elizabep\Documents\fieldwork_data\180724_Trondfj-T2\smiskaret\scan3';
imagefiles = dir(strcat(PathName,[filesep '*.png']));

nFrames = size(imagefiles,1);
frame_sz = imread([imagefiles(1).folder,filesep,imagefiles(1).name]);
nHeight = size(frame_sz,1);


% Get wavelength specfic values

matrix = zeros(nHeight,nFrames,length(band_num));
for bn = 1:length(band_num)
    for idx_fr = 1:nFrames
        frame = imread([imagefiles(idx_fr).folder,filesep,imagefiles(idx_fr).name]);
        for idx_ht = 1:nHeight
            matrix(idx_ht,idx_fr,bn) = frame(idx_ht,band_num(bn));
        end
    end
end

clear bn idx_fr idx_ht frame_sz


% Trim repeated-pixel edges of spatial reconstruction

[~,idx_peaks] = maxk(diff(matrix(ln,:,length(band_num))),tr_pk,'ComparisonMethod','abs');
cut1 = min(idx_peaks);
cut2 = max(idx_peaks);

matrix_trimmed = (matrix(:,cut1:cut2,:)/spec_max);

clear idx_peaks


% Visually select desired area to analyze

figure('Name','Selcting a Box')
image(matrix_trimmed(:,:,1:3))
title({'click and drag to draw desired box,','modify, double-click in box to finalize'},'FontSize',8)
box = imrect;
bx = wait(box);
pts_rect = [bx(1),bx(2);
            bx(1)+bx(3),bx(2);
            bx(1)+bx(3),bx(2)+bx(4);
            bx(1),bx(2)+bx(4)];

clear bx box


% Create mesh of possible pixels

size_trim = size(matrix_trimmed);
Height_trim = 1:size_trim(1);
Frames_trim = 1:size_trim(2);
[xx,yy] = meshgrid(Frames_trim,Height_trim);
pts_mesh = [xx(:),yy(:)];

clear size_trim Height_trim Frames_trim xx yy


% Plot points contained in desired rectangle

in = inpolygon(pts_mesh(:,1),pts_mesh(:,2),pts_rect(:,1),pts_rect(:,2));
idx_in = find(in);
pts_inrect = pts_mesh(idx_in,:);

hold on
plot(pts_inrect(:,1),pts_inrect(:,2),'r.')

clear in idx_in


% Get intensity values for each xy pair of points in selected rectangle
% @wl#1

lbnd_fr = min(pts_inrect(:,1));
ubnd_fr = max(pts_inrect(:,1));
lbnd_pix = min(pts_inrect(:,2));
ubnd_pix = max(pts_inrect(:,2));
sz_ptsrect = size(pts_inrect);

intens = zeros(sz_ptsrect(1),1);
cnt = 0;
for idx_fr = lbnd_fr:ubnd_fr
    frame2 = imread([imagefiles(idx_fr+(cut1-1)).folder,filesep,imagefiles(idx_fr+(cut1-1)).name]);
    idx_pixinfr = find(pts_inrect(:,1) == idx_fr);
    for idx_vpix = lbnd_pix:ubnd_pix
        cnt = cnt+1;
        intens(cnt,1) = frame2(idx_vpix,band_num(4));
    end
end

[x,y] = meshgrid(lbnd_fr:ubnd_fr,lbnd_pix:ubnd_pix);
grid_sz = size(x);
intens_grid = reshape(intens,[grid_sz(1),grid_sz(2)]);

figure
h1 = axes;
mesh(x,y,intens_grid)
title({['Intensity @' num2str(wlens(band_num(4))) ' nm'],['of selected square area']})
set(h1,'Ydir','reverse')
xlabel('Horizontal Pixel Number (Frame)')
ylabel('Vertical Pixel Number')
zlabel('Intensity [counts]')
zlim([3e3 8e3])
caxis([3e3 8e3])

clear intens cnt idx_fr frame2 idx_pixinfr idx_vpix x y grid_sz intens_grid h1


% Get intensity values for each xy pair of points in selected rectangle
% @wl#2

intens = zeros(sz_ptsrect(1),1);
cnt = 0;
for idx_fr = lbnd_fr:ubnd_fr
    frame2 = imread([imagefiles(idx_fr+(cut1-1)).folder,filesep,imagefiles(idx_fr+(cut1-1)).name]);
    idx_pixinfr = find(pts_inrect(:,1) == idx_fr);
    for idx_vpix = lbnd_pix:ubnd_pix
        cnt = cnt+1;
        intens(cnt,1) = frame2(idx_vpix,band_num(5));
    end
end

[x,y] = meshgrid(lbnd_fr:ubnd_fr,lbnd_pix:ubnd_pix);
grid_sz = size(x);
intens_grid = reshape(intens,[grid_sz(1),grid_sz(2)]);

figure
h1 = axes;
mesh(x,y,intens_grid)
title({['Intensity @' num2str(wlens(band_num(5))) ' nm'],['of selected square area']})
set(h1,'Ydir','reverse')
xlabel('Horizontal Pixel Number (Frame)')
ylabel('Vertical Pixel Number')
zlabel('Intensity [counts]')
zlim([3e3 8e3])
caxis([3e3 8e3])
% FIX THESE LIMITS TO BE AUTOMATIC
