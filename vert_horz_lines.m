% HSI v4 Data Processing: Trondheimsfjorden ground based measurements
% E.F. Prentice
% 7 Aug 2018

% Description: Choose a scene (folder) by specifying the PathName variable.
% Draw a (near) vertical line and (near) horizontal line through the scene
% to get intensity counts of pixels on that line and its two neighbors at
% two desired wavelengths (default is 675 and 683, peak chlorophyll
% flourescence wavelength.

% Research questions: How do the intensity counts from water pixels vary 
% spatially in chlorophyll-sensitive bands?

% Results: Peaks in intensity are pronounced, their origin is not clear.
% Could be chlorophyll or cloud cover or viewing angle or etc.  One trend
% seems to be increased intensity with distance- maybe from the sun?

% Next steps: Make a 3D mesh plot of a square of intensity values at peak
% chlorophyll-sensitive wavelength.

close all; clear all; clc

% Assumptions
load('HSI_wlens.mat')           % scaled wavelength values, called wlens (x-axis vals)
spec_max = 36855;               % max spectral values for scaling
ln = 322;                       % line number for trimming fucntinon [1 650]
tr_pk = 10;                     % number of peaks to consider when trimming
sf = 1800;                      % scale factor for image stretching
band_num = find_band_num([638,543,460,675,683]); % find band number of a desired wavelength
Cv = {'k','c','r',[.5 .6 .7],[.8 .2 .6]};       % color scales for plotting
Ch = {'b','g','m',[.5 .6 .7],[.8 .2 .6]};       % color scales for plotting


% Load data
%PathName = uigetdir('.','Select HSI folder');
%imagefiles = dir(strcat(PathName,[ filesep '*.png']));

PathName = 'C:\Users\elizabep\Documents\fieldwork_data\180724_Trondfj-T2\smiskaret\scan2';
imagefiles = dir(strcat(PathName,[filesep '*.png']));

nFrames = size(imagefiles,1);


% Get wavelength specfic values

frame_sz = imread([imagefiles(1).folder,filesep,imagefiles(1).name]);
nHeight = size(frame_sz,1);
matrix = zeros(nHeight,nFrames,length(band_num));

for idx_fr = 1:nFrames
    frame = imread([imagefiles(idx_fr).folder,filesep,imagefiles(idx_fr).name]);
    for idx_ht = 1:nHeight
        matrix(idx_ht,idx_fr,1) = frame(idx_ht,band_num(1));
        matrix(idx_ht,idx_fr,2) = frame(idx_ht,band_num(2));
        matrix(idx_ht,idx_fr,3) = frame(idx_ht,band_num(3));
        matrix(idx_ht,idx_fr,4) = frame(idx_ht,band_num(4));
        matrix(idx_ht,idx_fr,5) = frame(idx_ht,band_num(5));
    end
end


% Trim repeated-pixel edges of spatial reconstruction

[m,idx_peaks] = maxk(diff(matrix(ln,:,length(band_num))),tr_pk,'ComparisonMethod','abs');
cut1 = min(idx_peaks);
cut2 = max(idx_peaks);

matrix_trimmed = (matrix(:,cut1:cut2,:)/spec_max);


% Show spatial rep of scans

matrix_stretch = imresize(matrix_trimmed,[idx_ht sf]);

figure('Name','Spectral Signature of Selected Pixels')

subplot(3,1,1)
imshow(matrix_stretch(:,:,1:3))
title(['R = ' num2str(wlens(band_num(1))) ' nm, G = ' num2str(wlens(band_num(2))) ' nm, B = ' num2str(wlens(band_num(3))) ' nm'])

subplot(3,1,2)
imshow(matrix_stretch(:,:,4))
title([num2str(wlens(band_num(4))) ' nm'])

subplot(3,1,3)
imshow(matrix_stretch(:,:,5))
title([num2str(wlens(band_num(5))) ' nm'])


% Visually select desired line to analyze

figure('Name','Selcting Lines')

subplot(1,2,1)
imshow(matrix_trimmed(:,:,1:3))
title({'click and drag to draw approx VERTICAL line,','modify, double-click line to finalize'},'FontSize',8)
ln_vert = imline;
pts_vert = wait(ln_vert);

subplot(1,2,2)
imshow(matrix_trimmed(:,:,1:3))
title({'click and drag to draw approx HORIZONTAL line,','modify, double-click line to finalize'},'FontSize',8)
ln_horz = imline;
pts_horz = wait(ln_horz);


% Plot full spectral signature of selected points (pixels)

figure('Name','Spectral Signature of Selected Pixels')

trimmed_size = size(matrix_trimmed);
wl_scale = linspace(min(wlens),max(wlens),180); %%%%

x_vert = round(mean(pts_vert(:,1)))+(cut1-1); % frame
x_vert_sf = sf*(round(mean(pts_vert(:,1)))/trimmed_size(1,2)); % frame (scaled)
y1_vert = round(pts_vert(1,2));             % start vertical pixel (top)
y2_vert = round(pts_vert(2,2));             % end vertical pixel (bottom)

subplot(2,3,1:2) % spatial recon w/vertical line
imshow(matrix_stretch(:,:,1:3))
line([x_vert_sf-1 x_vert_sf-1],[y1_vert y2_vert],'Color',Cv{1},'LineWidth',1)
line([x_vert_sf x_vert_sf],[y1_vert y2_vert],'Color',Cv{2},'LineWidth',1)
line([x_vert_sf+1 x_vert_sf+1],[y1_vert y2_vert],'Color',Cv{3},'LineWidth',1)
title('Vertical Selection')

x1_horz_sf = sf*(round(pts_horz(1,1))/trimmed_size(1,2)); % start frame, scaled (left)
x2_horz_sf = sf*(round(pts_horz(2,1))/trimmed_size(1,2)); % end frame, scaled (right)
x1_horz = round(pts_horz(1,1));             % start frame (left)
x2_horz = round(pts_horz(2,1));             % end frame (right)
y_horz = round(mean(pts_horz(:,2)));        % vertical pixel

subplot(2,3,4:5) % spatial recon w/horizontal line
imshow(matrix_stretch(:,:,1:3))
line([x1_horz_sf x2_horz_sf],[y_horz-1 y_horz-1],'Color',Ch{1},'LineWidth',1)
line([x1_horz_sf x2_horz_sf],[y_horz y_horz],'Color',Ch{2},'LineWidth',1)
line([x1_horz_sf x2_horz_sf],[y_horz+1 y_horz+1],'Color',Ch{3},'LineWidth',1)
title('Horizontal Selection')

i = 0;
subplot(2,3,3) % spectral sig plot of vertical line pix
for idx_frmv = x_vert-1:x_vert+1
    i = i+1;
    framev = imread([imagefiles(idx_frmv).folder,filesep,imagefiles(idx_frmv).name]);
    for idx_vpixv = y1_vert:y2_vert
        plot(wl_scale,framev(idx_vpixv,:),'color',Cv{i})
        hold on
    end
end

title({'Spectral Signature of 3 Consecutive VERTICAL Lines',''},'FontSize',8)
xlabel('Wavelength [nm]')
ylabel('Intensity [counts]')
ylim([0 4e4])

j = 0;
subplot(2,3,6) % spectral sig plot of horizontal line pix
for idx_vpixh = y_horz-1:y_horz+1
    j = j+1;
    for idx_frmh = x1_horz+(cut1-1):x2_horz+(cut1-1)
        frameh = imread([imagefiles(idx_frmh).folder,filesep,imagefiles(idx_frmh).name]);
        plot(wl_scale,frameh(idx_vpixh,:),'color',Ch{j})
        hold on
    end
end

title({'Spectral Signature of 3 Consecutive HORIZONTAL Lines',''},'FontSize',8)
xlabel('Wavelength [nm]')
ylabel('Intensity [counts]')
ylim([0 4e4])

clear i j idx_vpixv idx_vpixh idx_frmv idx_frmh framev frameh
clear x_vert x_vert_sf y1_vert y2_vert
clear x1_horz_sf x2_horz_sf x1_horz x2_horz y_horz


% Plot intensity of selected pixels at chlorophyll flourescence wavelength
% VERTICAL

figure('Name','Intensity at chlorophyll flourescence wavelength (VERTICAL line)')

x_vert = round(mean(pts_vert(:,1)));        % frame
y1_vert = round(pts_vert(1,2));             % start vertical pixel (top)
y2_vert = round(pts_vert(2,2));             % end vertical pixel (bottom)

subplot(2,2,1:2) % spatial recon w/vertical line
image(matrix_trimmed(:,:,1:3))
line([x_vert-1 x_vert-1],[y1_vert y2_vert],'Color',Cv{1},'LineWidth',1)
line([x_vert x_vert],[y1_vert y2_vert],'Color',Cv{2},'LineWidth',1)
line([x_vert+1 x_vert+1],[y1_vert y2_vert],'Color',Cv{3},'LineWidth',1)
title('Vertical Selection')

subplot(2,2,3) % intensity plot @wv_1, vertical
ii = 0;
for idx_frm = (x_vert+(cut1-1))-1:(x_vert+(cut1-1))+1
    ii = ii+1;
    framev = imread([imagefiles(idx_frm).folder,filesep,imagefiles(idx_frm).name]);
    for idx_vpix = y1_vert:y2_vert
        x_ax = idx_vpix;
        y_ax = framev(idx_vpix,band_num(4));
        plot(x_ax,y_ax,'Color',Cv{ii},'Marker','.','MarkerSize',12)
        hold on
    end
end

title({['Intensity @' num2str(wlens(band_num(4))) ' nm'],['of 3 Consecutive VERTICAL Lines']})
xlabel('Vertical Pixel Number')
ylabel('Intensity [counts]')

false_plot = zeros(3,1);
false_plot(1) = plot(NaN,NaN,'Color',Cv{1},'Marker','.','MarkerSize',12);
false_plot(2) = plot(NaN,NaN,'Color',Cv{2},'Marker','.','MarkerSize',12);
false_plot(3) = plot(NaN,NaN,'Color',Cv{3},'Marker','.','MarkerSize',12);
legend(false_plot,'left line','mid line','right line')

clear ii idx_vpix idx_frm framev x_ax y_ax false_plot

subplot(2,2,4) % intensity plot @wv_2, vertical
iii = 0;
for idx_frm = (x_vert+(cut1-1))-1:(x_vert+(cut1-1))+1
    iii = iii+1;
    framev = imread([imagefiles(idx_frm).folder,filesep,imagefiles(idx_frm).name]);
    for idx_vpix = y1_vert:y2_vert
        x_ax = idx_vpix;
        y_ax = framev(idx_vpix,band_num(5));
        plot(x_ax,y_ax,'Color',Cv{iii},'Marker','.','MarkerSize',12)
        hold on
    end
end

title({['Intensity @' num2str(wlens(band_num(5))) ' nm'],['of 3 Consecutive VERTICAL Lines']})
xlabel('Vertical Pixel Number')
ylabel('Intensity [counts]')

false_plot = zeros(3,1);
false_plot(1) = plot(NaN,NaN,'Color',Cv{1},'Marker','.','MarkerSize',12);
false_plot(2) = plot(NaN,NaN,'Color',Cv{2},'Marker','.','MarkerSize',12);
false_plot(3) = plot(NaN,NaN,'Color',Cv{3},'Marker','.','MarkerSize',12);
legend(false_plot,'left line','mid line','right line')

clear iii idx_vpix idx_frm framev x_ax y_ax false_plot
clear x_vert y1_vert y2_vert


% Plot intensity of selected pixels at chlorophyll flourescence wavelength
% HORIZONTAL

figure('Name','Intensity at chlorophyll flourescence wavelength (HORIZONTAL line)')

x1_horz = round(pts_horz(1,1));         % start frame (left)
x2_horz = round(pts_horz(2,1));         % end frame (right)
y_horz = round(mean(pts_horz(:,2)));    % vertical pixel

subplot(2,2,1:2) % spatial recon w/horizontal line
image(matrix_trimmed(:,:,1:3))
line([x1_horz x2_horz],[y_horz-1 y_horz-1],'color',Ch{1},'LineWidth',1)
line([x1_horz x2_horz],[y_horz y_horz],'color',Ch{2},'LineWidth',1)
line([x1_horz x2_horz],[y_horz+1 y_horz+1],'color',Ch{3},'LineWidth',1)
title('Horizontal Selection')

subplot(2,2,3) % intensity plot @wv_1, horizontal
jj = 0;
for idx_vpix = y_horz-1:y_horz+1
    jj = jj+1;
    for idx_frm = x1_horz+(cut1-1):x2_horz+(cut1-1)
        frameh = imread([imagefiles(idx_frm).folder,filesep,imagefiles(idx_frm).name]);
        x_ax = idx_frm-(cut1-1);
        y_ax = frameh(idx_vpix,band_num(4));
        plot(x_ax,y_ax,'color',Ch{jj},'Marker','.','MarkerSize',12)
        hold on
    end
end

title({['Intensity @' num2str(wlens(band_num(4))) ' nm'],['of 3 Consecutive HORIZONTAL Lines']})
xlabel('Horizontal Pixel Number (Frame #)')
ylabel('Intensity [counts]')

false_plot = zeros(3,1);
false_plot(1) = plot(NaN,NaN,'Color',Ch{1},'Marker','.','MarkerSize',12);
false_plot(2) = plot(NaN,NaN,'Color',Ch{2},'Marker','.','MarkerSize',12);
false_plot(3) = plot(NaN,NaN,'Color',Ch{3},'Marker','.','MarkerSize',12);
legend(false_plot,'upper line','mid line','lower line')

clear jj idx_vpix idx_frm frameh x_ax y_ax false_plot

subplot(2,2,4) % intensity plot @wv_2, horizontal
jjj = 0;
for idx_vpix = y_horz-1:y_horz+1
    jjj = jjj+1;
    for idx_frm = x1_horz+(cut1-1):x2_horz+(cut1-1)
        frameh = imread([imagefiles(idx_frm).folder,filesep,imagefiles(idx_frm).name]);
        x_ax = idx_frm-(cut1-1);
        y_ax = frameh(idx_vpix,band_num(5));
        plot(x_ax,y_ax,'color',Ch{jjj},'Marker','.','MarkerSize',12)
        hold on
    end
end

title({['Intensity @' num2str(wlens(band_num(5))) ' nm'],['of 3 Consecutive HORIZONTAL Lines']})
xlabel('Horizontal Pixel Number (Frame #)')
ylabel('Intensity [counts]')

false_plot = zeros(3,1);
false_plot(1) = plot(NaN,NaN,'Color',Ch{1},'Marker','.','MarkerSize',12);
false_plot(2) = plot(NaN,NaN,'Color',Ch{2},'Marker','.','MarkerSize',12);
false_plot(3) = plot(NaN,NaN,'Color',Ch{3},'Marker','.','MarkerSize',12);
legend(false_plot,'upper line','mid line','lower line')

clear jjj idx_vpix idx_frm frameh x_ax y_ax false_plot
clear x1_horz x2_horz y_horz

