%% Correct lateral sample drift using cross correlation

% Cross correlstion software from:
% Wang, Yina, et al. "Localization events-based sample drift correction for localization microscopy with redundant cross-correlation algorithm." Optics express 22.13 (2014): 15982-15991.
% available at: http://huanglab.ucsf.edu/Resources.html

%% Load the data
clear, clc, close all

path            = 'Z:\data_HTP\2017-05-22_3D_Test_Mito\locResults\Tom20_A647_10ms_2';
filename        = 'Tom20_A647_10ms_2_MMStack_1_Localizations_Z';

cd(path)
locs = dlmread([filename '.csv'], ',',1,0);

file = fopen([filename '.csv']);
line = fgetl(file);
header1 = regexp( line, ',', 'split' );

xCol            = strmatch('x [nm]',header1);
yCol            = strmatch('y [nm]',header1);
framesCol       = strmatch('frame',header1);
LLCol           = strmatch('loglikelihood',header1);
photonsCol      = strmatch('intensity [photon]',header1);
sigmaCol        = strmatch('sigma [nm]',header1);
uncertaintyCol  = strmatch('uncertainty [nm]',header1);

fprintf(' -- Data Loaded -- ')

%% Navigate to the RCC server folder and generate coords variable

pxlsize = 106; % nm

cd('Z:\software\RCC')

coords(:,1) = locs(:,xCol)/pxlsize;
coords(:,2) = locs(:,yCol)/pxlsize;
coords(:,3) = locs(:,framesCol);

%% Run RCC on coord

% Input:    coords:             localization coordinates [x y frame], 
%           segpara:            segmentation parameters (time wimdow, frame)
%           imsize:             image size (pixel)
%           pixelsize:          camera pixel size (nm)
%           binsize:            spatial bin pixel size (nm)

segpara     = 3000;
imsize      = 390;
pixelsize 	= 106;
binsize     = 30;
rmax        = 0.2;

[coordscorr, finaldrift, A,b] = RCC(coords, segpara, imsize, pixelsize, binsize, rmax);

%% Plot final drift curves with corrected data

figure('Position',[100 600 900 400])
subplot(2,1,1)
plot(finaldrift(:,1))
title('x Drift')
subplot(2,1,2)
plot(finaldrift(:,2))
title('y Drift')

figure('Position',[100 400 400 400])
scatter(coords(:,1),coords(:,2),1);
title('Before RCC')

figure('Position',[500 400 400 400])
scatter(coordscorr(:,1),coordscorr(:,2),1);
title('After RCC')

%% Save corrected file as TS csv

locs_DC         = locs;
locs_DC(:,xCol) = coordscorr(:,1);
locs_DC(:,yCol) = coordscorr(:,2);

cd(locpath);tic;

NameCorrected = [filename '_DC_CC.csv'];

fileID = fopen(NameCorrected,'w');
fprintf(fileID,[line ' \n']);
dlmwrite(NameCorrected,locs_DC,'-append');
fclose('all');

fprintf('\n -- Data Saved in %f --\n',toc)





