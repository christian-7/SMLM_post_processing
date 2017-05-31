clear, clc, close all

locpath             = 'Z:\Christian-Sieben\data_HTP\2017-05-19_3D_Test_Centriole\locResults_Cent\DL755_COT_500mW_20ms_1';
locname             = 'DL755_COT_500mW_20ms_1_MMStack_1_Localizations_DC_corrected';

zCalibrationPath    = 'Z:\Christian-Sieben\data_HTP\2017-05-19_3D_Test_Centriole\Calibration';
zCal                = 'splineFit_Ch2.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(locpath)

locs_Ch1        = dlmread([locname '.csv'],',',1,0);

file = fopen(['DL755_COT_500mW_20ms_1_MMStack_1_Localizations_DC.csv']);
% file = fopen([locname '.csv']);
line = fgetl(file);
h = regexp( line, ',', 'split' );

xCol                = strmatch('x [nm]',h);
yCol                = strmatch('y [nm]',h);
frameCol            = strmatch('frame',h);
LLCol               = strmatch('loglikelihood',h);
photonsCol          = strmatch('intensity [photon]',h);
sigmaXCol           = strmatch('sigma_x [nm]',h);
sigmaYCol           = strmatch('sigma_y [nm]',h);
uncertaintyCol      = strmatch('uncertainty [nm]',h);
zCol                = size(locs_Ch1,2)+1;

cd(zCalibrationPath);
load(zCal);

fprintf('\n -- Data Loaded --\n')

%% Calculate delta Sigma and calculate Z position

deltaSigma = [];
deltaSigma = locs_Ch1(:,sigmaXCol) - locs_Ch1(:,sigmaYCol);

filter          = find(deltaSigma < maxD_Ch2 & deltaSigma > minD_Ch2);
locs_Zfiltered  = locs_Ch1(filter,1:end);

locs_Zfiltered(:,zCol) = deltaF_Ch2(deltaSigma(filter,1));

fprintf('\n -- Calculated Z position --\n')

%% Save Localization file

cd(locpath);tic;

NameCorrected = [locname '_Z.csv'];

fileID = fopen(NameCorrected,'w');
fprintf(fileID,[[line,',z [nm]'] ' \n']);
dlmwrite(NameCorrected,locs_Zfiltered,'-append');
fclose('all');

fprintf('\n -- Data Saved in %f --\n',toc)

%% Save for http://www.cake23.de/pointcloud-loader/

% The file should be composed as: 
% frame	x [nm]	y [nm]	z [nm]	uncertainty_xy [nm]	uncertainty_z [nm]

newVar = [];
newVar(:,1) = locs_Zfiltered(:,frameCol);
newVar(:,2) = locs_Zfiltered(:,xCol);
newVar(:,3) = locs_Zfiltered(:,yCol);
newVar(:,4) = locs_Zfiltered(:,zCol);
newVar(:,5) = locs_Zfiltered(:,uncertaintyCol);
newVar(:,6) = 40;

cd(locpath);
cd('C:\Users\sieben\Desktop')
tic;

NameCorrected = [locname '_openGL.csv'];

new_header = ['frame,x [nm],y [nm],z [nm],uncertainty_xy [nm],uncertainty_z [nm]'];

fileID = fopen(NameCorrected,'w');
fprintf(fileID,[new_header ' \n']);
dlmwrite(NameCorrected,newVar,'-append');
fclose('all');

fprintf('\n -- Data Saved in %f --\n',toc)
