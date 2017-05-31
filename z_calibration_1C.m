clear, clc, close all

locpath     = 'Z:\Christian-Sieben\data_HTP\2017-05-22_3D_Test_Mito\3D_Calibration';

bead_stack_Ch1     = 'z_stack_beads_40nm_Ch642_MMStack_Pos0_Localizations_BeadStack';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(locpath)

locs_Ch1        = dlmread([bead_stack_Ch1 '.csv'],',',1,0);

file = fopen([bead_stack_Ch1 '.csv']);
line = fgetl(file);
h = regexp( line, ',', 'split' );

xCol                = strmatch('x [nm]',h);
yCol                = strmatch('y [nm]',h);
zCol                = strmatch('z [nm]',h);
LLCol               = strmatch('loglikelihood',h);
photonsCol          = strmatch('intensity [photon]',h);
sigmaXCol           = strmatch('sigma_x [nm]',h);
sigmaYCol           = strmatch('sigma_y [nm]',h);
uncertaintyCol      = strmatch('uncertainty [nm]',h);

fprintf('\n -- Data Loaded --\n')

%% Compare with Calibration curve from Huang et al.

% Huang, Bo, et al. "Three-dimensional super-resolution imaging by stochastic optical reconstruction microscopy." Science 319.5864 (2008): 810-813.

load('Z:\Christian-Sieben\data_HTP\2017-05-16_3D_tests\huang.mat')
figure
scatter(huang(:,3),huang(:,2)-huang(:,1),30,'filled');hold on
scatter(locs_Ch1(:,zCol),locs_Ch1(:,sigmaYCol)-locs_Ch1(:,sigmaXCol),1)
legend('Own','Huang et al')

%% Select an single bead

% Plot the data
% Show 2D histogram 
% Select Au Fuducial using rectangular selection

close all
% 
% figure('Position',[100 400 500 500])
% scatter(locs(:,xcol),locs(:,ycol),1);

pxlsize = 50; 

heigth  = round((max(locs_Ch1(:,yCol))-min(locs_Ch1(:,yCol)))/pxlsize);
width   = round((max(locs_Ch1(:,xCol))-min(locs_Ch1(:,xCol)))/pxlsize);

figure('Position',[100 200 800 800])
im=hist3([locs_Ch1(:,xCol),locs_Ch1(:,yCol)],[width heigth]); % heigth x width
imagesc(imrotate(im,90),[0 5]);
% imagesc(im,[0 200]);
colormap('hot');
% colorbar
rect = getrect; % rect = [xmin ymin width height];

close all

fprintf('\n -- ROI selected --\n')

% Select ROI for both channels

xmin = min(locs_Ch1(:,xCol))+ rect(:,1)*pxlsize;
ymin = max(locs_Ch1(:,yCol)) - rect(:,2)*pxlsize - (rect(:,4)*pxlsize) ;
xmax = xmin + (rect(:,3)* pxlsize);
ymax = ymin + rect(:,4) * pxlsize;

vx=find(locs_Ch1(:,xCol)>xmin & locs_Ch1(:,xCol)<xmax);
subset1=locs_Ch1(vx,1:end);
vy=find(subset1(:,yCol)>ymin & subset1(:,yCol)<ymax);
subset2=subset1(vy,1:end);

figure('Position',[200 600 300 300])
scatter(subset2(:,zCol),subset2(:,sigmaXCol),2);hold on;
scatter(subset2(:,zCol),subset2(:,sigmaYCol),2);hold on;
title('Ch1 Ex: 642 nm')
box on;

single_bead_Ch1 = subset2;

fprintf('\n -- Plotted selected ROI  --\n')

%% Select Range to plot each bead

minZ_Ch1 = -800;
maxZ_Ch1 = 800;

figure('Position',[200 200 300 300])
scatter(single_bead_Ch1(single_bead_Ch1(:,zCol)>minZ_Ch1 & single_bead_Ch1(:,zCol)<maxZ_Ch1,zCol),single_bead_Ch1(single_bead_Ch1(:,zCol)>minZ_Ch1 & single_bead_Ch1(:,zCol)<maxZ_Ch1,sigmaXCol),5,'filled');hold on;
scatter(single_bead_Ch1(single_bead_Ch1(:,zCol)>minZ_Ch1 & single_bead_Ch1(:,zCol)<maxZ_Ch1,zCol),single_bead_Ch1(single_bead_Ch1(:,zCol)>minZ_Ch1 & single_bead_Ch1(:,zCol)<maxZ_Ch1,sigmaYCol),5,'filled');hold on;
legend('sigma X','sigma Y'); box on; title('Ch1 Ex: 642 nm');

%% Fit the curves for Channel 1

% Channel 1, Wx
fx  = polyfit(single_bead_Ch1(single_bead_Ch1(:,zCol)>minZ_Ch1 & single_bead_Ch1(:,zCol)<maxZ_Ch1,zCol),single_bead_Ch1(single_bead_Ch1(:,zCol)>minZ_Ch1 & single_bead_Ch1(:,zCol)<maxZ_Ch1,sigmaXCol),6);
x1_Ch1 = minZ_Ch1:50:maxZ_Ch1;
y1 = polyval(fx,x1_Ch1);

% Channel 1, Wy
fy  = polyfit(single_bead_Ch1(single_bead_Ch1(:,zCol)>minZ_Ch1 & single_bead_Ch1(:,zCol)<maxZ_Ch1,zCol),single_bead_Ch1(single_bead_Ch1(:,zCol)>minZ_Ch1 & single_bead_Ch1(:,zCol)<maxZ_Ch1,sigmaYCol),6);
y2 = polyval(fy,x1_Ch1);

figure('Position',[100 500 400 400])
scatter(single_bead_Ch1(single_bead_Ch1(:,zCol)>minZ_Ch1 & single_bead_Ch1(:,zCol)<maxZ_Ch1,zCol),single_bead_Ch1(single_bead_Ch1(:,zCol)>minZ_Ch1 & single_bead_Ch1(:,zCol)<maxZ_Ch1,sigmaXCol));hold on
plot(x1_Ch1,y1,'LineWidth',2);hold on;
scatter(single_bead_Ch1(single_bead_Ch1(:,zCol)>minZ_Ch1 & single_bead_Ch1(:,zCol)<maxZ_Ch1,zCol),single_bead_Ch1(single_bead_Ch1(:,zCol)>minZ_Ch1 & single_bead_Ch1(:,zCol)<maxZ_Ch1,sigmaYCol));hold on
plot(x1_Ch1,y2,'LineWidth',2);
box on;axis square;
axis([minZ_Ch1 maxZ_Ch1 100 500])
title('Polynomial fit of W_x and W_y');
xlabel('z position [nm]')
ylabel('W_x and W_y [nm]');
 
delta = transpose(polyval(fx,x1_Ch1)-polyval(fy,x1_Ch1));
deltaF = fit(transpose(x1_Ch1),delta,'smoothingspline','SmoothingParam',0.07);

figure('Position',[700 500 400 400])
scatter(x1_Ch1,delta)
axis square
plot(deltaF,x1_Ch1,delta);
title('Spline Fit of W_x - W_y');
ylabel('W_x - W_y [nm]');
xlabel('z position [nm]')

% Inverse fitting of selected region

minZ_Ch1_fit = -500;
maxZ_Ch1_fit = 500;

x1_Ch1 = minZ_Ch1_fit:50:maxZ_Ch1_fit;

delta = transpose(polyval(fx,x1_Ch1)-polyval(fy,x1_Ch1));
deltaF_Ch1 = fit(delta,transpose(x1_Ch1),'smoothingspline','SmoothingParam',0.07);

figure('Position',[1200 500 400 400])
scatter(delta,x1_Ch1); hold on;
axis square
plot(deltaF_Ch1,delta,x1_Ch1);
title('Inverted spline fit of W_x - W_y');
xlabel('W_x - W_y [nm]');
ylabel('z position [nm]');
box on;

%% Save fit

close all

minD = min(delta);
maxD = max(delta);

save('splineFit_Ch1.mat','deltaF_Ch1','minD','maxD');

% Use: z = deltaF(Wx-Wy);
