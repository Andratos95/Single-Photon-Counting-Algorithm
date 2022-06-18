clear all, close all

 %% -------------------------------------------------------------------------------------------------------
date_image = "31-03-2022"; % The date of the measurements
initial_filename_background = "0001_0090_X-ray.tif"; % The file name of the first image to analyze
initial_filename_ross = "0001_0105_X-ray.tif";

xmin = 2e3; % The minimum number of eV to be calculated for the histogram
xmax = 15e3; % The maximum number of eV to be calculated for the histogram
hist_bin_size = 11.7; % The bin size for the histogram in eV (rms of noise * 11.75) 
%% -------------------------------------------------------------------------------------------------------

filenames_background = getFilenames(date_image,initial_filename_background,10);
filenames_ross = getFilenames(date_image,initial_filename_ross,100);
xray_batch_background = filenames2xrayMat(filenames_background);
xray_batch_ross = filenames2xrayMat(filenames_ross);


t1 = tic;

parfor i=1:size(xray_batch_background,2)
    
    coordClusters_background{i}=findCoordClusters(xray_batch_background(i));
    xray_batch_background(i).coordClusters=coordClusters_background{i};
end
parfor i=1:size(xray_batch_ross,2)
    
    coordClusters_ross{i}=findCoordClusters(xray_batch_ross(i));
    xray_batch_ross(i).coordClusters=coordClusters_ross{i};
end
t_findCoordClusters = toc(t1);
fprintf("\nAll clusters coordinates have been found. %d background images, %d ross filter images.\n Elapsed time: %f seconds.\n",size(xray_batch_background,2),size(xray_batch_ross,2),t_findCoordClusters);

%%
%rect_wind_test = [1068 495 1418 845; 1068 1280 1418 1630; 304 890 654 1240; 723 884 1073 1234; 303 1292 653 1642; 309 495 659 845; 1350 261 1700 611 ];
sideRes = 128;
window_size = 256;
rect_windows = genRectWindCoord(2048,sideRes,window_size);
t2 = tic;

[sum_hist_background,~]=fastCoordClustersHist(xray_batch_background,xmin,xmax,hist_bin_size,[1 1 2048 2048]);
%[sum_hist_ross,x_hist]=fastCoordClustersHist(xray_batch_ross,xmin,xmax,hist_bin_size,rect_wind_test);
[sum_hist_ross,x_hist]=fastCoordClustersHist(xray_batch_ross,xmin,xmax,hist_bin_size,rect_windows);
t_coordClustersHist = toc(t2);

fprintf("\nHistograms completed. %d rectangular windows were used.\nElapsed time: %f seconds.\n",size(rect_windows,1),t_coordClustersHist);
%%
t3 = tic;
metal_guess = detectMetals(sum_hist_ross,sum_hist_background,x_hist);
metalChannels = genMultispectralImage(metal_guess,rect_windows);
t_detectMetalsAndGenImages = toc(t3);
fprintf("\nMetals detected.\nElapsed time: %f seconds.\n",t_detectMetalsAndGenImages);
%%
%-----------------------------------------------------------
sourceSize = 2048;
binShift = floor((sourceSize-window_size)/(sideRes-1));
minCropLim = floor(window_size/2);
maxCropLim = floor((binShift*(sideRes-1)+window_size))-minCropLim;
%-----------------------------------------------------------
fig1=figure(1) 
imshow(metalChannels(:,:,1),[],'InitialMagnification', 10000)
title("Titanium channel")
axis equal
axis on
colorbar
xticklabels = minCropLim:254:maxCropLim;
xticks = linspace(1, size(metalChannels(:,:,1), 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = minCropLim:254:maxCropLim;
yticks = linspace(1, size(metalChannels(:,:,1), 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)

fig2=figure(2)
imshow(metalChannels(:,:,2),[],'InitialMagnification', 12000)
title("Copper channel")
axis equal
axis on
colorbar
xticklabels = minCropLim:254:maxCropLim;
xticks = linspace(1, size(metalChannels(:,:,1), 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = minCropLim:254:maxCropLim;
yticks = linspace(1, size(metalChannels(:,:,1), 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)

fig3=figure(3)
imshow(metalChannels(:,:,3),[],'InitialMagnification', 12000)
title("Zinc channel")
axis equal
axis on
colorbar
xticklabels = minCropLim:254:maxCropLim;
xticks = linspace(1, size(metalChannels(:,:,1), 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = minCropLim:254:maxCropLim;
yticks = linspace(1, size(metalChannels(:,:,1), 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
%%
fig4=figure(4)
imshow(metalChannels,'InitialMagnification', 12000)
title("RGB image")
axis equal
axis on
xticklabels = minCropLim:254:maxCropLim;
xticks = linspace(1, size(metalChannels(:,:,1), 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = minCropLim:254:maxCropLim;
yticks = linspace(1, size(metalChannels(:,:,1), 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)

%%
% figure(5)
% subplot(2,1,1);
% plot(x_hist/1000,sum_hist_background/max(sum_hist_background),'.','MarkerSize',10)
% hold on
% active_hist_ross = sum_hist_ross(4,:); % <- CHANGE TO DISPLAY ANOTHER WINDOW!!!
% % 1: Copper
% % 2: Titanium
% % 3: Zinc
% % 4: Background
% % 5: Titanium + zinc overlap
% % 6 : Copper + zinc overlap
% % 7: Copper + Air (half cover test)
% plot(x_hist/1000,active_hist_ross/max(active_hist_ross),'.','MarkerSize',10)
% sum_hist_background_smooth = smoothdata(sum_hist_background,'gaussian',60);
% sum_hist_ross_smooth = smoothdata(active_hist_ross,'gaussian',60);
% plot(x_hist/1000,sum_hist_background_smooth./max(sum_hist_background),'--','LineWidth',1)
% plot(x_hist/1000,sum_hist_ross_smooth./max(active_hist_ross),'b--','LineWidth',1)
% legend ("background","ROI")
% title("10 background sum + 100 ROI sum ")
% 
% subplot(2,1,2);
% ross_tx = (sum_hist_ross_smooth/max(sum_hist_ross_smooth))./(sum_hist_background_smooth/max(sum_hist_background_smooth));
% plot(x_hist/1000,ross_tx);
% title("ROI/background")
% 
% figure(6)
% hist_diff = diff(ross_tx); % Differentiate the measured filter tranmission to see where the K-edge is
% hist_diff(end+1)=hist_diff(end);
% plot(x_hist/1000,hist_diff)
% title("Tx derivative")
