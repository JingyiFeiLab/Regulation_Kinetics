function totalheatmaps_XY(heatmapbinx,heatmapbiny,timelapse,heatmapanalysis_on_or_off)
close all;
% Example of Usage :
% totalheatmaps_XY(30,30,10,'on');

% heatmapbinx : (30 is the bin size in X when I am making the
% histograms of the points projection.....this is now parameterizable);
% heatmapbiny : (30 is the bin size in Y when I am making the histograms of
% the points projection.....this is now parameterizable);
% timelapse : 10 seconds, The time the script allows you to adjust the colormap
% according to your needs
% heatmapanalysis_on_or_off: This allows you all the fancy stuff about adjusting the colormap and settings the upper and lower limit on bins plus checking the 
% actual number of points in each bins, Turn it on if you want such fancy analysis or just turn it off to speeden up the process.

% IMPORTANT EDIT Now this new code includes three major routines
% 1> To show the number of points in each bin/cell along with the colormap
% editor that
% will be on for seconds allowing you to check values in the bin and check
% the minimum and maximum value in the colormap editor, you can also
% manually edit the colormap( BE SURE TO CHECK IMMEDIATE APPLY in the
% colormap editor) this would save the colormap and it will be used to
% further displaing and saving the heatmap.
% 
% 2> Following which the script will ask you to input the lower and upper
% threshold, any bin carring values lower than lowerlimit will be set to
% zero and any bin carrying value greater than upper threshold limit will
% be set to MAXIMUM value of the matrix, this is for the smoothening.
% 
% 2> ColorMap Editor would be evident again this time with the final
% heatmap, you can further edit it for few seconds before the heatmap image
% will be saved.
% 
% This process will be repeated for each XY,XY and XY.
    
%%%% XY

prompt = {'Enter the name that you want to give this total collective set that you are going to combine'};
dlg_title = 'Input';
num_lines = 1;
def = {'Enter any name that you like'};
datasetname = inputdlg(prompt,dlg_title,num_lines,def);
datasetname = datasetname{1};
finalxprojected_StormsXYPlane=[];
finalyprojected_StormsXYPlane=[];
finalzprojected_StormsXYPlane=[];

listing = dir('*.dat');
filenames=listing.name;
listing=listing(1:end);
for j=1:numel(listing)
    
dataset=listing(j,1).name;
data=textread(dataset);
data=data(:,[1 2 3]);
x=data(:,1);
y=data(:,2);
z=data(:,3);
finalxprojected_StormsXYPlane=[finalxprojected_StormsXYPlane ;x];    
finalyprojected_StormsXYPlane=[finalyprojected_StormsXYPlane ;y];  
finalzprojected_StormsXYPlane=[finalzprojected_StormsXYPlane ;z];  
end
close all;

projectedpoints_write=[finalxprojected_StormsXYPlane finalyprojected_StormsXYPlane finalzprojected_StormsXYPlane];
projectedpoints_write_StormsXYPlane=sprintf('%s_Combined_XYPlane.dat',datasetname);
dlmwrite(projectedpoints_write_StormsXYPlane,projectedpoints_write,' ');
close all;

imagewindow = figure('visible','on');
plot(finalxprojected_StormsXYPlane,finalyprojected_StormsXYPlane,'ro');
saveas(imagewindow,'finalprojectpoints_on_StormsXYplane.jpg');
imagewindow = figure('visible','on');
[n,c] = hist3([finalxprojected_StormsXYPlane finalyprojected_StormsXYPlane],[heatmapbinx,heatmapbiny]);
n=rot90(n,1);
n = flipud(n);
noriginal=flipud(n);
if strcmp(heatmapanalysis_on_or_off,'on')==1
close all;
heatmaptext(noriginal); title('Horizontal axis is X and Vertical axis is Y')
colormapeditor
tic;pause(timelapse);toc;
h = gcf;
mycmap = get(h,'Colormap'); 
prompt = {'Enter the value of the cells, anything below which would be set to ZERO for smoothening'};
dlg_title = 'Input the Lower Limit';
num_lines = 1;
def = {'Enter the lower limit number'};
LowerLimitNumber = inputdlg(prompt,dlg_title,num_lines,def);
LowerLimitNumber = LowerLimitNumber{1};
LowerLimitNumber = str2num(LowerLimitNumber);
n(n<LowerLimitNumber)=0;
prompt = {'Enter the value of the cells, anything above which would be set to MAX value of the matrix for smoothening'};
dlg_title = 'Input the Upper Limit';
num_lines = 1;
def = {'Enter the upper limit number'};
UpperLimitNumber = inputdlg(prompt,dlg_title,num_lines,def);
UpperLimitNumber = UpperLimitNumber{1};
UpperLimitNumber = str2num(UpperLimitNumber);
n(n>UpperLimitNumber)=max(max(n));
close all;
imagewindow = figure('visible','on');
imagesc(c{1},c{2},n);
set(gca, 'YDir', 'normal');
xlabel('X Axis');
ylabel('Y Axis');
colorbar;
new_fig=gcf;
set(new_fig,'Colormap',mycmap)
colormapeditor
tic;pause(timelapse);toc;
end
    
if strcmp(heatmapanalysis_on_or_off,'on')==0
imagesc(c{1},c{2},n);
xlabel('X Axis');
ylabel('Y Axis');
set(gca, 'YDir', 'normal');
end
    
projectpoints_on_StormsXYplaneHeatMap=strcat(datasetname,'_TotalFinal_projectpoinmsXYplane_HeatMap.jpg');
saveas(imagewindow,projectpoints_on_StormsXYplaneHeatMap);
projectpoints_on_StormsXYplaneHeatMap=strcat(datasetname,'_TotalFinal_projectpoinmsXYplane_HeatMap.fig');
saveas(imagewindow,projectpoints_on_StormsXYplaneHeatMap);
end