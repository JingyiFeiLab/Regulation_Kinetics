function projectpoints_along_cell_axis_oldd(pixelingfactor,front_distancepercentage,back_distancepercentage,histogrambinsize,heatmapbinx,heatmapbiny,timelapse,heatmapanalysis_on_or_off)
% Example of Usage : projectpoints_along_cell_axis_oldd(1,20,20,30,30,30,10,'off');
% This function will read all the .txt files in the folder and project
% their points and overlay them It can also include the fourth column for
% the color codes.

% pixelingfactor :1, If you want to convert between nm/Angstrom and pixel
% cooordinates, then you can mention the appropriate pixelingfactor.

% front_distancepercentage :(20 indicates 20% of the points from the front)
% The initial few points may be discarded..so you can mention the distance
% that needs to be truncated from the beginning. 
% back_distancepercentage% :(20 indicates 20% of the points from the front) 
% The initial few points may be discarded..so you can mention the distance  
% that needs to be truncated from the end. 
% histogrambinsize : (30 is the bin size for histogram); 
% heatmapbinx : (30 is the bin size in X when I am making the
% histograms of the points projection.....this is now parameterizable);
% heatmapbiny : (30 is the bin size in Y when I am making the histograms of
% the points projection.....this is now parameterizable);
% timelapse : The time the script allows you to adjust the colormap
% according to your needs
% heatmapanalysis_on_or_off: This allows you all the fancy stuff about adjusting the colormap and settings the upper and lower limit on bins plus checking the 
% actual number of points in each bins, Turn it on if you want such fancy analysis or just turn it off to speeden up the process.


% Important Notes : This function stores the projection points (all the
% three planes) for each individual cell files in a folder called
% EACHINDIVIDUALCELL in the same working directory, It has sub folders
% called XY,YZ,XZ.

% These individual files then can be combined for combined heatmaps with
% the totalheatmaps function I had written. In addition to allowing you to
% keep track of each individual cell files, this program also asks to you
% the name the current dataset so you can keep track of the data set as
% whole as well and the files from each data set can be combined in a
% similar way using totalheatmaps script to generate combined heat.


% IMPORTANT EDIT Now this new code includes three major routines
%1> To show the number of points in each bin/cell along with the colormap
%  editor that
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
% This process will be repeated for each XY,YZ and XZ.



% Send me a message at dgvjay@illinois.edu if there are any questions.

% *******Reading in the input datafile
%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating EachIndividualCell Folder which stores the projection points (all
%the three planes) for each individual cell files in a folder called
%EACHINDIVIDUALCELL in the same
% working directory, It has sub folders called XY,YZ,XZ.


prompt = {'Enter the name that you want to give this dataset you just analyszed'};
dlg_title = 'Input';
num_lines = 1;
def = {'Enter any name that you like'};
datasetname = inputdlg(prompt,dlg_title,num_lines,def);
datasetname = datasetname{1};

mkdir('EachIndividualCell');  
fileattrib('./EachIndividualCell','+w');
cd('./EachIndividualCell');
mkdir('XY');
fileattrib('./XY','+w');
mkdir('YZ');
fileattrib('./YZ','+w');
mkdir('XZ');
fileattrib('./XZ','+w');
cd ..

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


listing = dir('*.txt');
listing=listing(1:end);
Ro_StormsXZPlane=[];
Ro_StormsYZPlane=[];
Ro_StormsXYPlane=[];

%%%%%%%%%%Defining DataSets that would append the projected points for each
%%%%%%%%%%individual cell in each three different planes.
xprojected_StormsXZPlane=[];
xprojected_StormsYZPlane=[];
xprojected_StormsXYPlane=[];
yprojected_StormsXZPlane=[];
yprojected_StormsYZPlane=[];
yprojected_StormsXYPlane=[];
zprojected_StormsXZPlane=[];
zprojected_StormsYZPlane=[];
zprojected_StormsXYPlane=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SelectedCells=0;
for j=1:numel(listing)
    
  dataset=listing(j,1).name;
  data=textread(dataset);
  data=data(:,[1 2 3]);
  x=data(:,1)./pixelingfactor;
  y=data(:,2)./pixelingfactor;
  z=data(:,3)./pixelingfactor;
       
  plot(x,y,'k.','LineWidth',0.4);
  axis equal;
  [axis_x axis_y]=ginput(2);    % Selecting two points along which we expect the helix axis to run 
  plot3(x,y,z,'k.','LineWidth',0.4);
  axis_z=[0 0] ;


  % ********Translate the points to make the first point of Helix axis the origin
  axis_x(1)
  axis_y(1)

  for i=1:numel(x)
    x=x-axis_x(1);
    y=y-axis_y(1);
    z=z;
    axis_x=axis_x-axis_x(1);
    axis_y=axis_y-axis_y(1);
    axis_z=[0 0] ;
  end


  % *********Rotate the axis so that the whole setup is now alinged on one of the primary y or x axis
  plot3(x,y,z,'k.','LineWidth',0.4);
  hold;
  plot3(axis_x,axis_y,axis_z,'r-','LineWidth',2);
  axis equal;
  theta=atan(axis_y(2)/axis_x(2));
  theta=-(pi/2-theta);
  x1=x*cos(theta) +  y*sin(theta);
  y1=y*cos(theta) -  x*sin(theta);
  z1=z;
  axis_x1=axis_x*cos(theta) +  axis_y*sin(theta);
  axis_y1=axis_y*cos(theta) -  axis_x*sin(theta);
  axis_z1=axis_z;
  plot3(x1,y1,z1,'r.');
  plot3(axis_x1,axis_y1,axis_z1,'k-','LineWidth',2);
  axis equal;
  hold;
  tic;pause(3);toc;


  % *********Arranging points in the order of their increasing Y ( as the helix axis is actually alinged along the Y axis)
  list_of_points=[];
  for i=1:numel(x1)
    distance_along_y_from_neworigin=abs(y1(i));
    list_of_points=[list_of_points; x1(i) y1(i) z1(i) distance_along_y_from_neworigin];
  end
  new_list=sortrows(list_of_points,4);
  x1=new_list(:,1);
  y1=new_list(:,2);
  z1=new_list(:,3);
  LengthofDistribution=max(y1)-min(y1);
  front_distance=(front_distancepercentage/100)*LengthofDistribution;
  back_distance=(back_distancepercentage/100)*LengthofDistribution;
  y11=[];
  x11=[];
  z11=[];
  Finaloutput=[];
  Originalpoints=numel(x1)
  Finaloutput=[Finaloutput Originalpoints];
  y1_abs=abs(y1);

  for i=1:numel(x1)
    if (y1_abs(i) >= front_distance)
        if (y1_abs(i) <= (max(y1_abs)-back_distance))
         x11=[x11 x1(i)];
         y11=[y11 y1(i)];
         z11=[z11 z1(i)];
        end
    end
  end
  Truncatedpoints=numel(x11)
  Finaloutput=[Finaloutput Truncatedpoints];
  x1=x11';
  y1=y11';
  z1=z11';
  points=[x1 y1 z1];
  pt1=[0 0 1];
  pt2=[0 0 0];
  pt3=[1 0 0];
  plane_StormsXZPlane=createPlane(pt1,pt2,pt3);
  projectedpoints_StormsXZPlane= projPointOnPlane(points, plane_StormsXZPlane);
  points=[x1 y1 z1];
  pt1=[0 1 0];
  pt2=[0 0 0];
  pt3=[1 0 0];
  plane_StormsXYPlane=createPlane(pt1,pt2,pt3);
  projectedpoints_StormsXYPlane= projPointOnPlane(points, plane_StormsXYPlane);
  pt1=[0 1 0];
  pt2=[0 0 0];
  pt3=[0 0 1];
  plane_StormsYZPlane=createPlane(pt1,pt2,pt3);
  projectedpoints_StormsYZPlane= projPointOnPlane(points, plane_StormsYZPlane);
  plot(projectedpoints_StormsXZPlane(:,1)',projectedpoints_StormsXZPlane(:,3)','bo');
  XZPoints=[projectedpoints_StormsXZPlane(:,1) projectedpoints_StormsXZPlane(:,2) projectedpoints_StormsXZPlane(:,3) ];
  title('Projections on all the planes.....decision coming up !');
  hold all;
  plot(projectedpoints_StormsXYPlane(:,1)',projectedpoints_StormsXYPlane(:,2)','ro');
  XYPoints=[projectedpoints_StormsXYPlane(:,1) projectedpoints_StormsXYPlane(:,2) projectedpoints_StormsXYPlane(:,3) ];
  plot(projectedpoints_StormsYZPlane(:,2)',projectedpoints_StormsYZPlane(:,3)','go');
  YZPoints=[projectedpoints_StormsYZPlane(:,1) projectedpoints_StormsYZPlane(:,2) projectedpoints_StormsYZPlane(:,3) ];
  temporaryname=datasetname;
  button = questdlg('Want to include this cell into your analysis ?','Inclusion Decision');

  if strcmp(button,'Yes')==1
     close all;    
     SelectedCells=SelectedCells+1;
     SelectedCells_String = sprintf('SelectedCell%d', SelectedCells);
     % Writing out the individual cell files now that you have indicated to
     % include it
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     projectedpoints_filenameXZ=sprintf('./EachIndividualCell/XZ/%s_%s_XZ.dat',temporaryname,SelectedCells_String);
     projectedpoints_filenameXY=sprintf('./EachIndividualCell/XY/%s_%s_XY.dat',temporaryname,SelectedCells_String);
     projectedpoints_filenameYZ=sprintf('./EachIndividualCell/YZ/%s_%s_YZ.dat',temporaryname,SelectedCells_String);
     dlmwrite(projectedpoints_filenameXZ,XZPoints,' ');
     dlmwrite(projectedpoints_filenameXY,XYPoints,' ');
     dlmwrite(projectedpoints_filenameYZ,YZPoints,' ');
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     % Appending the projected coordinates onto a bigger data structures
     % which would be used for the entire data set that you are analysing in
     % the current folder.
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     xprojected_StormsXZPlane=[xprojected_StormsXZPlane projectedpoints_StormsXZPlane(:,1)'];
     yprojected_StormsXZPlane=[yprojected_StormsXZPlane projectedpoints_StormsXZPlane(:,2)'];
     zprojected_StormsXZPlane=[zprojected_StormsXZPlane projectedpoints_StormsXZPlane(:,3)'];
     xprojected_StormsYZPlane=[xprojected_StormsYZPlane projectedpoints_StormsYZPlane(:,1)'];
     yprojected_StormsYZPlane=[yprojected_StormsYZPlane projectedpoints_StormsYZPlane(:,2)'];
     zprojected_StormsYZPlane=[zprojected_StormsYZPlane projectedpoints_StormsYZPlane(:,3)'];
     xprojected_StormsXYPlane=[xprojected_StormsXYPlane projectedpoints_StormsXYPlane(:,1)'];
     yprojected_StormsXYPlane=[yprojected_StormsXYPlane projectedpoints_StormsXYPlane(:,2)'];
     zprojected_StormsXYPlane=[zprojected_StormsXYPlane projectedpoints_StormsXYPlane(:,3)'];
   end

 if strcmp(button,'No')==1
      close all;
 end

 if strcmp(button,'Cancel')==1
      msgbox('You have decided to cancel the entire operation hence no further analysis generated...closing the entire program in 2 seconds','Cancel Decision Taken');
      tic;pause(2);toc;
      close all;
      clear all;
      return;
 end

end

if isempty(xprojected_StormsXZPlane)==1
    msgbox('You have selected 0 cells hence no further analysis generated...closing the entire program in 2 seconds','Final Selection');
    tic;pause(2);toc;
    close all;
    clear all;
    return;
end

if isempty(xprojected_StormsXZPlane)==0
    msgscript=sprintf('You have selected %d/%d cells',SelectedCells,numel(listing)); 
    msgbox(msgscript,'Final Selection'); 
    
    
    %Distance Calculation...Right now all of them being calculated as
    %Euclidean distance from the center for each projection because they
    %are all centered, Please make changes in any distance criteria here
    
    Ro_StormsXZPlane=[((xprojected_StormsXZPlane.^2 +zprojected_StormsXZPlane.^2).^(0.5))];
    Ro_StormsYZPlane=[((xprojected_StormsYZPlane.^2 +zprojected_StormsYZPlane.^2).^(0.5))];
    Ro_StormsXYPlane=[((xprojected_StormsXYPlane.^2 +zprojected_StormsXYPlane.^2).^(0.5))];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    HistogramEdges=0:histogrambinsize:max(Ro_StormsXZPlane);   % Here I am creating a series of edges for histc from minimum value of Radius(min(Ro) all upto the maximum value of Radius(max(Ro));
    HistogramValues=[HistogramEdges' histc(Ro_StormsXZPlane,HistogramEdges)'];
    HistogramValues_StormsXZPlane=strcat(datasetname,'_HistogramValues_StormsXZPlane.dat');
    dlmwrite(HistogramValues_StormsXZPlane,HistogramValues,' ');
    projectedpoints_write=[xprojected_StormsXZPlane' yprojected_StormsXZPlane' zprojected_StormsXZPlane'];
    projectedpoints_write_StormsXZPlane=strcat(datasetname,'_projectedpoints_write_StormsXZPlane.dat');
    dlmwrite(projectedpoints_write_StormsXZPlane,projectedpoints_write,' ');
    imagewindow = figure('visible','on');
    plot(xprojected_StormsXZPlane,zprojected_StormsXZPlane,'bo');
    projectpoints_on_StormsXZplane=strcat(datasetname,'_projectpoints_on_StormsXZplane.jpg');
    saveas(imagewindow,projectpoints_on_StormsXZplane);
    [n,c] = hist3([xprojected_StormsXZPlane' zprojected_StormsXZPlane'],[heatmapbinx,heatmapbiny]);
    close all;
    
    if strcmp(heatmapanalysis_on_or_off,'on')==1
    heatmaptext(n);
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
    xlabel('X Axis');
    ylabel('Z Axis');
    colorbar;
    new_fig=gcf;
    set(new_fig,'Colormap',mycmap)
    colormapeditor
    tic;pause(timelapse);toc;
    end
    
    if strcmp(heatmapanalysis_on_or_off,'on')==0
    imagewindow = figure('visible','on');
    heatmaptext(n);
    figure ;
    imagesc(c{1},c{2},n);
    xlabel('X Axis');
    ylabel('Z Axis');
    end
    
    projectpoints_on_StormsXZplaneHeatMap=strcat(datasetname,'_projectpoints_on_StormsXZplane_HeatMap.jpg');
    saveas(imagewindow,projectpoints_on_StormsXZplaneHeatMap);
    projectpoints_on_StormsXZplaneHeatMap=strcat(datasetname,'_projectpoints_on_StormsXZplane_HeatMap.fig');
    saveas(imagewindow,projectpoints_on_StormsXZplaneHeatMap);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    HistogramEdges=0:histogrambinsize:max(Ro_StormsYZPlane);   % Here I am creating a series of edges for histc from minimum value of Radius(min(Ro) all upto the maximum value of Radius(max(Ro));
    HistogramValues=[HistogramEdges' histc(Ro_StormsYZPlane,HistogramEdges)'];
    HistogramValues_StormsYZPlane=strcat(datasetname,'_HistogramValues_StormsYZPlane.dat');
    dlmwrite(HistogramValues_StormsYZPlane,HistogramValues,' ');
    projectedpoints_write=[xprojected_StormsYZPlane' yprojected_StormsYZPlane' zprojected_StormsYZPlane'];
    projectedpoints_write_StormsYZPlane=strcat(datasetname,'_projectedpoints_write_StormsYZPlane.dat');
    dlmwrite(projectedpoints_write_StormsYZPlane,projectedpoints_write,' ');
    imagewindow = figure('visible','on');
    plot(yprojected_StormsYZPlane,zprojected_StormsYZPlane,'ro');
    projectpoints_on_StormsYZplane=strcat(datasetname,'_projectpoints_on_StormsYZplane.jpg'); 
    saveas(imagewindow,projectpoints_on_StormsYZplane);
    [n,c] = hist3([yprojected_StormsYZPlane' zprojected_StormsYZPlane'],[heatmapbinx,heatmapbiny]);   
    close all;
    if strcmp(heatmapanalysis_on_or_off,'on')==1
    heatmaptext(n);
    colormapeditor
    h = gcf;
    tic;pause(timelapse);toc;
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
    xlabel('Y Axis');
    ylabel('Z Axis');
    colorbar;
    new_fig=gcf;
    set(new_fig,'Colormap',mycmap)
    colormapeditor
    tic;pause(timelapse);toc;
    end
    if strcmp(heatmapanalysis_on_or_off,'on')==0
    imagewindow = figure('visible','on');
    heatmaptext(n);
    figure ;
    imagesc(c{1},c{2},n);
    xlabel('Y Axis');
    ylabel('Z Axis');
    end
    projectpoints_on_StormsYZplaneHeatMap=strcat(datasetname,'_projectpoints_on_StormsYZplane_HeatMap.jpg');
    saveas(imagewindow,projectpoints_on_StormsYZplaneHeatMap);
    projectpoints_on_StormsYZplaneHeatMap=strcat(datasetname,'_projectpoints_on_StormsYZplane_HeatMap.fig');
    saveas(imagewindow,projectpoints_on_StormsYZplaneHeatMap);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    HistogramEdges=0:histogrambinsize:max(Ro_StormsXYPlane);   % Here I am creating a series of edges for histc from minimum value of Radius(min(Ro) all upto the maximum value of Radius(max(Ro));
    HistogramValues=[HistogramEdges' histc(Ro_StormsXYPlane,HistogramEdges)'];
    HistogramValues_StormsXYPlane=strcat(datasetname,'_HistogramValues_StormsXYPlane.dat');
    dlmwrite(HistogramValues_StormsXYPlane,HistogramValues,' ');
    projectedpoints_write=[xprojected_StormsXYPlane' yprojected_StormsXYPlane' zprojected_StormsXYPlane'];
    projectedpoints_write_StormsXYPlane=strcat(datasetname,'_projectedpoints_write_StormsXYPlane.dat');
    dlmwrite(projectedpoints_write_StormsXYPlane,projectedpoints_write,' ');
    imagewindow = figure('visible','on'); 
    plot(xprojected_StormsXYPlane,yprojected_StormsXYPlane,'go');
    projectpoints_on_StormsXYplane=strcat(datasetname,'_projectpoints_on_StormsXYplane.jpg'); 
    saveas(imagewindow,projectpoints_on_StormsXYplane);
    [n,c] = hist3([xprojected_StormsXYPlane' yprojected_StormsXYPlane'],[heatmapbinx,heatmapbiny]);
    close all;
    
    if strcmp(heatmapanalysis_on_or_off,'on')==1
    heatmaptext(n);
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
    xlabel('X Axis');
    ylabel('Y Axis');
    colorbar;
    new_fig=gcf;
    set(new_fig,'Colormap',mycmap)
    colormapeditor   
    tic;pause(timelapse);toc;
    end
    
    if strcmp(heatmapanalysis_on_or_off,'on')==0
    imagewindow = figure('visible','on');
    heatmaptext(n);
    figure;
    imagesc(c{1},c{2},n);
    xlabel('X Axis');
    ylabel('Y Axis');
    end
    
    projectpoints_on_StormsXYplaneHeatMap=strcat(datasetname,'_projectpoints_on_StormsXYplane_HeatMap.jpg');
    saveas(imagewindow,projectpoints_on_StormsXYplaneHeatMap);
    projectpoints_on_StormsXYplaneHeatMap=strcat(datasetname,'_projectpoints_on_StormsXYplane_HeatMap.fig');
    saveas(imagewindow,projectpoints_on_StormsXYplaneHeatMap);
end
close all;
clear all;
