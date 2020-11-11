function projections(pixelingfactor,histogrambinsize,heatmapbinx,heatmapbiny,timelapse,heatmapanalysis_on_or_off)
% Example of Usage : projections(1,30,30,30,10,'off');
% This function will read all the .txt files in the folder and project
% their points and overlay them It can also include the fourth column for
% the color codes.

% pixelingfactor :1, If you want to convert between nm/Angstrom and pixel
% cooordinates, then you can mention the appropriate pixelingfactor.

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

% front_distancepercentage :(20 indicates 20% of the points from the front)
% The initial few points may be discarded..so you can mention the distance
% that needs to be truncated from the beginning. 
% back_distancepercentage% :(20 indicates 20% of the points from the front) 
% The initial few points may be discarded..so you can mention the distance  
% that needs to be truncated from the end. 

% prompt = {'Enter the FRONT distance percentage that you want to truncate for the points for projection on XY Plane'};
% dlg_title = 'FRONT DISTANCE PERCENTAGE TRUNCATION_XY';
% num_lines = 1;
% def = {'FRONT DISTANCE PERCENTAGE TRUNCATION_XY'};
% front_distancepercentage_XY = inputdlg(prompt,dlg_title,num_lines,def);
% front_distancepercentage_XY = front_distancepercentage_XY{1};
% front_distancepercentage_XY = str2num(front_distancepercentage_XY);
%     
% prompt = {'Enter the END distance percentage that you want to truncate for the points for projection on XY Plane'};
% dlg_title = 'END DISTANCE PERCENTAGE TRUNCATION_XY';
% num_lines = 1;
% def = {'END DISTANCE PERCENTAGE TRUNCATION_XY'};
% back_distancepercentage_XY = inputdlg(prompt,dlg_title,num_lines,def);
% back_distancepercentage_XY = back_distancepercentage_XY{1};
% back_distancepercentage_XY = str2num(back_distancepercentage_XY);


% prompt = {'Enter the FRONT distance percentage that you want to truncate for the points for projection on YZ Plane'};
% dlg_title = 'FRONT DISTANCE PERCENTAGE TRUNCATION_YZ';
% num_lines = 1;
% def = {'FRONT DISTANCE PERCENTAGE TRUNCATION_YZ'};
% front_distancepercentage_YZ = inputdlg(prompt,dlg_title,num_lines,def);
% front_distancepercentage_YZ = front_distancepercentage_YZ{1};
% front_distancepercentage_YZ = str2num(front_distancepercentage_YZ);
%     
% prompt = {'Enter the END distance percentage that you want to truncate for the points for projection on YZ Plane'};
% dlg_title = 'END DISTANCE PERCENTAGE TRUNCATION_YZ';
% num_lines = 1;
% def = {'END DISTANCE PERCENTAGE TRUNCATION_YZ'};
% back_distancepercentage_YZ = inputdlg(prompt,dlg_title,num_lines,def);
% back_distancepercentage_YZ = back_distancepercentage_YZ{1};
% back_distancepercentage_YZ = str2num(back_distancepercentage_YZ);
% 
% 
% 
prompt = {'Enter the FRONT distance percentage that you want to truncate for the points for projection on XZ Plane'};
dlg_title = 'FRONT DISTANCE PERCENTAGE TRUNCATION_XZ';
num_lines = 1;
def = {'FRONT DISTANCE PERCENTAGE TRUNCATION_XZ'};
front_distancepercentage_XZ = inputdlg(prompt,dlg_title,num_lines,def);
front_distancepercentage_XZ = front_distancepercentage_XZ{1};
front_distancepercentage_XZ = str2num(front_distancepercentage_XZ);
    
prompt = {'Enter the END distance percentage that you want to truncate for the points for projection on XZ Plane'};
dlg_title = 'END DISTANCE PERCENTAGE TRUNCATION_XZ';
num_lines = 1;
def = {'END DISTANCE PERCENTAGE TRUNCATION_XZ'};
back_distancepercentage_XZ = inputdlg(prompt,dlg_title,num_lines,def);
back_distancepercentage_XZ = back_distancepercentage_XZ{1};
back_distancepercentage_XZ = str2num(back_distancepercentage_XZ);
    
  front_distancepercentage_XY=0;
  back_distancepercentage_XY=0;
  front_distancepercentage_YZ=0;
  back_distancepercentage_YZ=0;
  

mkdir('EachIndividualCell');  
fileattrib('./EachIndividualCell','+w');
cd('./EachIndividualCell');
mkdir('XY');
fileattrib('./XY','+w');
mkdir('YZ');
fileattrib('./YZ','+w');
mkdir('XZ');
fileattrib('./XZ','+w');
mkdir('YZPositive');
fileattrib('./YZ','+w');
mkdir('YZNegative');
fileattrib('./YZ','+w');
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
xprojected_StormsYZPlane_NegativeHalf= [];
yprojected_StormsYZPlane_NegativeHalf= [];
zprojected_StormsYZPlane_NegativeHalf= [];
xprojected_StormsYZPlane_PositiveHalf= [];
yprojected_StormsYZPlane_PositiveHalf= [];
zprojected_StormsYZPlane_PositiveHalf= [];

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
  %y1=y1- (max(y1)-min(y1))/2;
  [YS IS]=sort(y1);
  YS1=y1(IS);
  YS1=YS1(1:100);
  [YS IS]=sort(y1,'descend');
  YS2=y1(IS);
  YS2=YS2(1:100);
  YStotal=[YS1;YS2];
  
  y1dash=y1-mean(YStotal);
  y1=y1dash;
  z1=new_list(:,3);
  xoriginal=x1;
  yoriginal=y1;
  zoriginal=z1;
  LengthofDistribution=max(y1)-min(y1);
     
  %%%%%%%%%%%% XY Set of Points Truncation according to their respective
  %%%%%%%%%%%% front and end percentages
  front_distance_XY=(front_distancepercentage_XY/100)*LengthofDistribution;
  back_distance_XY=(back_distancepercentage_XY/100)*LengthofDistribution;
  y11_XY=[];
  x11_XY=[];
  z11_XY=[];
  Originalpoints_XY=numel(xoriginal)
  y1_abs_XY=abs(yoriginal);

  for i=1:numel(xoriginal)
    if (y1_abs_XY(i) >=   front_distance_XY)
        if (y1_abs_XY(i) <= (max(y1_abs_XY)-back_distance_XY))
         x11_XY=[x11_XY xoriginal(i)];
         y11_XY=[y11_XY yoriginal(i)];
         z11_XY=[z11_XY zoriginal(i)];
        end
    end
  end
  Truncatedpoints_XY=numel(x11_XY)
  points_XY=[x11_XY' y11_XY' z11_XY'];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%%%%%%% XZ Set of Points Truncation according to their respective
  %%%%%%%%%%%% front and end percentages
  front_distance_XZ=(front_distancepercentage_XZ/100)*LengthofDistribution;
  back_distance_XZ=(back_distancepercentage_XZ/100)*LengthofDistribution;
  y11_XZ=[];
  x11_XZ=[];
  z11_XZ=[];
  Originalpoints_XZ=numel(xoriginal)
  y1_abs_XZ=abs(yoriginal);
  for i=1:numel(xoriginal)
    if (y1_abs_XZ(i) >=   front_distance_XZ)
        if (y1_abs_XZ(i) <= (max(y1_abs_XZ)-back_distance_XZ))
         x11_XZ=[x11_XZ xoriginal(i)];
         y11_XZ=[y11_XZ yoriginal(i)];
         z11_XZ=[z11_XZ zoriginal(i)];
        end
    end
  end
  Truncatedpoints_XZ=numel(x11_XZ)
  points_XZ=[x11_XZ' y11_XZ' z11_XZ'];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%%%%%%% YZ Set of Points Truncation according to their respective
  %%%%%%%%%%%% front and end percentages
  front_distance_YZ=(front_distancepercentage_YZ/100)*LengthofDistribution;
  back_distance_YZ=(back_distancepercentage_YZ/100)*LengthofDistribution;
  y11_YZ=[];
  x11_YZ=[];
  z11_YZ=[];
  Originalpoints_YZ=numel(xoriginal)
  y1_abs_YZ=abs(yoriginal);

  for i=1:numel(xoriginal)
    if (y1_abs_YZ(i) >=   front_distance_YZ)
        if (y1_abs_YZ(i) <= (max(y1_abs_YZ)-back_distance_YZ))
         x11_YZ=[x11_YZ xoriginal(i)];
         y11_YZ=[y11_YZ yoriginal(i)];
         z11_YZ=[z11_YZ zoriginal(i)];
        end
    end
  end
  Truncatedpoints_YZ=numel(x11_YZ)
  points_YZ=[x11_YZ' y11_YZ' z11_YZ'];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    for gg=1:numel(xprojected_StormsXZPlane)
   points_YZ_PositiveHalf=points_YZ(points_YZ(:,1)>0,:);
   points_YZ_NegativeHalf=points_YZ(points_YZ(:,1)<=0,:);
   points_XZ_PositiveHalf=points_XZ(points_XZ(:,1)>0,:);
   points_XZ_NegativeHalf=points_XZ(points_XZ(:,1)<=0,:);
   points_XY_PositiveHalf=points_XY(points_XY(:,1)>0,:);
   points_XY_NegativeHalf=points_XY(points_XY(:,1)<=0,:);
   
  pt1=[0 0 1];
  pt2=[0 0 0];
  pt3=[1 0 0];
  plane_StormsXZPlane=createPlane(pt1,pt2,pt3);
  projectedpoints_StormsXZPlane= projPointOnPlane(points_XZ, plane_StormsXZPlane);
  projectedpoints_StormsXZPlane_NegativeHalf= projPointOnPlane(points_XZ_NegativeHalf, plane_StormsXZPlane);
  projectedpoints_StormsXZPlane_PositiveHalf= projPointOnPlane(points_XZ_PositiveHalf, plane_StormsXZPlane);
  points=[x1 y1 z1];
  pt1=[0 1 0];
  pt2=[0 0 0];
  pt3=[1 0 0];
  plane_StormsXYPlane=createPlane(pt1,pt2,pt3);
  projectedpoints_StormsXYPlane= projPointOnPlane(points_XY, plane_StormsXYPlane);
  projectedpoints_StormsXYPlane_NegativeHalf= projPointOnPlane(points_XY_NegativeHalf, plane_StormsXYPlane);
  projectedpoints_StormsXYPlane_PositiveHalf= projPointOnPlane(points_XY_PositiveHalf, plane_StormsXYPlane);
  pt1=[0 1 0];
  pt2=[0 0 0];
  pt3=[0 0 1];
  plane_StormsYZPlane=createPlane(pt1,pt2,pt3);
  projectedpoints_StormsYZPlane_NegativeHalf= projPointOnPlane(points_YZ_NegativeHalf, plane_StormsYZPlane);
  projectedpoints_StormsYZPlane_PositiveHalf= projPointOnPlane(points_YZ_PositiveHalf, plane_StormsYZPlane);
  projectedpoints_StormsYZPlane= projPointOnPlane(points_YZ, plane_StormsYZPlane);
  plot(projectedpoints_StormsXZPlane(:,1)',projectedpoints_StormsXZPlane(:,3)','bo');
  XZPoints=[projectedpoints_StormsXZPlane(:,1) projectedpoints_StormsXZPlane(:,2) projectedpoints_StormsXZPlane(:,3) ];
  title('Projections on all the planes.....decision coming up !');
  hold all;
  %plot((projectedpoints_StormsXYPlane(:,1)-5*(max(projectedpoints_StormsXYPlane(:,1))-min(projectedpoints_StormsXYPlane(:,1))))',(projectedpoints_StormsXYPlane(:,2)-5*(max(projectedpoints_StormsXYPlane(:,2))-min(projectedpoints_StormsXYPlane(:,2))))','ro');
  XYPoints=[projectedpoints_StormsXYPlane(:,1) projectedpoints_StormsXYPlane(:,2) projectedpoints_StormsXYPlane(:,3) ];
  %plot((projectedpoints_StormsYZPlane(:,2)+5*(max(projectedpoints_StormsYZPlane(:,2))-min(projectedpoints_StormsYZPlane(:,2))))',(projectedpoints_StormsYZPlane(:,3)+5*(max(projectedpoints_StormsYZPlane(:,2))-min(projectedpoints_StormsYZPlane(:,2))))','go');
  YZPoints=[projectedpoints_StormsYZPlane(:,1) projectedpoints_StormsYZPlane(:,2) projectedpoints_StormsYZPlane(:,3) ];
  YZPointsPositive=[projectedpoints_StormsYZPlane_PositiveHalf(:,1) projectedpoints_StormsYZPlane_PositiveHalf(:,2) projectedpoints_StormsYZPlane_PositiveHalf(:,3) ];
  YZPointsNegative=[projectedpoints_StormsYZPlane_NegativeHalf(:,1) projectedpoints_StormsYZPlane_NegativeHalf(:,2) projectedpoints_StormsYZPlane_NegativeHalf(:,3) ];
  temporaryname=datasetname;
  button = questdlg('Want to include this cell into your analysis ?','Inclusion Decision');

  if strcmp(button,'Yes')==1
     close all;    
     % Writing out the individual cell files now that you have indicated to
     % include it
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     projectedpoints_filenameXZ=sprintf('./EachIndividualCell/XZ/%s_%s_XZ.dat',temporaryname,dataset);
     projectedpoints_filenameXY=sprintf('./EachIndividualCell/XY/%s_%s_XY.dat',temporaryname,dataset);
     projectedpoints_filenameYZ=sprintf('./EachIndividualCell/YZ/%s_%s_YZ.dat',temporaryname,dataset);
     projectedpoints_filenameYZPositive=sprintf('./EachIndividualCell/YZPositive/%s_%s_YZ_Positive.dat',temporaryname,dataset);
     projectedpoints_filenameYZNegative=sprintf('./EachIndividualCell/YZNegative/%s_%s_YZ_Negative.dat',temporaryname,dataset);
     dlmwrite(projectedpoints_filenameXZ,XZPoints,' ');
     dlmwrite(projectedpoints_filenameXY,XYPoints,' ');
     dlmwrite(projectedpoints_filenameYZ,YZPoints,' ');
     dlmwrite(projectedpoints_filenameYZPositive,YZPointsPositive,' ');
     dlmwrite(projectedpoints_filenameYZNegative,YZPointsNegative,' ');
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
     xprojected_StormsYZPlane_NegativeHalf= [xprojected_StormsYZPlane_NegativeHalf projectedpoints_StormsYZPlane_NegativeHalf(:,1)'];
     yprojected_StormsYZPlane_NegativeHalf= [yprojected_StormsYZPlane_NegativeHalf projectedpoints_StormsYZPlane_NegativeHalf(:,2)'];
     zprojected_StormsYZPlane_NegativeHalf= [zprojected_StormsYZPlane_NegativeHalf projectedpoints_StormsYZPlane_NegativeHalf(:,3)'];
     xprojected_StormsYZPlane_PositiveHalf= [xprojected_StormsYZPlane_PositiveHalf projectedpoints_StormsYZPlane_PositiveHalf(:,1)'];
     yprojected_StormsYZPlane_PositiveHalf= [yprojected_StormsYZPlane_PositiveHalf projectedpoints_StormsYZPlane_PositiveHalf(:,2)'];
     zprojected_StormsYZPlane_PositiveHalf= [zprojected_StormsYZPlane_PositiveHalf projectedpoints_StormsYZPlane_PositiveHalf(:,3)'];
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
end