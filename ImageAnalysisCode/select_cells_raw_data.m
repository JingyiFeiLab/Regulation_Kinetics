function select_cells_raw_data(coordinatefile,Columnx,Columny,Columnz,ColorLineColumn,referenceimage,pixelscaling)
% A script for selecting the points corresponding to cells in the reference image.
% please send me a message at dgvjay@illinois.edu if there is any query about these scripts- Digvijay Singh

% *****PARAMETERS EXPLAINED
%coordinatefile                                  The input data file name which has x y columns plus a column containing colorcode of the points
%Columnx                                         The column from the input file which has to be read as X coordinates
%Columny                                         The column from the input file which has to be read as Y coordinates
%Columnz                                         The column from the input file which has to be read as Z coordinate
%Column_No_Spots                                 The column from the input file which has to be read as Number of Spots
%ColumnR                                         The column from the input file which has to be read as AverageR
%ColorLineColumn                                 The Column Contains the ColorCode of the points
%referenceimage                                  The reference image file name on which x,y coordinates have to be overlapped
%pixelscaling                                    Pixel Scaling is to convert coordinates in nm,Angstrom to Pixel Coordinates


%******EXAMPLE OF USAGE
% select_cells('sample.txt',1,2,4,'Frame14.jpg',100);
% Where 
% sample.txt     == is the input file containing points and their colorcode in 4th Column
% 1              == means X coordinates is in column1
% 2              == means Y coordinates is in column2
% 4              == means ColorCoordinates is in column2
% Frame14.jpg    == Image on which overlap has to be done
% 100            == The sample input file has nm coordinates which when divided by 100 yield pixel coordinates

% NOTE: while running and if you close the figure file/stop the script before selecting the mentioned number of cell, 
% it will show an error but ignore it, all your cells files are stored in the same directoy


if exist('newxy.dat','file')==2
    delete('newxy.dat');
end

particles=textread(coordinatefile);              % Read in the input file
x=particles(:,Columnx);                          % Select the column containing X coordinates in Pixel form
y=particles(:,Columny);                          % Select the column containing Y coordinates in Pixel form
z=particles(:,Columnz);                          % Select the column containing Y coordinates in Pixel form
%No_Spots=particles(:,Column_No_Spots);                  % Select the column containing Number of Spots for each cluster center
%AverageR=particles(:,ColumnR);                          % Select the column containing Average R for each cluster center
Color=particles(:,ColorLineColumn);              % Select the column containing ColorNumbers for each point

%*********Scaling the X Y Coordinates as per the pixel of the input image
image = imread(referenceimage);                  % Put the reference image based on which you want to select the points
X = size(image,1) ;                              % CapitalX is the number of X pixels in the reference image
Y = size(image,2) ;                              % CapitalY is the number of Y pixels in the reference image
scalex=X/256;                                    % Getting the scaling factors to match X and Y pixels from the image to the ones from the coordinate file
scaley=Y/256;                                    % We are assuming that the reference image is 256 by 256 in pixel
x=(scalex.*x)./pixelscaling;                     % Scale it up...Pixel Scaling is to convert coordinates in nm,Angstrom to Pixel Coordinates
y=(scaley.*y)./pixelscaling;                     % Scale it up...Pixel Scaling is to convert coordinates in nm,Angstrom to Pixel Coordinates
z=(scaley.*z)./pixelscaling;                     % Scale it up...Pixel Scaling is to convert coordinates in nm,Angstrom to Pixel Coordinates


overlap(referenceimage,x,y,Color);
continuetranslation=1;
while continuetranslation==1
prompt={'X_Translation:','Y_Translation:','Continue the translation(1 == Yes && 0== NO) ?'};   % A box will take in the values for the X/Ytranslation
title='Translation';                             % The title of the box
answer=inputdlg(prompt,title);
Xtranslation = str2num(answer{1}); 
Ytranslation = str2num(answer{2});
continuetranslation = str2num(answer{3});
x=x-Xtranslation;                                % Translated
y=y-Ytranslation;
overlap(referenceimage,x,y,Color);
end
points=[x y z];
filename=sprintf('newxy.dat');                   % The Newly Translated Points are saved as newxy.dat
dlmwrite(filename,points);                       % Saving the coordinates
save('referenceimage.mat','referenceimage');
save('Color.mat','Color');

                                       
% NOW THE SELECTION OF POINTS WITHIN ANY SELECTED POLYGONAL AREA
%clear 
load('newxy.dat');
load('referenceimage.mat');
load('Color.mat');
x=newxy(:,1);
y=newxy(:,2);
overlap(referenceimage,x,y,Color);
prompt={'No. of Cells to select:'};              % A box will take in the value for the number of cells that you want to select
title='Cell Counts';                             % The title of the box
answer=inputdlg(prompt,title);
cellcount = str2num(answer{1});                  % Mention the number of cells that you want to select from the given image
count=0;
%CellAnalysis=[];
while count < cellcount                          % Allows you to select upto "cellcount" cells in the image through polygonal selections of the area
    count=count+1;
   % Total_Numberof_Spots_1=0;
   % Total_Numberof_Spots_2=0;
   % Total_Numberof_Spots_3=0;
   % Total_Numberof_Spots_4=0;
   % Total_Numberof_Spots_5=0;
    h = impoly();
    nodes = getPosition(h);
    Xpoly=nodes(:,1);
    Ypoly=nodes(:,2);
    Area = polyarea(Xpoly,Ypoly);
    selected_indices = inpoly(newxy(:, 1:2),nodes);      % Using script inpoly.m ( From MATLAB CENTRAL)
    temporarypoints=[]; 
    ColumnCell=[];
    for i=1:numel(selected_indices)
      if selected_indices(i)==1
        temporarypoints=[temporarypoints; [newxy(i,1) newxy(i,2) newxy(i,3) Color(i)]];  % Making a list of the points within the selected area
        ColumnCell=[ColumnCell Color(i)];
        
       % if Color(i)==1
       % Total_Numberof_Spots_1=Total_Numberof_Spots_1+newxy(i,5);
       % end
        
       % if Color(i)==2
       % Total_Numberof_Spots_2=Total_Numberof_Spots_2+newxy(i,5);
       % end
        
       % if Color(i)==3
       % Total_Numberof_Spots_3=Total_Numberof_Spots_3+newxy(i,5);
       % end
        
       % if Color(i)==4
       % Total_Numberof_Spots_4=Total_Numberof_Spots_4+newxy(i,5);
       % end
        
       % if Color(i)==5
       % Total_Numberof_Spots_5=Total_Numberof_Spots_5+newxy(i,5);
       % end
        
      end
    end
% CellAnalysis=[CellAnalysis; [count length(find(ColumnCell==1)) length(find(ColumnCell==2)) length(find(ColumnCell==3)) length(find(ColumnCell==4)) ...
%    length(find(ColumnCell==5)) Total_Numberof_Spots_1 Total_Numberof_Spots_2 Total_Numberof_Spots_3 Total_Numberof_Spots_4 polyarea Total_Numberof_Spots_5 Area]]; %check the color number
filename=sprintf('cell_%d.txt',count);
dlmwrite(filename,temporarypoints,' ');              % Saving the coordinates with X,Y coordinates of the points and their corresponding Number of Spots,Average Radius and the ColorCode.

%filename2=sprintf('TotalCells_Analysis.dat');
%dlmwrite(filename2,CellAnalysis,' ');              % Saving the TotalCellAnalysis


end
delete('referenceimage.mat');                    % Deleting some *.mat files which were created during the execution
delete('Color.mat');
close all;