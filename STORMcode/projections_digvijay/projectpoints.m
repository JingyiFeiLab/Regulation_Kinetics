function projectpoints(dataset,pixelingfactor,front_distancepercentage,back_distancepercentage,Window)

% Example of Usage :
% fit_to_helix('110712_DJ480_10min_1_noTrail_1.txt',1,0.2,0.2,100);

% dataset:110712_DJ480_10min_1_noTrail_1.txt 
% This is the file containing the x y z coordinates of the points which are to be fit to a helix.
% It can also include the fourth column for the color codes.

% pixelingfactor :1, If you want to convert between nm/Angstrom and pixel cooordinates, then you can mention the appropriate pixelingfactor.

% front_distance : The initial few points may be discarded..so you can mention the distance  that needs to be truncated from the beginning.
% end_distance :The initial few points may be discarded..so you can mention the distance  that needs to be truncated from the end.
% end : 0.2 The final few points may be discarded..so you can mention the fraction of points that needs to be truncated from the end.
% Window : Small Slabs on the Unwrapped Points Plot on which you want to
% fit straight lines, The idea is to play around with different window to
% get the distribution that makes closest sense to a helical distribution.


% *******Reading in the input datafile
data=textread(dataset);
data=data(:,[1 2 3]);
x=data(:,1)./pixelingfactor;
y=data(:,2)./pixelingfactor;
z=data(:,3)./pixelingfactor;

% front_number=ceil(front*numel(x));
% back_number=ceil(back*numel(x));
% 
% x=x(front_number:(numel(x)-back_number));
% y=y(front_number:(numel(y)-back_number));
% z=z(front_number:(numel(z)-back_number));

       
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
% front_number=ceil(front*numel(x1));
% back_number=ceil(back*numel(x1));
% 
% x1=x1(front_number:(numel(x1)-back_number));
% y1=y1(front_number:(numel(y1)-back_number));
% z1=z1(front_number:(numel(z1)-back_number));

% *********Projecting all the points on the XZ plane passing through origin and we will take the average of their Radius values from the origin which would
% be sort of an expected Radius of the Helical distribution

points=[x1 y1 z1];
pt1=[0 0 1];
pt2=[0 0 0];
pt3=[1 0 0];
plane=createPlane(pt1,pt2,pt3);
projectedpoints= projPointOnPlane(points, plane);
save('projectedpoints.mat','projectedpoints');
xprojected=projectedpoints(:,1);
yprojected=projectedpoints(:,2);
zprojected=projectedpoints(:,3);

Ro=0;
for i=1:numel(xprojected)
 Ro=Ro+sqrt( xprojected(i).^2+ zprojected(i).^2);
     
end
imagewindow = figure('visible','on');
plot(xprojected,zprojected,'bo');
saveas(imagewindow,'ProjectedPoints.jpg');

disp('Average Expected Radius of the Helix ::'); Ro=Ro./numel(xprojected) 
Finaloutput=[Finaloutput Ro];
longitudional_dist=[];
circumferential_dist=[];
RcList=[];
for i=1:numel(x1)
    rc=sqrt(x1(i).^2+z1(i).^2);
    RcList=[RcList rc];
    
    if rc <1.1*Ro && rc > 0.9*Ro   % Only taking in the points which are within 20% region of average Radius calculated above (Ro)....
                                   % This is to avoid including in points which are otherwise wayoff from the helical distribution
    if x1(i)==-0
        x1(i)=0;
    end
    if z1(i)==-0
        z1(i)=0;
    end
    
    if x1(i)>=0 && z1(i) >=0
        angle=asin(abs(x1(i))/rc);
    end
    if x1(i)>=0 && z1(i) < 0
        angle=pi/2 + acos(abs(x1(i))/rc);
    end
    if x1(i)<0 && z1(i) < 0
        angle=pi+ asin(abs(x1(i))/rc);
    end
    if x1(i)<0 && z1(i) >= 0
        angle=3*pi/2 + acos(abs(x1(i))/rc);
    end  
    cdist=rc*angle;
    circumferential_dist=[circumferential_dist cdist];
    longitudional_dist=[longitudional_dist abs(y1(i))];
    end
end
Rclist_Line=RcList';
Rcname=sprintf('RcList_%s.dat',dataset);

dlmwrite(Rcname,Rclist_Line,' ');
imagewindow = figure('visible','off');
edges = 5:10:round(max(Rclist_Line));
[Xhist Yhist]=hist(Rclist_Line,edges);
Histout=[Yhist' Xhist'];
dlmwrite('Histogram.txt',Histout, ' ' );
saveas(imagewindow,'R_Histogram.jpg');

% imagewindow = figure('visible','on');
Rclist=sort(RcList,'descend');
Rclist_Top100=Rclist(1:100);
Rmax=sum(Rclist_Top100)/100;
Finaloutput=[Finaloutput Rmax];
Finaloutput=Finaloutput';
dlmwrite('Finaloutput.txt',Finaloutput);
% longitudional_dist=longitudional_dist';
% circumferential_dist=circumferential_dist';
% plot(longitudional_dist,circumferential_dist,'b.');
% tic;pause(5);toc;
% dfunctionvalues=[longitudional_dist circumferential_dist];
% save('dfunctionvalues.mat','dfunctionvalues');
% 
% 
% %Rolling Function
% NWindows=max(longitudional_dist)./Window
% NWindows=round(NWindows);
% close all;
% imagewindow = figure('visible','on');
% for j=1:NWindows
%     Xset=[];
%     Yset=[];
%     for i=1:numel(longitudional_dist)
%         if (longitudional_dist(i)./Window <j && longitudional_dist(i)./Window >=j-1)
%             Xset=[Xset longitudional_dist(i)];
%             Yset=[Yset circumferential_dist(i)];
%         end
%     end
%     save('Xset.mat','Xset');
%     save('Yset.mat','Yset');
%     [answer garbage]=fminsearch('line',[1 1]);
%     m=answer(1);
%     c=answer(2);
%     plot(Xset,Yset,'ro','LineWidth',2);
%     hold on;
%     Ynewset=m.*Xset+c;
%     plot(Xset,Ynewset,'b-','LineWidth',2);
% end
% saveas(imagewindow,'UnwrappedPointsFit.jpg');
% tic;pause(5);toc;
% end

   

