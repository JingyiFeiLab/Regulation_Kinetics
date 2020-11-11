% This is a script which does the analysis of the data generated from DBSCAN script run in the master_cluster script 
% Please send me a message at dgvjay@illinois.edu if you have any doubts about these scripts- Digvijay
                                                        
rmdir('ClusterFiles','s')                                               % Please make an empty folder called ClusterFiles which basically stores .mat file for each cluster but there is rmdir term here to remove any previous set from a previous data set and make a new one for new dataset
mkdir('ClusterFiles');
Radius_Gyration=[];                                                     % Stores the 3D(x,y,z considered)    Radius/Size of each cluster formed
Radius_Gyration2D=[];                                                   % Stores the 2D(only x,y considered) Radius/Size of each cluster formed
cluster_centers=[];                                                     % Stores the x,y,z center of each cluster formed
clusterpointsnumber=0;

for i=1:max(class)
     Clusterx=[];
     Clustery=[];
     Clusterz=[];
     CurrentClusterNumber=i;                                            % Keeps track of current cluster number as we go through all the points and find out about their cluster number generated from DBSCAN data
     Cluster=[];
     for j=1:numel(class)
       if class(j)==CurrentClusterNumber    
       Clusterx=[Clusterx data(j,1)];
       Clustery=[Clustery data(j,2)];
       Clusterz=[Clusterz data(j,3)];
       end
     end
     vector=CurrentClusterNumber*ones(numel(Clusterx),1);               % TotalClusters is a data file which has 4 columns..first 3 is x,y,z cooordinate and the 4th one is the column containing the cluster number of each data point
     Clusterx=Clusterx';
     Clustery=Clustery';
     Clusterz=Clusterz';
     clusterpointsnumber=clusterpointsnumber+numel(Clusterx);
     TotalClusters=[TotalClusters; [Clusterx Clustery Clusterz vector]];% TotalClusters file as described above
     Cluster      =[Cluster      ; [Clusterx Clustery Clusterz]];
     meanclusterx=mean(Clusterx);
     meanclustery=mean(Clustery);
     meanclusterz=mean(Clusterz);
     Rg=0;
     Rg2D=0;
     for k=1:numel(Clusterx)
         Rg=Rg+(Clusterx(k)-meanclusterx).^2+(Clustery(k)-meanclustery).^2+(Clusterz(k)-meanclusterz).^2;   % Calculating the Radius of each cluster(both 2D and 3D form as described above
         Rg2D=Rg2D+(Clusterx(k)-meanclusterx).^2+(Clustery(k)-meanclustery).^2;
     end
     if numel(meanclusterx)>=1
     Radius_Gyration=[Radius_Gyration;[sqrt(Rg/numel(Clusterx)) numel(Clusterx)]];                          % Stores the Radius of Gyration and the number of points in the particule cluster...this is later stored as a SHEETFILE
     Radius_Gyration2D=[Radius_Gyration2D;[sqrt(Rg2D/numel(Clusterx)) numel(Clusterx)]];                    % Stores the Radius of Gyration_2D and the number of points in the particule cluster...this is later stored as a SHEETFILE
     cluster_centers=[cluster_centers;[meanclusterx meanclustery meanclusterz sqrt(Rg/numel(Clusterx)) numel(Clusterx)]];
     clusterfile = sprintf('./ClusterFiles/Cluster_%d.mat',CurrentClusterNumber);
     save(clusterfile,'Cluster');
     end
end
   

if Radius_Gyration~=0      % i.e if there are any clusters formed....there could be a case where all the points are noise.
sheetfile2=sprintf('%s_CSVFiles_ColorCode%d/Sheet_%d_%d_%d.csv',dataset,ColorCode,Npts,Eps,ColorCode);
sheetfile2D=sprintf('%s_CSVFiles_ColorCode%d/Sheet_2D_%d_%d_%d.csv',dataset,ColorCode,Npts,Eps,ColorCode);
datfile=sprintf('%s_Coordinates_ColorCode%d/Coordinates_%d_%d_%d.dat',dataset,ColorCode,Npts,Eps,ColorCode);
datfile_centers=sprintf('%s_Cluster_Centers_ColorCode%d_%d_%d.dat',dataset,ColorCode,Npts,Eps);
% xlswrite(sheetfile2, Radius_Gyration) ;    %Stored the Radius of Gyration and the number of points in the particule cluster
% xlswrite(sheetfile2D, Radius_Gyration2D) ; %Stored the Radius of Gyration2D and the number of points in the particule cluster
dlmwrite(datfile_centers,cluster_centers,' ');
dlmwrite( datfile,TotalClusters,' ');
averagesheet=[averagesheet;[Npts Eps mean(Radius_Gyration(:,[1])) mean(Radius_Gyration2D(:,[1])) mean(Radius_Gyration(:,[2])) max(class) (3*clusterpointsnumber/numel(data))]];
%averegasheetfile stores : Npts Eps Average_Size_of_Cluster_at_the_given_parameters Average_Size_of_Cluster_(2D)_at_the_given_parameters Average_Number_of_Points_in_each_Cluster_at_the_given_parameters Number_of_Cluster_Formed_atthe_given_parameters %Fraction_of_points_which_are_in_clusters_atthe_givenparameters
end

