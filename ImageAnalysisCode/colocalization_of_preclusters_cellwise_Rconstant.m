function summ=colocalization_of_preclusters_cellwise_Rconstant(...
    color2clusterfile,AlreadyClusteredColorCode,ColorCode2,...
    Rconstant,Scale_up_Cellular_Files_to_nm_by_this_Factor)

ClusterThatFoundOtherColocalizingCluster=0;
TotalNumberofPointsin_ClusterThatFoundOtherColocalizingCluster=0;
dataset=color2clusterfile;                                         
M=textread(dataset); 
ExpansivePoints=[M(:,1) M(:,2) M(:,3)];
OnlyClusterNumbers=M(:,4);
s = dir('*.dat'); 
file_list = {s.name}'; 
TotalSummary=[];
for i=1:numel(file_list)  
    
  List_ColorCode2_Cluster_Index=[];  
% color2clusterfile data file that has XYZ of all the clustered points of 
%colorcode2 along with their colorcode
% VERY IMPORTANT
% 
% All the .dat files in the folder will be taken as Cell Files so please
% make sure that the only .dat files in your folder are only Cell Files 
% 
% VERY IMPORTANT
% 
% The color2clusterfile has to be different from any other file this I 
%  create by changing it's extension to .dattt so please make this change in the name
% of the olor2clusterfile file before running this function
% AlreadyClusteredColorCode is the value of the "parameter colorcode" that has already been preclustered
% ColorCode2 is the value of the colorcode that needs to be colocalized 
% around the already pre clustered color code
% Rconstant is the extra bit of R allowed to accomodate more ColorCode 2 Points.
% Start_Npts,End_Npts,Increment_Npts,Start_Eps,End_Eps,Increment_Eps : 
% All these are the parameters for refinining the second colorcode which is being colocalized around the 
% already pre clustered colorcode

% Example of Usage:
% colocalization_of_preclusters_cellwise_Rconstant('1_Coordinates_2_15_1.dattt',1,3,40,130);
% ExpansivePoints=[];                                                   
% for t=1:numel(M)/4
% if M(t,4)==ColorCode2
% ExpansivePoints=[ExpansivePoints ;[M(t,1) M(t,2) M(t,3)]];                       
% end
% end
    dataset=file_list{i};
    CellNumber_Array= sscanf(dataset,'%c%c%c%c%c%d%c%c%c%c'); %why 10 components?
    CellNumber=CellNumber_Array(6);  %this is cell number from file name
    AlreadyClusteredPoints=textread(dataset);  %individual cell data are clustered already
    AlreadyClusteredPoints=[...
        AlreadyClusteredPoints(:,1)*Scale_up_Cellular_Files_to_nm_by_this_Factor...
        AlreadyClusteredPoints(:,2)*Scale_up_Cellular_Files_to_nm_by_this_Factor...
        AlreadyClusteredPoints(:,3)*Scale_up_Cellular_Files_to_nm_by_this_Factor...
        AlreadyClusteredPoints(:,4) AlreadyClusteredPoints(:,5) AlreadyClusteredPoints(:,6)];
        %the raw x,y,z data are in pixel values, so here changing back to nm scales
    TotalEntititiesin_AlreadyClusteredPoints=numel(AlreadyClusteredPoints)/6;
    % divided by 6, since there are six columns. This gives the total # of
    % clusters in ith cell

    TotalNumberofClustersfor_AlreadyClusteredPoints=0;
    Centers=[]; %centers of the first color (color1) clusters
    NumberofPoints=[];
    Radii=[];
    Entity_6thColumn=[];
    for kkk=1:TotalEntititiesin_AlreadyClusteredPoints
        if AlreadyClusteredPoints(kkk,6)==AlreadyClusteredColorCode;
            TotalNumberofClustersfor_AlreadyClusteredPoints=...
            TotalNumberofClustersfor_AlreadyClusteredPoints+1;  
            %counting the number of clusters in the first color, having the right parameter code.
            Radii=[Radii ;[AlreadyClusteredPoints(kkk,4)]];
            Centers=[Centers;[AlreadyClusteredPoints(kkk,1)...
            AlreadyClusteredPoints(kkk,2) AlreadyClusteredPoints(kkk,3)]];
            NumberofPoints=[NumberofPoints ; [AlreadyClusteredPoints(kkk,5)]] ;
            Entity_6thColumn=[Entity_6thColumn; [AlreadyClusteredPoints(kkk,6)]] ;
        end
    end
    SummaryList_Homo1=[];
    SummaryList_Homo2=[];
    SummaryList_Hetero=[];
    ReWrittenCell=[];
    summary_color_AlreadyPreClusteredPoints=[];
    summary_color_ExpansiveClusteredPoints=[];
    for NumberofClusters=1:TotalNumberofClustersfor_AlreadyClusteredPoints
        CurrentCenter=[Centers(NumberofClusters,1)...
            Centers(NumberofClusters,2) Centers(NumberofClusters,3)];

        % The most important line so what I am doing here is that finding the
        % distance of all the EXPANSIVEPOINTS from one of the cluster centers in
        % the cell_1.dat or cell_2.dat ( whatever file is being processed). The
        % EXPANSIVE POINTS are the points that are present in
        % '1_Coordinates_2_15_1.dattt'

        % After finding these distances, I find out all the points which lie within
        % the distance of sum of (radius of that particular cluster + Rconstant)

        DistanceColumn=pdist2(CurrentCenter,ExpansivePoints); 
            %here "CurrentCenter" is only a 3X1 vector (single cluster)
        DistanceColumn=DistanceColumn';
        Cell_Number_File=i;
        ClusterNumberofThatCell=NumberofClusters;
        Indexes1=find(DistanceColumn<(Rconstant));  %this runs up to the dimension of "ExpansivePoints"
        % I find the indexes/list of all the points that satisfy distance criterio
        Number_of_Points_of_ColorCode2_That_are_found_within_Colocalizing_Distance= numel(Indexes1);
        Fraction_of_Points_of_ColorCode2_That_are_found_within_Colocalizing_Distance=...
            numel(Indexes1)./((numel(ExpansivePoints)/3));
       
        if isempty(Indexes1)==0  % If no points were found none of the analysis ..no further analysis takes place
            ClusterThatFoundOtherColocalizingCluster=ClusterThatFoundOtherColocalizingCluster+1;
            summary_color_AlreadyPreClusteredPoints=...
                [summary_color_AlreadyPreClusteredPoints ;[NumberofClusters NumberofPoints(NumberofClusters)]];
            TotalNumberofPointsin_ClusterThatFoundOtherColocalizingCluster=...
                NumberofPoints(NumberofClusters)+TotalNumberofPointsin_ClusterThatFoundOtherColocalizingCluster;
            ColorCode2_Cluster_Index=OnlyClusterNumbers(Indexes1);  %% I find the cluster numbers of all the points indexed in line 81
            List_ColorCode2_Cluster_Index=[List_ColorCode2_Cluster_Index ;ColorCode2_Cluster_Index];  
            % I keep on adding the cluster numbers of color2 for each cluster center of color1 in an array
            Unique_ColorCode2_Cluster_Index=unique(ColorCode2_Cluster_Index);  
            % here I have made the list of all the unique cluster numbers ( from all the points that was found to colocalize
            %and then uniqued them to get the unique number of all the cluster numbers that are colocalizing .

            Indexes=[];
            ColoredPoints_Count=zeros(10,1);  %why ten components?
            for jj=1:numel(Unique_ColorCode2_Cluster_Index)
                Indexes=[Indexes (find(OnlyClusterNumbers==Unique_ColorCode2_Cluster_Index(jj)))'];  
                % As discussed...once the unique number of cluster numbers are found that partially or
                %fully colocalize...I collect all the points irrespective of whether they were originally 
                %colocalizing based on distance criterion or not since some other points in that
                %cluster family was found to colocalize based on the distance criterion

                ColoredPoints_Count(jj)=numel((find(OnlyClusterNumbers==Unique_ColorCode2_Cluster_Index(jj)))');
            end

            if isempty(Unique_ColorCode2_Cluster_Index)==0
                Tempo=[Centers(NumberofClusters,1) Centers(NumberofClusters,2)...
                    Centers(NumberofClusters,3) Radii(NumberofClusters)...
                    NumberofPoints(NumberofClusters) AlreadyClusteredPoints(NumberofClusters)];
                    %the last element is strange. I think it should be Entity_6thColumn(NumberofClusters)
                ReWrittenCell=[ ReWrittenCell; [Tempo numel(Unique_ColorCode2_Cluster_Index) sum(ColoredPoints_Count)]];
            end


            % NumberofPoints is an array which has the number of totalpoints that
            % belong to any cluster center mentioned in the cell_1.dat or cell_2.dat ( whatever file is being processed).
            % numel(Indexes) will give u the number of points that were found to colocalize at the given (Rconstant+ Radius of that particular cluster)
            % and based on these two numbers you can decide if it is homogeneous or
            % heteregenous cluster and that is what is being written down here.

            if NumberofPoints(NumberofClusters)/(NumberofPoints(NumberofClusters)+numel(Indexes))==0
                SummaryList_Homo1=[SummaryList_Homo1;[CurrentCenter Radii(NumberofClusters)...
                    NumberofPoints(NumberofClusters)/(NumberofPoints(NumberofClusters)+numel(Indexes))...
                    numel(Indexes)/(NumberofPoints(NumberofClusters)+numel(Indexes))...
                    numel(Unique_ColorCode2_Cluster_Index) ColoredPoints_Count(1) ColoredPoints_Count(2)...
                    ColoredPoints_Count(3) ColoredPoints_Count(4) ColoredPoints_Count(5) ]];
            elseif numel(Indexes)/(NumberofPoints(NumberofClusters)+numel(Indexes))==0
                SummaryList_Homo2=[SummaryList_Homo2;[CurrentCenter Radii(NumberofClusters)...
                    NumberofPoints(NumberofClusters)/(NumberofPoints(NumberofClusters)+numel(Indexes))...
                    numel(Indexes)/(NumberofPoints(NumberofClusters)+numel(Indexes))...
                    numel(Unique_ColorCode2_Cluster_Index) ColoredPoints_Count(1) ColoredPoints_Count(2)...
                    ColoredPoints_Count(3) ColoredPoints_Count(4) ColoredPoints_Count(5)]];       
            else
                SummaryList_Hetero=[SummaryList_Hetero;[CurrentCenter Radii(NumberofClusters)...
                    NumberofPoints(NumberofClusters)/(NumberofPoints(NumberofClusters)...
                    +numel(Indexes)) numel(Indexes)/(NumberofPoints(NumberofClusters)+numel(Indexes))...
                    numel(Unique_ColorCode2_Cluster_Index) ColoredPoints_Count(1) ColoredPoints_Count(2)...
                    ColoredPoints_Count(3) ColoredPoints_Count(4) ColoredPoints_Count(5)]];
            end
            TotalColorCode2Points=[];
            TotalColorCode2Points=[M(Indexes,1) M(Indexes,2) M(Indexes,3) M(Indexes,4)];
            colorcode2file= sprintf('%s_RConst%d_ALRClusCol%d__Col2%d_CluNum_%4.4d.datt',dataset,Rconstant,AlreadyClusteredColorCode,ColorCode2,NumberofClusters);
            % dlmwrite(colorcode2file,TotalColorCode2Points, ' ' );
        end
    end

    SummaryFileHomo1    = sprintf('%s_Summ_RConst%d_Homg_1_ALRClusCol%d_Col2%d.txt',dataset,Rconstant,AlreadyClusteredColorCode,ColorCode2);
    dlmwrite(SummaryFileHomo1,SummaryList_Homo1, ' ' );

    SummaryFileHomo2    = sprintf('%s_Summ_RConst%d_Homg_2_ALRClusCol%d_Col2%d.txt',dataset,Rconstant,AlreadyClusteredColorCode,ColorCode2);
    dlmwrite(SummaryFileHomo2,SummaryList_Homo2, ' ' );

    SummaryFileHetero    = sprintf('%s_Summ_RConst%d_Hetr_ALRClusCol%d_Col2%d.txt',dataset,Rconstant,AlreadyClusteredColorCode,ColorCode2);
    dlmwrite(SummaryFileHetero,SummaryList_Hetero, ' ' );


    summaryFile_color_AlreadyPreClusteredPoints    = sprintf('%s_Summ_RConst%d_summary_color_AlreadyPreClusteredPoints _ALRClusCol%d_Col2%d.txt',dataset,Rconstant,AlreadyClusteredColorCode,ColorCode2);
    dlmwrite(summaryFile_color_AlreadyPreClusteredPoints,summary_color_AlreadyPreClusteredPoints, ' ' );


    unique_List_ColorCode2_Cluster_Index=unique(List_ColorCode2_Cluster_Index); 
    % Pullup the array created in line 89 ..THIS HAS list of all the cluster numbers of color2...UNIQUE it to get rid of repetitive cluster numbers in the list


    % This loop is to go each of the cluster number in the list above and find out how many points belong to this cluster and put it in an array for to print later
    for jjj=1:numel(unique_List_ColorCode2_Cluster_Index)
    
        kkk=unique_List_ColorCode2_Cluster_Index(jjj); %unique cluster# for the 2nd image having overlap with 1st image clusters
        Indexes_kkk=find(OnlyClusterNumbers==kkk); %# finding spots having these clusters numbers
        summary_color_ExpansiveClusteredPoints=...
            [summary_color_ExpansiveClusteredPoints; [kkk numel(Indexes_kkk)]];
        %two column data: first column: unique cluster# from 2nd image
        %colocalizing with red clusters, 
        %2nd column: # of spots in each of these clusters
    end
    summaryFile_color_ExpansiveClusteredPoints=...
        sprintf('%s_Summ_RConst%d__color_ExpansiveClusteredPoints_ALRClusCol%d_Col2%d.txt',...
        dataset,Rconstant,AlreadyClusteredColorCode,ColorCode2);
    dlmwrite(summaryFile_color_ExpansiveClusteredPoints ,summary_color_ExpansiveClusteredPoints, ' ' );

    CellRewritten    = sprintf('Rewritten_%s_Rconstant_%d_ALRClusCol%d.txt',dataset,Rconstant,AlreadyClusteredColorCode);
    dlmwrite(CellRewritten,ReWrittenCell, ' ' );
    % zipfile_nonrefinedset=sprintf('CellNumber_%d_NonRef_ALRClusCol%d_Col2',CellNumber,AlreadyClusteredColorCode,ColorCode2);
    % zip(zipfile_nonrefinedset,'*.datt');
    % delete('*.datt');
    TotalSummary=[TotalSummary; [CellNumber TotalNumberofClustersfor_AlreadyClusteredPoints... % #of the first image (cell_i.dat) clusters having the right parameter code
        ClusterThatFoundOtherColocalizingCluster... % #of the first image clusters having colocalization with the 2nd clusters
        sum(NumberofPoints)... % #total # of points in the first image clusters
        TotalNumberofPointsin_ClusterThatFoundOtherColocalizingCluster... % #total # of points in the first image clusters, having colocalization
        ClusterThatFoundOtherColocalizingCluster/TotalNumberofClustersfor_AlreadyClusteredPoints... %the ratio between the first & second elements
        TotalNumberofPointsin_ClusterThatFoundOtherColocalizingCluster/sum(NumberofPoints)]]; %the ratio between the 3rd & 4th elements
    ClusterThatFoundOtherColocalizingCluster=0;
    TotalNumberofPointsin_ClusterThatFoundOtherColocalizingCluster=0;
end  %the end of checking each cell_i.dat
TotalSummaryname   = sprintf('TotalSummary_AlreadyClustered%d_Color2_%d_Rconstant',AlreadyClusteredColorCode,ColorCode2,Rconstant);
dlmwrite('Summmmary.txt',TotalSummary , ' ' );
summ=zeros(4,1);
summ(1,1)=mean(TotalSummary(:,6))*100; %percentage of colocalization by counting colocalizing clusters (color1 to color2 background)
summ(2,1)=mean(TotalSummary(:,7))*100; %very similar above, but it is the ratio of spots in colocalizing clusters
summ(3,1)=mean(TotalSummary(:,2)); %the total number of color1 clusters
summ(4,1)=mean(TotalSummary(:,3)); %the number of clusters in color1 that colocalizes with color2 background



end 