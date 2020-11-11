%% Minus sRNA a_m and kx fitting, 
close all

folderTitle = '/Users/reyer/Data/SingleCellEpi/MR156/plus6520/';
parentDir = '/Users/reyer/Data/SingleCellEpi/';
strains = [156];
time = {};
DatesV = {'June_5_2020_1','June_5_2020_2'};
graph_title = ' WT ptsG, + SgrS, 6120';
timestamps = {[1,2,3,4,5,6,7],[1,2,3,4,5,6,7]};

GFP = 1;
PTSG = 1; % 1 if you want to make ptsG comparisons
SGRS = 1; %1 if you want to make SgrS comparisons

sigma_check = 0; % If you want to include standard Deviation

mean_gfp = {};
mean_blue = {};
mean_ptsg = {};
mean_green = {};
mean_sgrs = {};
mean_red = {};
mean_violet = {};

std_gfp = {};
std_blue = {};
std_ptsg = {};
std_green = {};
std_sgrs = {};
std_red = {};
std_violet = {};

se_blue = {};
se_green = {};
se_red = {};


for i = 1:length(DatesV)
    strain = strcat(['MR',num2str(strains)]);
        

    strainDir = strcat([parentDir,strain]);
    if exist(strainDir,'dir')~=7
        mkdir(parentDir,strain{i})
    end
    
    if exist(folderTitle,'dir')~=7
        mkdir(folderTitle)
    end
   
    
    dataDir = strcat([parentDir,strain,'/',DatesV{i},'/']);
    compareDir = strcat([dataDir,'/Fits']);
    if exist(compareDir,'dir')~=7
        mkdir(dataDir,strcat(['Fits']))
    end
    dataFile = strcat([dataDir,'VThresh_BackgroundSummaryData.mat']);
    load(dataFile)
    time{i} = int_table{timestamps{1},1};
    
    if GFP == 1
        mean_gfp{i} = copy_table{timestamps{1},2};
        mean_blue{i} = int_table{timestamps{1},2};
        std_gfp{i} = copy_table{timestamps{1},3};
        std_blue{i} = int_table{timestamps{1},3};
        se_blue{i} = int_table{timestamps{1},4};
    end
    
    if SGRS == 1
        mean_sgrs{i} = copy_table{timestamps{1},6};
        mean_red{i} = int_table{timestamps{1},8};
        std_sgrs{i} = copy_table{timestamps{1},7};
        std_red{i} = int_table{timestamps{1},9};
        se_red{i} = int_table{timestamps{1},10};
    end
    
    if PTSG == 1
        mean_ptsg{i} = copy_table{timestamps{1},4};
        mean_green{i} = int_table{timestamps{1},5};
        std_ptsg{i} = copy_table{timestamps{1},5};
        std_green{i} = int_table{timestamps{1},6};
        se_green{i} = int_table{timestamps{1},7};
    end
    time
end

filename = strcat([compareDir,'/fit.mat']);

compare_green = zeros(length(DatesV),length(timestamps{1,1}));
compare_blue = zeros(length(DatesV),length(timestamps{1,1}));
compare_red = zeros(length(DatesV),length(timestamps{1,1}));

error_green = zeros(1,length(time{1,1}));
error_blue = zeros(1,length(time{1,1}));
error_red = zeros(1,length(time{1,1}));

replicate_green = zeros(1,length(time{1,1}));
replicate_blue = zeros(1,length(time{1,1}));
replicate_red = zeros(1,length(time{1,1}));

for i = 1:length(DatesV)
    %for j = 1:length(time{1,1})
    for j = 1:length(timestamps{1})
        if sum(j == [9] > 0) && i == 4
            compare_green(i,j) = mean_green{1,2}(j);
            compare_blue(i,j) = mean_blue{1,2}(j);
        else
            compare_green(i,j) = mean_green{1,i}(j);
            compare_blue(i,j) = mean_blue{1,i}(j);
        end
        
%         compare_green(i,j) = mean_green{1,i}(j);
%         compare_blue(i,j) = mean_blue{1,i}(j);
        
        if SGRS == 1
            compare_red(i,j) = mean_red{1,i}(j);
        else
            compare_red(i,j) = 0;
        end
    
    end
end


for i = 1:length(timestamps{1})
    replicate_green(i) = mean(compare_green(:,i));
    replicate_blue(i) = mean(compare_blue(:,i));
    replicate_red(i) = mean(compare_red(:,i));
    
    error_green(i) = std(compare_green(:,i))/1.414;
    error_blue(i) = std(compare_blue(:,i))/1.414;
    error_red(i) = std(compare_red(:,i))/2;
    
end



time = time{1}'*60; % update time axis, unit in second 
t_m = 360; % GFP maturation time, assuming 10min, unit in second
time_shift = time - t_m;

Legend=[];
for i = 1:length(DatesV)
    Legend= strcat([Legend DatesV{i}]);
end


figure(2)
errorbar(time,replicate_green,error_green,'k','LineWidth',3);
hold on
for i = 1:length(DatesV)
    
    
    
    txt = DatesV{i};
    plot(time,compare_green(i,:),'LineWidth',2,'DisplayName',txt)
    hold on
    
end
    
hold off
legend show
title(strcat([graph_title,' mRNA Reps']),'FontSize',24)
xlabel('Time after Induction (s)','FontSize',24)
set(gca,'YLim',[0 5.5E6])
set(gca,'XLim',[0,1800])
set(gca,'FontSize',18)
set(gcf,'position',[835,883,868,667]) 

file1 = strcat([folderTitle,'mRNA_reps']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'mRNA_reps.fig']);
savefig(gcf,file1_fig)
    
figure(4)
errorbar(time,replicate_blue,error_blue,'k','LineWidth',3);
hold on
for i = 1:length(DatesV)
    
    txt = DatesV{i};
    plot(time,compare_blue(i,:),'LineWidth',2,'DisplayName',txt)
    hold on
    
end
    
hold off
legend show
title(strcat([graph_title,' Protein Reps']),'FontSize',24)
xlabel('Time after Induction (s)','FontSize',24)
set(gca,'YLim',[0 1.5E7])
set(gca,'XLim',[0 1800])
set(gca,'FontSize',18)
set(gcf,'position',[835,883,868,667])

file1 = strcat([folderTitle,'protein_reps']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'protein_reps.fig']);
savefig(gcf,file1_fig)
    
figure(3)
errorbar(time,replicate_red,error_red,'k','LineWidth',3);
hold on
for i = 1:length(DatesV)
    
    txt = DatesV{i};
    plot(time,compare_red(i,:),'LineWidth',2,'DisplayName',txt)
    hold on
    
end
    
hold off
legend show
title(strcat([graph_title,' sRNA Reps']),'FontSize',24)
xlabel('Time after Induction (min)','FontSize',24)
set(gca,'YLim',[0 5.5E6])
set(gca,'XLim',[0 1800])
set(gca,'FontSize',18)
set(gcf,'position',[835,883,868,667])

file1 = strcat([folderTitle,'sRNA_reps']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([compareDir,'/',strain,'sRNA_fit.fig']);
savefig(gcf,file1_fig)

save(filename)


