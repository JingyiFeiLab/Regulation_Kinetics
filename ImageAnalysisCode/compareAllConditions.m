parentDir = '/Users/reyer/Data/SingleCellEpi/';
strain_wt = [187];
strain_rne = [184];
time = [0,1,3,6,12,18,24,30];
Dates_wt_plus = {'September_19_2018'};
Dates_wt_minus = {'September_13_2018'};
Dates_rne_plus = {'September_20_2018'};
Dates_rne_minus = {'September_15_2018'};
graph_title = ' manX Replicates ';

file = '/Users/reyer/Data/SingleCellEpi/manXcomparisons/';


GFP = 1;
SGRS = 1; %1 if you want to make SgrS comparisons
PTSG = 1; % 1 if you want to make ptsG comparisons

sigma_check = 1; % If you want to include standard Deviation


mean_blue_wt_plus = {};
mean_green_wt_plus = {};
mean_red_wt_plus = {};

mean_blue_wt_minus = {};
mean_green_wt_minus = {};
mean_red_wt_minus = {};

mean_blue_rne_plus = {};
mean_green_rne_plus = {};
mean_red_rne_plus = {};

mean_blue_rne_minus = {};
mean_green_rne_minus = {};
mean_red_rne_minus = {};



for i = 1:length(Dates_wt_minus)
    strain = strcat(['MR',num2str(strain_wt)]);
        
    dataDir = strcat([parentDir,strain,'/',Dates_wt_minus{i},'/']);
    
    dataFile = strcat([dataDir,'VNorm_BackgroundSummaryData.mat']);
    load(dataFile)
    
    for j = 1:length(time)
        
        if sum(int_table{:,1}==time(j))
            
            if GFP == 1
                mean_blue_wt_minus{i,j} = int_table{int_table{:,1}==time(j),2};
            end
            
            if PTSG == 1
                mean_green_wt_minus{i,j} = int_table{int_table{:,1}==time(j),5};
            end
            
            if SGRS == 1
                mean_red_wt_minus{i,j} = int_table{int_table{:,1}==time(j),8};
            end
        end
    end

end

for i = 1:length(Dates_wt_plus)
    strain = strcat(['MR',num2str(strain_wt)]);
        
    dataDir = strcat([parentDir,strain,'/',Dates_wt_plus{i},'/']);
    
    dataFile = strcat([dataDir,'VNorm_BackgroundSummaryData.mat']);
    load(dataFile)
    
    for j = 1:length(time)
        
        if sum(int_table{:,1}==time(j))
            
            if GFP == 1
                mean_blue_wt_plus{i,j} = int_table{int_table{:,1}==time(j),2};
            end
            
            if PTSG == 1
                mean_green_wt_plus{i,j} = int_table{int_table{:,1}==time(j),5};
            end
            
            if SGRS == 1
                mean_red_wt_plus{i,j} = int_table{int_table{:,1}==time(j),8};
            end
        end
    end

end

for i = 1:length(Dates_rne_plus)
    strain = strcat(['MR',num2str(strain_rne)]);
        
    dataDir = strcat([parentDir,strain,'/',Dates_rne_plus{i},'/']);
    
    dataFile = strcat([dataDir,'VNorm_BackgroundSummaryData.mat']);
    load(dataFile)
    
    for j = 1:length(time)
        
        if sum(int_table{:,1}==time(j))
            
            if GFP == 1
                mean_blue_rne_plus{i,j} = int_table{int_table{:,1}==time(j),2};
            end
            
            if PTSG == 1
                mean_green_rne_plus{i,j} = int_table{int_table{:,1}==time(j),5};
            end
            
            if SGRS == 1
                mean_red_rne_plus{i,j} = int_table{int_table{:,1}==time(j),8};
            end
        end
    end

end

for i = 1:length(Dates_rne_minus)
    strain = strcat(['MR',num2str(strain_rne)]);
        
    dataDir = strcat([parentDir,strain,'/',Dates_rne_minus{i},'/']);
    
    dataFile = strcat([dataDir,'VNorm_BackgroundSummaryData.mat']);
    load(dataFile)
    
    for j = 1:length(time)
        
        if sum(int_table{:,1}==time(j))
            
            if GFP == 1
                mean_blue_rne_minus{i,j} = int_table{int_table{:,1}==time(j),2};
            end
            
            if PTSG == 1
                mean_green_rne_minus{i,j} = int_table{int_table{:,1}==time(j),5};
            end
            
            if SGRS == 1
                mean_red_rne_minus{i,j} = int_table{int_table{:,1}==time(j),8};
            end
        end
    end

end

error_green_wt_plus = zeros(1,length(time));
error_blue_wt_plus = zeros(1,length(time));
error_red_wt_plus = zeros(1,length(time));

replicate_green_wt_plus = zeros(1,length(time));
replicate_blue_wt_plus = zeros(1,length(time));
replicate_red_wt_plus = zeros(1,length(time));

for i = 1:length(time)
    replicate_green_wt_plus(i) = mean([mean_green_wt_plus{:,i}]);
    replicate_blue_wt_plus(i) = mean([mean_blue_wt_plus{:,i}]);
    replicate_red_wt_plus(i) = mean([mean_red_wt_plus{:,i}]);
    
    error_green_wt_plus(i) = std([mean_green_wt_plus{:,i}])/2;
    error_blue_wt_plus(i) = std([mean_blue_wt_plus{:,i}])/2;
    error_red_wt_plus(i) = std([mean_red_wt_plus{:,i}])/2;
    
end

error_green_wt_minus = zeros(1,length(time));
error_blue_wt_minus = zeros(1,length(time));
error_red_wt_minus = zeros(1,length(time));

replicate_green_wt_minus = zeros(1,length(time));
replicate_blue_wt_minus = zeros(1,length(time));
replicate_red_wt_minus = zeros(1,length(time));

for i = 1:length(time)
    replicate_green_wt_minus(i) = mean([mean_green_wt_minus{:,i}]);
    replicate_blue_wt_minus(i) = mean([mean_blue_wt_minus{:,i}]);
    replicate_red_wt_minus(i) = mean([mean_red_wt_minus{:,i}]);
    
    error_green_wt_minus(i) = std([mean_green_wt_minus{:,i}])/2;
    error_blue_wt_minus(i) = std([mean_blue_wt_minus{:,i}])/2;
    error_red_wt_minus(i) = std([mean_red_wt_minus{:,i}])/2;
    
end

error_green_rne_plus = zeros(1,length(time));
error_blue_rne_plus = zeros(1,length(time));
error_red_rne_plus = zeros(1,length(time));

replicate_green_rne_plus = zeros(1,length(time));
replicate_blue_rne_plus = zeros(1,length(time));
replicate_red_rne_plus = zeros(1,length(time));

for i = 1:length(time)
    replicate_green_rne_plus(i) = mean([mean_green_rne_plus{:,i}]);
    replicate_blue_rne_plus(i) = mean([mean_blue_rne_plus{:,i}]);
    replicate_red_rne_plus(i) = mean([mean_red_rne_plus{:,i}]);
    
    error_green_rne_plus(i) = std([mean_green_rne_plus{:,i}])/2;
    error_blue_rne_plus(i) = std([mean_blue_rne_plus{:,i}])/2;
    error_red_rne_plus(i) = std([mean_red_rne_plus{:,i}])/2;
    
end

error_green_rne_minus = zeros(1,length(time));
error_blue_rne_minus = zeros(1,length(time));
error_red_rne_minus = zeros(1,length(time));

replicate_green_rne_minus = zeros(1,length(time));
replicate_blue_rne_minus = zeros(1,length(time));
replicate_red_rne_minus = zeros(1,length(time));

for i = 1:length(time)
    replicate_green_rne_minus(i) = mean([mean_green_rne_minus{:,i}]);
    replicate_blue_rne_minus(i) = mean([mean_blue_rne_minus{:,i}]);
    replicate_red_rne_minus(i) = mean([mean_red_rne_minus{:,i}]);
    
    error_green_rne_minus(i) = std([mean_green_rne_minus{:,i}])/2;
    error_blue_rne_minus(i) = std([mean_blue_rne_minus{:,i}])/2;
    error_red_rne_minus(i) = std([mean_red_rne_minus{:,i}])/2;
    
end


if PTSG == 1 && sigma_check == 1
    close all
    figure(1)
    
    errorbar(time,replicate_green_wt_plus,error_green_wt_plus,'r', 'LineWidth',3);
    hold on
    errorbar(time,replicate_green_wt_minus,error_green_wt_minus,'g', 'LineWidth',3);
    hold on
    errorbar(time,replicate_green_rne_plus,error_green_rne_plus,'b', 'LineWidth',3);
    hold on
    errorbar(time,replicate_green_rne_minus,error_green_rne_minus,'k', 'LineWidth',3);
    hold on
    
    set(gca,'XLim',[0 30])
    xlabel('Time','FontSize',24)
    ylabel('GFP','FontSize',24)
    title(strcat([graph_title,'mRNA']),'FontSize',28)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 18)
    lgd = legend('WT + SgrS','WT -SgrS','rne701 +SgrS','rne701 -SgrS');
    lgd.Location = 'northeast';
    set(gcf,'position',[835,883,868,667])
    file1 = strcat([file,'mRNA']);
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    print(file1,'-painters','-depsc','-r0')
    file1_fig = strcat([file,'mRNA.fig']);
    savefig(gcf,file1_fig)
    
    
end

if SGRS == 1 && sigma_check == 1
    
    close all
    figure(2)
    
    errorbar(time,replicate_red_wt_plus,error_red_wt_plus,'r', 'LineWidth',3);
    hold on
    errorbar(time,replicate_red_wt_minus,error_red_wt_minus,'g', 'LineWidth',3);
    hold on
    errorbar(time,replicate_red_rne_plus,error_red_rne_plus,'b', 'LineWidth',3);
    hold on
    errorbar(time,replicate_red_rne_minus,error_red_rne_minus,'k', 'LineWidth',3);
    hold on
    
    set(gca,'XLim',[0 30])
    xlabel('Time','FontSize',24)
    ylabel('GFP','FontSize',24)
    title(strcat([graph_title,'SgrS']),'FontSize',28)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 18)
    lgd = legend('WT + SgrS','WT -SgrS','rne701 +SgrS','rne701 -SgrS');
    lgd.Location = 'northeast';
    set(gcf,'position',[835,883,868,667])
    file1 = strcat([file,'SgrS']);
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    print(file1,'-painters','-depsc','-r0')
    file1_fig = strcat([file,'SgrS.fig']);
    savefig(gcf,file1_fig)
    
    close all
        
end

if GFP == 1 && sigma_check == 1
    close all
    
    figure(3)
    
    errorbar(time,replicate_blue_wt_plus,error_blue_wt_plus,'r', 'LineWidth',3);
    hold on
    errorbar(time,replicate_blue_wt_minus,error_blue_wt_minus,'g', 'LineWidth',3);
    hold on
    errorbar(time,replicate_blue_rne_plus,error_blue_rne_plus,'b', 'LineWidth',3);
    hold on
    errorbar(time,replicate_blue_rne_minus,error_blue_rne_minus,'k', 'LineWidth',3);
    hold on
    
    set(gca,'XLim',[0 30])
    xlabel('Time','FontSize',24)
    ylabel('GFP','FontSize',24)
    title(strcat([graph_title,'protein']),'FontSize',28)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 18)
    lgd = legend('WT + SgrS','WT -SgrS','rne701 +SgrS','rne701 -SgrS');
    lgd.Location = 'northeast';
    set(gcf,'position',[835,883,868,667])
    file1 = strcat([file,'protein']);
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    print(file1,'-painters','-depsc','-r0')
    file1_fig = strcat([file,'protein.fig']);
    savefig(gcf,file1_fig)
    close all
end



