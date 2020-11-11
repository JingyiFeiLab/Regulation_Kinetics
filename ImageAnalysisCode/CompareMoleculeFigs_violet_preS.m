parentDir = '/Users/reyer/Data/SingleCellEpi/';
strains = [208];
strain = {};
time = {};
Dates = {'November_23_2018'};
graph_title = ' rne701 asd1, pre-Induced sRNA ';



GFP = 1;
SGRS = 1; %1 if you want to make SgrS comparisons
PTSG = 1; % 1 if you want to make ptsG comparisons

sigma_check = 1; % If you want to include standard Deviation

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



for i = 1:length(strains)
    strain{i} = strcat(['MR',num2str(strains(i))]);
        

    strainDir = strcat([parentDir,strain{i}]);
    if exist(strainDir,'dir')~=7
        mkdir(parentDir,strain{i})
    end
    
    dataDir = strcat([parentDir,strain{i},'/',Dates{i},'/']);
    backDir = strcat([dataDir,'/VioletBackground']);
    if exist(backDir,'dir')~=7
        mkdir(dataDir,'VioletBackground')
    end
    dataFile = strcat([dataDir,'VNorm_BackgroundSummaryData.mat']);
    load(dataFile)
    time{i} = int_table{:,1};
    
    if GFP == 1
        mean_gfp{i} = copy_table{:,2};
        mean_blue{i} = int_table{:,2};
        std_gfp{i} = copy_table{:,3};
        std_blue{i} = int_table{:,3};
        se_blue{i} = int_table{:,4};
    end
    
    if SGRS == 1
        mean_sgrs{i} = copy_table{:,6};
        mean_red{i} = int_table{:,8};
        std_sgrs{i} = copy_table{:,7};
        std_red{i} = int_table{:,9};
        se_red{i} = int_table{:,10};
    end
    
    if PTSG == 1
        mean_ptsg{i} = copy_table{:,4};
        mean_green{i} = int_table{:,5};
        std_ptsg{i} = copy_table{:,5};
        std_green{i} = int_table{:,6};
        se_green{i} = int_table{:,7};
    end
    
end

if PTSG == 1 && SGRS == 1 && sigma_check == 0
    close all
    figure(1)
   
    for i = 1:length(strains)
        
        strainDir = strcat([parentDir,strain{i}]);
        if exist(strainDir,'dir')~=7
            mkdir(parentDir,strain{i})
        end
        
        dataDir = strcat([parentDir,strain{i},'/',Dates{i},'/']);
        backDir = strcat([dataDir,'/VioletBackground']);
        if exist(backDir,'dir')~=7
            mkdir(dataDir,'VioletBackground')
        end
        
        ylabels{1}='mRNA Intensity';
        ylabels{2}='SgrS Intensity';
        [ax,hlines(1),hlines(2)] = plotyy(time{i},mean_green{i},time{i},mean_red{i});
        cfig = get(gcf,'color');
        pos = [0.1  0.1  0.7  0.8];
        offset = pos(3)/5.5;
        hlines(1).LineWidth = 6;
        hlines(2).LineWidth = 6;
        hlines(1).Color = 'g';
        hlines(2).Color = 'r';
        set(get(ax(1),'ylabel'),'string',ylabels{1})
        set(get(ax(2),'ylabel'),'string',ylabels{2})
        set(ax,{'ycolor'},{'g';'r'})
        set(ax(1),'YLim',[1.86E4 2.89E5])
        set(ax(2),'YLim',[4.33E5 2.83E6])
        
        
        
        title(strcat([graph_title]),'FontSize',32)
        xlabel('Time after Induction (min)','FontSize',24)
        ax(1).FontSize = 18;
        ax(2).FontSize = 18;
        
        set(gcf,'position',[835,883,868,667])
        lgd2 = legend('mRNA','SgrS');
        lgd2.FontSize = 14;
        file1 = strcat([backDir,'/',strain{i},'mRNA_sgrs_intensity']);
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-painters','-depsc','-r0')
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-dpng','-r0')
        file1_fig = strcat([backDir,'/',strain{i},'mRNA_sgrs_intensity.fig']);
        savefig(gcf,file1_fig)
        close all
    
    end
end

if GFP == 1  && SGRS == 1 && sigma_check == 0
    close all
    figure(1)
   
    for i = 1:length(strains)
        
        strainDir = strcat([parentDir,strain{i}]);
        if exist(strainDir,'dir')~=7
            mkdir(parentDir,strain_i)
        end
        
        dataDir = strcat([parentDir,strain{i},'/',Dates{i},'/']);
        backDir = strcat([dataDir,'/VioletBackground']);
        if exist(backDir,'dir')~=7
            mkdir(dataDir,'VioletBackground')
        end
        
        ylabels{1}='GFP Intensity';
        ylabels{2}='SgrS Intensity';
        [ax,hlines(1),hlines(2)] = plotyy(time{i},mean_blue{i},time{i},mean_red{i});
        cfig = get(gcf,'color');
        pos = [0.1  0.1  0.7  0.8];
        offset = pos(3)/5.5;
        hlines(1).LineWidth = 6;
        hlines(2).LineWidth = 6;
        hlines(1).Color = 'b';
        hlines(2).Color = 'r';
        set(get(ax(1),'ylabel'),'string',ylabels{1})
        set(get(ax(2),'ylabel'),'string',ylabels{2})
        set(ax,{'ycolor'},{'b';'r'})
        set(ax(1),'YLim',[4.79E3 1.10E6])
        set(ax(2),'YLim',[4.33E5 2.83E6])
        
        
        title(strcat([graph_title]),'FontSize',32)
        xlabel('Time after Induction (min)','FontSize',24)
        ax(1).FontSize = 18;
        ax(2).FontSize = 18;
        
        set(gcf,'position',[835,883,868,667])
        lgd2 = legend('GFP','SgrS');
        lgd2.FontSize = 14;
        file1 = strcat([backDir,'/',strain{i},'gfp_sgrs_intensity']);
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-painters','-depsc','-r0')
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-dpng','-r0')
        file1_fig = strcat([backDir,'/',strain{i},'gfp_sgrs_intensity.fig']);
        savefig(gcf,file1_fig)
        close all
    
    end
end

if GFP == 1 && PTSG == 1 && sigma_check == 0
    close all
    figure(1)
   
    for i = 1:length(strains)
        
        strainDir = strcat([parentDir,strain{i}]);
        if exist(strainDir,'dir')~=7
            mkdir(parentDir,strain{i})
        end
        
        dataDir = strcat([parentDir,strain{i},'/',Dates{i},'/']);
        backDir = strcat([dataDir,'/VioletBackground']);
        if exist(backDir,'dir')~=7
            mkdir(dataDir,'VioletBackground')
        end
        
        ylabels{1}='mRNA Intensity';
        ylabels{2}='GFP Intensity';
        [ax,hlines(1),hlines(2)] = plotyy(time{i},mean_green{i},time{i},mean_blue{i});
        cfig = get(gcf,'color');
        pos = [0.1  0.1  0.7  0.8];
        offset = pos(3)/5.5;
        hlines(1).LineWidth = 6;
        hlines(2).LineWidth = 6;
        hlines(1).Color = 'g';
        hlines(2).Color = 'b';
        set(get(ax(1),'ylabel'),'string',ylabels{1})
        set(get(ax(2),'ylabel'),'string',ylabels{2})
        set(ax,{'ycolor'},{'g';'b'})
        set(ax(1),'YLim',[1.86E4 2.89E5])
        set(ax(2),'YLim',[4.79E3 1.10E6])
        
        
        title(strcat([graph_title]),'FontSize',32)
        xlabel('Time after Induction (min)','FontSize',24)
        ax(1).FontSize = 18;
        ax(2).FontSize = 18;
        
        set(gcf,'position',[835,883,868,667])
        lgd2 = legend('mRNA','GFP');
        lgd2.FontSize = 14;
        file1 = strcat([backDir,'/',strain{i},'mRNA_gfp_intensity']);
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-painters','-depsc','-r0')
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-dpng','-r0')
        file1_fig = strcat([backDir,'/',strain{i},'mRNA_gfp_intensity.fig']);
        savefig(gcf,file1_fig)
        close all
    
    end
end

if PTSG == 1 && SGRS == 1 && sigma_check == 1
    close all
    figure(1)
   
    for i = 1:length(strains)
        
        strainDir = strcat([parentDir,strain{i}]);
        if exist(strainDir,'dir')~=7
            mkdir(parentDir,strain{i})
        end
        
        dataDir = strcat([parentDir,strain{i},'/',Dates{i},'/']);
        backDir = strcat([dataDir,'/VioletBackground']);
        if exist(backDir,'dir')~=7
            mkdir(dataDir,'VioletBackground')
        end
        
        ylabels{1}='mRNA Intensity';
        ylabels{2}='SgrS Intensity';
        [ax,hlines(1),hlines(2)] = plotyy(time{i},mean_green{i},time{i},mean_red{i});
        cfig = get(gcf,'color');
        pos = [0.1  0.1  0.7  0.8];
        offset = pos(3)/5.5;
        hlines(1).LineWidth = 6;
        hlines(2).LineWidth = 6;
        hlines(1).Color = 'g';
        hlines(2).Color = 'r';
        set(get(ax(1),'ylabel'),'string',ylabels{1})
        set(get(ax(2),'ylabel'),'string',ylabels{2})
        set(ax,{'ycolor'},{'g';'r'})
        set(ax(1),'YLim',[0E5 4E6])
        set(ax(2),'YLim',[0 3E6])
        hold(ax(1),'on')
        scatter(ax(1),time{i},mean_green{i},200,'g','filled')
        errorbar(ax(1),time{i},mean_green{i},se_green{i},'g','LineStyle','none');
        hold(ax(2),'on')
        scatter(ax(2),time{i},mean_red{i},200,'r','filled')
        errorbar(ax(2),time{i},mean_red{i},se_red{i},'r','LineStyle','none');
        
        
        
        
        title(strcat([graph_title]),'FontSize',32)
        xlabel('Time after Induction (min)','FontSize',24)
        ax(1).FontSize = 18;
        ax(2).FontSize = 18;
        
        set(gcf,'position',[835,883,868,667])
        lgd2 = legend('mRNA','SgrS');
        lgd2.FontSize = 14;
        file1 = strcat([backDir,'/',strain{i},'mRNA_sgrs_intensity']);
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-painters','-depsc','-r0')
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-dpng','-r0')
        file1_fig = strcat([backDir,'/',strain{i},'mRNA_sgrs_intensity.fig']);
        savefig(gcf,file1_fig)
        close all
    
    end
end

if GFP == 1  && SGRS == 1 && sigma_check == 1
    close all
    figure(1)
   
    for i = 1:length(strains)
        
        strainDir = strcat([parentDir,strain{i}]);
        if exist(strainDir,'dir')~=7
            mkdir(parentDir,strain_i)
        end
        
        dataDir = strcat([parentDir,strain{i},'/',Dates{i},'/']);
        backDir = strcat([dataDir,'/VioletBackground']);
        if exist(backDir,'dir')~=7
            mkdir(dataDir,'VioletBackground')
        end
        
        ylabels{1}='GFP Intensity';
        ylabels{2}='SgrS Intensity';
        [ax,hlines(1),hlines(2)] = plotyy(time{i},mean_blue{i},time{i},mean_red{i});
        cfig = get(gcf,'color');
        pos = [0.1  0.1  0.7  0.8];
        offset = pos(3)/5.5;
        hlines(1).LineWidth = 6;
        hlines(2).LineWidth = 6;
        hlines(1).Color = 'b';
        hlines(2).Color = 'r';
        set(get(ax(1),'ylabel'),'string',ylabels{1})
        set(get(ax(2),'ylabel'),'string',ylabels{2})
        set(ax,{'ycolor'},{'b';'r'})
        set(ax(1),'YLim',[1E5 3E7])
        set(ax(2),'YLim',[0 3E6])
        hold(ax(1),'on')
        scatter(ax(1),time{i},mean_blue{i},200,'b','filled')
        errorbar(ax(1),time{i},mean_blue{i},se_blue{i},'b','LineStyle','none');
        hold(ax(2),'on')
        scatter(ax(2),time{i},mean_red{i},200,'r','filled')
        errorbar(ax(2),time{i},mean_red{i},se_red{i},'r','LineStyle','none');
        
        
        title(strcat([graph_title]),'FontSize',32)
        xlabel('Time after Induction (min)','FontSize',24)
        ax(1).FontSize = 18;
        ax(2).FontSize = 18;
        
        set(gcf,'position',[835,883,868,667])
        lgd2 = legend('GFP','SgrS');
        lgd2.FontSize = 14;
        file1 = strcat([backDir,'/',strain{i},'gfp_sgrs_intensity']);
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-painters','-depsc','-r0')
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-dpng','-r0')
        file1_fig = strcat([backDir,'/',strain{i},'gfp_sgrs_intensity.fig']);
        savefig(gcf,file1_fig)
        close all
    
    end
end

if GFP == 1 && PTSG == 1 && sigma_check == 1
    close all
    figure(1)
   
    for i = 1:length(strains)
        
        strainDir = strcat([parentDir,strain{i}]);
        if exist(strainDir,'dir')~=7
            mkdir(parentDir,strain{i})
        end
        
        dataDir = strcat([parentDir,strain{i},'/',Dates{i},'/']);
        backDir = strcat([dataDir,'/VioletBackground']);
        if exist(backDir,'dir')~=7
            mkdir(dataDir,'VioletBackground')
        end
        
        ylabels{1}='mRNA Intensity';
        ylabels{2}='GFP Intensity';
        [ax,hlines(1),hlines(2)] = plotyy(time{i},mean_green{i},time{i},mean_blue{i});
        cfig = get(gcf,'color');
        pos = [0.1  0.1  0.7  0.8];
        offset = pos(3)/5.5;
        hlines(1).LineWidth = 6;
        hlines(2).LineWidth = 6;
        hlines(1).Color = 'g';
        hlines(2).Color = 'b';
        set(get(ax(1),'ylabel'),'string',ylabels{1})
        set(get(ax(2),'ylabel'),'string',ylabels{2})
        set(ax,{'ycolor'},{'g';'b'})
        set(ax(1),'YLim',[0E5 4E6])
        set(ax(2),'YLim',[1E5 3E7])
        hold(ax(1),'on')
        scatter(ax(1),time{i},mean_green{i},200,'g','filled')
        errorbar(ax(1),time{i},mean_green{i},se_green{i},'g','LineStyle','none');
        hold(ax(2),'on')
        scatter(ax(2),time{i},mean_blue{i},200,'b','filled')
        errorbar(ax(2),time{i},mean_blue{i},se_blue{i},'b','LineStyle','none');
        
        
        title(strcat([graph_title]),'FontSize',32)
        xlabel('Time after Induction (min)','FontSize',24)
        ax(1).FontSize = 18;
        ax(2).FontSize = 18;
        
        set(gcf,'position',[835,883,868,667])
        lgd2 = legend('mRNA','GFP');
        lgd2.FontSize = 14;
        file1 = strcat([backDir,'/',strain{i},'mRNA_gfp_intensity']);
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-painters','-depsc','-r0')
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-dpng','-r0')
        file1_fig = strcat([backDir,'/',strain{i},'mRNA_gfp_intensity.fig']);
        savefig(gcf,file1_fig)
        close all
    
    end
end





