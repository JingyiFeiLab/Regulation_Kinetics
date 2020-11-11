parentDir = '/Users/reyer/Data/SingleCellEpi/';
strains = [187];
strain = {};
time = {};
Dates = {'May_16_2018'};
graph_title = 'manX Degradation ';



GFP = 1;
SGRS = 0; %1 if you want to make SgrS comparisons
PTSG = 1; % 1 if you want to make ptsG comparisons
VIOLET = 0;

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



for i = 1:length(strains)
    strain{i} = strcat(['MR',num2str(strains(i))]);
        

    strainDir = strcat([parentDir,strain_i]);
    if exist(strainDir,'dir')~=7
        mkdir(parentDir,strain_i)
    end
    
   
    
    dataDir = strcat([parentDir,strain{i},'/',Dates{i},'/']);
    backDir = strcat([dataDir,'/Background']);
    if exist(backDir,'dir')~=7
        mkdir(dataDir,'Background')
    end
    dataFile = strcat([dataDir,'BackgroundSummaryData.mat']);
    load(dataFile)
    time{i} = int_table{:,1};
    
    if GFP == 1
        mean_gfp{i} = copy_table{:,2};
        mean_blue{i} = int_table{:,2};
        std_gfp{i} = copy_table{:,3};
        std_blue{i} = int_table{:,3};
    end
    
    if SGRS == 1
        mean_sgrs{i} = copy_table{:,6};
        mean_red{i} = int_table{:,6};
        std_sgrs{i} = copy_table{:,7};
        std_red{i} = int_table{:,7};
    end
    
    if PTSG == 1
        mean_ptsg{i} = copy_table{:,4};
        mean_green{i} = int_table{:,4};
        std_ptsg{i} = copy_table{:,5};
        std_green{i} = int_table{:,5};
    end
    
    if VIOLET == 1
        
        mean_violet{i} = int_table{:,8};
        
        std_violet{i} = int_table{:,9};
    end
    
end

if PTSG == 1 && SGRS == 1 && sigma_check == 0
    close all
    figure(1)
   
    for i = 1:length(strains)
        
        strainDir = strcat([parentDir,strain_i]);
        if exist(strainDir,'dir')~=7
            mkdir(parentDir,strain_i)
        end
        
        dataDir = strcat([parentDir,strain{i},'/',Dates{i},'/']);
        backDir = strcat([dataDir,'/Background']);
        if exist(backDir,'dir')~=7
            mkdir(dataDir,'Background')
        end
        
        ylabels{1}='ptsG Intensity';
        ylabels{2}='SgrS Intensity';
        [ax,hlines(1),hlines(2)] = plotyy(time{i},mean_green{i},time{i},mean_red{i});
        cfig = get(gcf,'color');
        pos = [0.1  0.1  0.7  0.8];
        offset = pos(3)/5.5;
        hlines(1).LineWidth = 6;
        hlines(2).LineWidth = 6;
        set(get(ax(1),'ylabel'),'string',ylabels{1})
        set(get(ax(2),'ylabel'),'string',ylabels{2})
        
        
        
        title(strcat([graph_title,' ',strain{i}]),'FontSize',32)
        xlabel('Time after Induction (min)','FontSize',24)
        ax(1).FontSize = 18;
        ax(2).FontSize = 18;
        
        set(gcf,'position',[835,883,868,667])
        lgd2 = legend('ptsG','SgrS');
        lgd2.FontSize = 14;
        file1 = strcat([backDir,'/',strain{i},'ptsg_sgrs_intensity']);
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-painters','-depsc','-r0')
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-dpng','-r0')
        file1_fig = strcat([backDir,'/',strain{i},'ptsg_sgrs_intensity.fig']);
        savefig(gcf,file1_fig)
        close all
    
    end
end

if  PTSG == 1 && SGRS == 1 && sigma_check == 0
    close all
    figure(1)
   
    for i = 1:length(strains)
        
        strainDir = strcat([parentDir,strain_i]);
        if exist(strainDir,'dir')~=7
            mkdir(parentDir,strain_i)
        end
        
        dataDir = strcat([parentDir,strain{i},'/',Dates{i},'/']);
        backDir = strcat([dataDir,'/Background']);
        if exist(backDir,'dir')~=7
            mkdir(dataDir,'Background')
        end
        
        ylabels{1}='ptsG Intensity';
        ylabels{2}='SgrS Intensity';
        [ax,hlines(1),hlines(2)] = plotyy(time{i},mean_green{i},time{i},mean_red{i});
        cfig = get(gcf,'color');
        pos = [0.1  0.1  0.7  0.8];
        offset = pos(3)/5.5;
        hlines(1).LineWidth = 6;
        hlines(2).LineWidth = 6;
        set(get(ax(1),'ylabel'),'string',ylabels{1})
        set(get(ax(2),'ylabel'),'string',ylabels{2})
        
        
        
        title(strcat([graph_title,' ',strain{i}]),'FontSize',32)
        xlabel('Time after Induction (min)','FontSize',24)
        ax(1).FontSize = 18;
        ax(2).FontSize = 18;
        
        set(gcf,'position',[835,883,868,667])
        lgd2 = legend('ptsG','SgrS');
        lgd2.FontSize = 14;
        file1 = strcat([backDir,'/',strain{i},'ptsg_sgrs_intensity']);
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-painters','-depsc','-r0')
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-dpng','-r0')
        file1_fig = strcat([backDir,'/',strain{i},'ptsg_sgrs_intensity.fig']);
        savefig(gcf,file1_fig)
        close all
    
    end
end

if GFP == 1  && SGRS == 1 && sigma_check == 0
    close all
    figure(1)
   
    for i = 1:length(strains)
        
        strainDir = strcat([parentDir,strain_i]);
        if exist(strainDir,'dir')~=7
            mkdir(parentDir,strain_i)
        end
        
        dataDir = strcat([parentDir,strain{i},'/',Dates{i},'/']);
        backDir = strcat([dataDir,'/Background']);
        if exist(backDir,'dir')~=7
            mkdir(dataDir,'Background')
        end
        
        ylabels{1}='GFP Intensity';
        ylabels{2}='SgrS Intensity';
        [ax,hlines(1),hlines(2)] = plotyy(time{i},mean_blue{i},time{i},mean_red{i});
        cfig = get(gcf,'color');
        pos = [0.1  0.1  0.7  0.8];
        offset = pos(3)/5.5;
        hlines(1).LineWidth = 6;
        hlines(2).LineWidth = 6;
        set(get(ax(1),'ylabel'),'string',ylabels{1})
        set(get(ax(2),'ylabel'),'string',ylabels{2})
        
        
        
        title(strcat([graph_title,' ',strain{i}]),'FontSize',32)
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
        
        strainDir = strcat([parentDir,strain_i]);
        if exist(strainDir,'dir')~=7
            mkdir(parentDir,strain_i)
        end
        
        dataDir = strcat([parentDir,strain{i},'/',Dates{i},'/']);
        backDir = strcat([dataDir,'/Background']);
        if exist(backDir,'dir')~=7
            mkdir(dataDir,'Background')
        end
        
        ylabels{1}='ptsG Intensity';
        ylabels{2}='GFP Intensity';
        [ax,hlines(1),hlines(2)] = plotyy(time{i},mean_green{i},time{i},mean_blue{i});
        cfig = get(gcf,'color');
        pos = [0.1  0.1  0.7  0.8];
        offset = pos(3)/5.5;
        hlines(1).LineWidth = 6;
        hlines(2).LineWidth = 6;
        set(get(ax(1),'ylabel'),'string',ylabels{1})
        set(get(ax(2),'ylabel'),'string',ylabels{2})
        
        
        
        title(strcat([graph_title,' ',strain{i}]),'FontSize',32)
        xlabel('Time after Induction (min)','FontSize',24)
        ax(1).FontSize = 18;
        ax(2).FontSize = 18;
        
        set(gcf,'position',[835,883,868,667])
        lgd2 = legend('ptsG','GFP');
        lgd2.FontSize = 14;
        file1 = strcat([backDir,'/',strain{i},'ptsg_gfp_intensity']);
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-painters','-depsc','-r0')
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-dpng','-r0')
        file1_fig = strcat([backDir,'/',strain{i},'ptsg_gfp_intensity.fig']);
        savefig(gcf,file1_fig)
        close all
    
    end
end

if VIOLET == 1
    close all
    figure(1)
   
    for i = 1:length(strains)
        
        strainDir = strcat([parentDir,strain_i]);
        if exist(strainDir,'dir')~=7
            mkdir(parentDir,strain_i)
        end
        
        dataDir = strcat([parentDir,strain{i},'/',Dates{i},'/']);
        backDir = strcat([dataDir,'/Background']);
        if exist(backDir,'dir')~=7
            mkdir(dataDir,'Background')
        end
        
        ylabels{1}='rRNA Intensity';
        
        plot(time{i},mean_violet{i},'m','LineWidth',6);
        cfig = get(gcf,'color');
        pos = [0.1  0.1  0.7  0.8];
        offset = pos(3)/5.5;
        
        title(strcat([graph_title,' ',strain{i}]),'FontSize',32)
        xlabel('Time after Induction (min)','FontSize',24)
        
        
        
        set(gcf,'position',[835,883,868,667])
        lgd2 = legend('rRNA');
        lgd2.FontSize = 14;
        file1 = strcat([backDir,'/',strain{i},'rRNA_intensity']);
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-painters','-depsc','-r0')
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-dpng','-r0')
        file1_fig = strcat([backDir,'/',strain{i},'rRNA_intensity.fig']);
        savefig(gcf,file1_fig)
        close all
    
    end
end





