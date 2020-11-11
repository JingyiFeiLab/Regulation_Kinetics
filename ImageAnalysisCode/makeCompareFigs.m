parentDir = '/Users/reyer/Data/SingleCellEpi/';
strains = [60,94,97,105,148];
strain = {};
time = {};
Dates = {'March_11_2018','March_11_2018','March_11_2018','March_11_2018','March_11_2018',};
newFile = 'RBS_mutants_combined_5ng';
graph_title = 'RBS Mutants, Degradation';

compareDir = strcat([parentDir,newFile]);
if exist(compareDir,'dir')~=7
    mkdir(parentDir, newFile)
end

GFP = 1; % 1 if you want to make GFP comparisons
SGRS = 0; %1 if you want to make SgrS comparisons
PTSG = 1; % 1 if you want to make ptsG comparisons

sigma_check = 0; % If you want to include standard Deviation

mean_gfp = {};
mean_blue = {};
mean_ptsg = {};
mean_green = {};
mean_sgrs = {};
mean_red = {};

std_gfp = {};
std_blue = {};
std_ptsg = {};
std_green = {};
std_sgrs = {};
std_red = {};


for i = 1:length(strains)
    strain{i} = strcat(['MR',num2str(strains(i))]);
    dataDir = strcat([parentDir,strain{i},'/',Dates{i},'/']);
    dataFile = strcat([dataDir,'summaryData.mat']);
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
    
end

if GFP == 1 && sigma_check == 0
    close all
    figure(1)
    for i = 1:length(strains)
        plot(time{i},mean_gfp{i},'LineWidth',6)
        hold on
    end
    
    
    title(strcat([graph_title,' GFP']),'FontSize',32)
    xlabel('Time after Induction (min)','FontSize',24)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 18)
    ylabel('GFP Copy Number','FontSize',24)
    set(gcf,'position',[835,883,868,667])
    lgd = legend(strain);
    lgd.FontSize = 14;
    file1 = strcat([compareDir,'/','GFP']);
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([compareDir,'/GFP.fig']);
    savefig(gcf,file1_fig)
    
    figure(2)
    for i = 1:length(strains)
        plot(time{i},mean_blue{i},'LineWidth',6)
        hold on
    end
    
    
    title(strcat([graph_title,' Blue']),'FontSize',32)
    xlabel('Time after Induction (min)','FontSize',24)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 18)
    ylabel('GFP Whole Cell Intensity','FontSize',24)
    set(gcf,'position',[835,883,868,667])
    lgd2 = legend(strain);
    lgd2.FontSize = 14;
    file1 = strcat([compareDir,'/BLUE']);
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([compareDir,'BLUE.fig']);
    savefig(gcf,file1_fig)

elseif GFP == 1 && sigma_check == 1
    figure(1)
    for i = 1:length(strains)
        errorbar(time{i},mean_gfp{i},std_gfp{i},'LineWidth',6)
        hold on
    end
    
    
    title(strcat([graph_title,' GFP']),'FontSize',32)
    xlabel('Time after Induction (min)','FontSize',24)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 18)
    ylabel('GFP Copy Number','FontSize',24)
    set(gcf,'position',[835,883,868,667])
    lgd = legend(strain);
    lgd.FontSize = 14;
    file1 = strcat([compareDir,'/','GFP']);
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([compareDir,'/GFP.fig']);
    savefig(gcf,file1_fig)
    
    figure(2)
    for i = 1:length(strains)
        errorbar(time{i},mean_blue{i},std_gfp{i},'LineWidth',6)
        hold on
    end
    
    
    title(strcat([graph_title,' Blue']),'FontSize',32)
    xlabel('Time after Induction (min)','FontSize',24)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 18)
    ylabel('GFP Whole Cell Intensity','FontSize',24)
    set(gcf,'position',[835,883,868,667])
    lgd = legend(strain);
    lgd.FontSize = 14;
    file1 = strcat([compareDir,'/BLUE']);
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([compareDir,'BLUE.fig']);
    savefig(gcf,file1_fig)
    
end

if SGRS == 1 && sigma_check == 0
    figure(1)
    for i = 1:length(strains)
        plot(time{i},mean_sgrs{i},'LineWidth',6)
        hold on
    end
    
    
    title(strcat([graph_title,' SgrS']),'FontSize',32)
    xlabel('Time after Induction (min)','FontSize',24)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 18)
    ylabel('SgrS Copy Number','FontSize',24)
    set(gcf,'position',[835,883,868,667])
    lgd = legend(strain);
    lgd.FontSize = 14;
    file1 = strcat([compareDir,'/','SGRS']);
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([compareDir,'/SGRS.fig']);
    savefig(gcf,file1_fig)
    
    figure(2)
    for i = 1:length(strains)
        plot(time{i},mean_red{i},'LineWidth',6)
        hold on
    end
    
    
    title(strcat([graph_title,' Red']),'FontSize',32)
    xlabel('Time after Induction (min)','FontSize',24)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 18)
    ylabel('SgrS Whole Cell Intensity','FontSize',24)
    set(gcf,'position',[835,883,868,667])
    lgd = legend(strain);
    lgd.FontSize = 14;
    file1 = strcat([compareDir,'/RED']);
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([compareDir,'RED.fig']);
    savefig(gcf,file1_fig)

elseif SGRS == 1 && sigma_check == 1
    figure(1)
    for i = 1:length(strains)
        errorbar(time{i},mean_sgrs{i},std_sgrs{i},'LineWidth',6)
        hold on
    end
    
    
    title(strcat([graph_title,' SgrS']),'FontSize',32)
    xlabel('Time after Induction (min)','FontSize',24)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 18)
    ylabel('SgrS Copy Number','FontSize',24)
    set(gcf,'position',[835,883,868,667])
    lgd = legend(strain);
    lgd.FontSize = 14;
    file1 = strcat([compareDir,'/','SgrS']);
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([compareDir,'/SgrS.fig']);
    savefig(gcf,file1_fig)
    
    figure(2)
    for i = 1:length(strains)
        plot(time{i},mean_red{i},std_sgrs{i},'LineWidth',6)
        hold on
    end
    
    
    title(strcat([graph_title,' Red']),'FontSize',32)
    xlabel('Time after Induction (min)','FontSize',24)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 18)
    ylabel('SgrS Whole Cell Intensity','FontSize',24)
    set(gcf,'position',[835,883,868,667])
    lgd = legend(strain);
    lgd.FontSize = 14;
    file1 = strcat([compareDir,'/Red']);
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([compareDir,'Red.fig']);
    savefig(gcf,file1_fig)
    
end

if PTSG == 1 && sigma_check == 0
    close all
    figure(1)
    for i = 1:length(strains)
        plot(time{i},mean_ptsg{i},'LineWidth',6)
        hold on
    end
    
    
    title(strcat([graph_title,' ptsG']),'FontSize',32)
    xlabel('Time after Induction (min)','FontSize',24)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 18)
    ylabel('ptsG Copy Number','FontSize',24)
    set(gcf,'position',[835,883,868,667])
    lgd = legend(strain);
    lgd.FontSize = 14;
    file1 = strcat([compareDir,'/','PTSG']);
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([compareDir,'/PTSG.fig']);
    savefig(gcf,file1_fig)
    
    figure(2)
    for i = 1:length(strains)
        plot(time{i},mean_green{i},'LineWidth',6)
        hold on
    end
    
    
    title(strcat([graph_title,' Green']),'FontSize',32)
    xlabel('Time after Induction (min)','FontSize',24)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 18)
    ylabel('ptsG Whole Cell Intensity','FontSize',24)
    set(gcf,'position',[835,883,868,667])
    lgd = legend(strain);
    lgd.FontSize = 14;
    file1 = strcat([compareDir,'/GREEN']);
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([compareDir,'GREEN.fig']);
    savefig(gcf,file1_fig)

elseif PTSG == 1 && sigma_check == 1
    figure(1)
    for i = 1:length(strains)
        errorbar(time{i},mean_ptsg{i},std_ptsg{i},'LineWidth',6)
        hold on
    end
    
    
    title(strcat([graph_title,' ptsG']),'FontSize',32)
    xlabel('Time after Induction (min)','FontSize',24)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 18)
    ylabel('ptsG Copy Number','FontSize',24)
    set(gcf,'position',[835,883,868,667])
    lgd = legend(strain);
    lgd.FontSize = 14;
    file1 = strcat([compareDir,'/','PTSG']);
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([compareDir,'/PTSG.fig']);
    savefig(gcf,file1_fig)
    
    figure(2)
    for i = 1:length(strains)
        plot(time{i},mean_green{i},std_green{i},'LineWidth',6)
        hold on
    end
    
    
    title(strcat([graph_title,' Green']),'FontSize',32)
    xlabel('Time after Induction (min)','FontSize',24)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 18)
    ylabel('ptsG Whole Cell Intensity','FontSize',24)
    set(gcf,'position',[835,883,868,667])
    lgd = legend(strain);
    lgd.FontSize = 14;
    file1 = strcat([compareDir,'/GREEN']);
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([compareDir,'GREEN.fig']);
    savefig(gcf,file1_fig)
    
end






