parentDir = '/Users/reyer/Data/SingleCellEpi/';
strains = [156,153,157];
strain = {};
time = {};
Dates = {'April_12_2018','April_12_2018','April_12_2018'};
newFile = 'Plus_SgrS_mRNA_pre_Induced';
graph_title = '+SgrS, pre-induced mRNA,';

compareDir = strcat([parentDir,newFile]);
if exist(compareDir,'dir')~=7
    mkdir(parentDir, newFile)
end

GFP = 0;
SGRS = 1; %1 if you want to make SgrS comparisons
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

if GFP == 0 && PTSG == 1 && SGRS == 1 && sigma_check == 0
    close all
    figure(1)
    for i = 1:length(strains)
        
        
        scatter(time{i},mean_green{i},time{i},mean_red{i});
        cfig = get(gcf,'color');
        pos = [0.1  0.1  0.7  0.8];
        
        set(get(ax(1),'ylabel'),'string',ylabels{1})
        set(get(ax(2),'ylabel'),'string',ylabels{2})
        
        
        
        title(strcat([graph_title,' ',strain{i}]),'FontSize',32)
        xlabel('Time after Induction (min)','FontSize',24)
        ax(1).FontSize = 18;
        ax(2).FontSize = 18;
        
        set(gcf,'position',[835,883,868,667])
        lgd2 = legend('ptsG','SgrS');
        lgd2.FontSize = 14;
        file1 = strcat([compareDir,'/',strain{i},'_intensity']);
        set(gcf,'PaperPositionMode','auto')
        set(gcf,'PaperPositionMode','auto')
       
        close all
    
    end
end






