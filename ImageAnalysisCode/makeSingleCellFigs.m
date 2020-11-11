close all
clear all
time = [0,10,20];
blue_back = 472590;
figures = 1;

strain = [60,156]; %Also, MR style
green_laser = [20];

sample = {{[3],[1,2,3],[1,2,3]},{[1,2],[1,2,3],[1,2]}};
Date = {'April_21_2018_10'};

for s_l = 1:length(strain)
    strain_i = strcat(['MR',num2str(strain(s_l))]);
    dateDir = strcat(['/Users/reyer/Data/SingleCellEpi/',strain_i,'/',Date{s_l},'/']);
    red_channel = 1;
    green_channel = 2;
    blue_channel = 3;
    violet_channel = 4;
    green_back = 6334.84*green_laser(s_l) + 367401.27;
    
    
    mean_SGRS = zeros(1,length(time));
    mean_PTSG = zeros(1,length(time));
    mean_GFP = zeros(1,length(time));
    mean_RIB = zeros(1,length(time));
    mean_RED = zeros(1,length(time));
    mean_GREEN = zeros(1,length(time));
    mean_BLUE = zeros(1,length(time));
    mean_VIOLET = zeros(1,length(time));
    
    std_SGRS = zeros(1,length(time));
    std_PTSG = zeros(1,length(time));
    std_GFP = zeros(1,length(time));
    std_RIB = zeros(1,length(time));
    std_RED = zeros(1,length(time));
    std_GREEN = zeros(1,length(time));
    std_BLUE = zeros(1,length(time));
    std_VIOLET = zeros(1,length(time));
    
    for j = 1:length(time)
        t = strcat(['t',num2str(time(j))]);
        close all
        d = strcat([dateDir,t]);
        d2 = dir([d, '/*.mat']);
        %samples = length(d2);
        samples = sample{s_l}{j};
        file = strcat([d,'/']);
        if red_channel==1
            red = [];
            red_back_pixel = 0;
            for k = samples
                s = strcat([d,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                red_back_image = ((total_cells_one==0).*stack_one);
                for l = 1:length(part4)
                    red = [red part4(l).Intensity_One - 740000];
                end
            end
            red(red<0) = 0;
            
            sgrs = [];
            for m = 1:length(red)
                sgrs = [sgrs cell2RNA(red(m),'sgrs')];
            end
        elseif red_channel == 2
            red = [];
            samples = sample{s_l}{j};
            for k = samples
                s = strcat([d,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                for l = 1:length(part4)
                    red = [red part4(l).Intensity_Two - 740000];
                end
            end
            red(red<0) = 0;
            
            sgrs = [];
            for m = 1:length(red)
                sgrs = [sgrs cell2RNA(red(m),'sgrs')];
            end
        elseif red_channel == 3
            red = [];
            samples = sample{s_l}{j};
            for k = samples
                s = strcat([d,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                for l = 1:length(part4)
                    red = [red part4(l).Intensity_Three - 740000];
                end
            end
            red(red<0) = 0;
            
            sgrs = [];
            for m = 1:length(red)
                sgrs = [sgrs cell2RNA(red(m),'sgrs')];
            end
            
        elseif red_channel == 4
            red = [];
            samples = sample{s_l}{j};
            for k = samples
                s = strcat([d,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                for l = 1:length(part4)
                    red = [red part4(l).Intensity_Four - 740000];
                end
            end
            red(red<0) = 0;
            
            sgrs = [];
            for m = 1:length(red)
                sgrs = [sgrs cell2RNA(red(m),'sgrs')];
            end
        end
        
        if green_channel==1
            green = [];
            samples = sample{s_l}{j};
            for k = samples
                s = strcat([d,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                for l = 1:length(part4)
                    green = [green part4(l).Intensity_One - green_back];
                end
            end
            green(green<0) = 0;
            
            ptsg = [];
            for m = 1:length(green)
                ptsg = [ptsg cell2RNA(green(m),'ptsg')];
            end
        elseif green_channel == 2
            green = [];
            samples = sample{s_l}{j};
            
            for k = samples
                s = strcat([d,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                for l = 1:length(part4)
                    green = [green part4(l).Intensity_Two ];
                end
            end
            green(green<0) = 0;
            
            ptsg = [];
            for m = 1:length(green)
                ptsg = [ptsg cell2RNA(green(m),'ptsg')];
            end
        elseif green_channel == 3
            green = [];
            samples = sample{s_l}{j};
            
            for k = samples
                s = strcat([d,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                for l = 1:length(part4)
                    green = [green part4(l).Intensity_Three - green_back];
                end
            end
            green(green<0) = 0;
            
            ptsg = [];
            for m = 1:length(green)
                ptsg = [ptsg cell2RNA(green(m),'ptsg')];
            end
            
        elseif green_channel == 4
            green = [];
            samples = sample{s_l}{j};
            
            for k = samples
                s = strcat([d,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                for l = 1:length(part4)
                    green = [green part4(l).Intensity_Four ];
                end
            end
            green(green<0) = 0;
            
            ptsg = [];
            for m = 1:length(green)
                ptsg = [ptsg cell2RNA(green(m),'ptsg')];
            end
        end
        
        if blue_channel==1
            blue = [];
            samples = sample{s_l}{j};
            
            for k = samples
                s = strcat([d,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                for l = 1:length(part4)
                    blue = [blue part4(l).Intensity_One - blue_back];
                end
            end
            blue(blue<0)=0;
            
            gfp = [];
            for m = 1:length(blue)
                gfp = [gfp cell2molecule(blue(m))];
            end
        elseif blue_channel == 2
            blue = [];
            
            samples = sample{s_l}{j};
            for k = samples
                s = strcat([d,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                for l = 1:length(part4)
                    blue = [blue part4(l).Intensity_Two - blue_back];
                end
            end
            blue(blue<0)=0;
            
            gfp = [];
            for m = 1:length(blue)
                gfp = [gfp cell2molecule(blue(m))];
            end
        elseif blue_channel == 3
            blue = [];
            samples = sample{s_l}{j};
            
            for k = samples
                s = strcat([d,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                for l = 1:length(part4)
                    blue = [blue part4(l).Intensity_Three - blue_back];
                end
            end
            blue(blue<0)=0;
            
            gfp = [];
            for m = 1:length(blue)
                gfp = [gfp cell2molecule(blue(m))];
            end
            
        elseif blue_channel == 4
            blue = [];
            
            samples = sample{s_l}{j};
            for k = samples
                s = strcat([d,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                for l = 1:length(part4)
                    blue = [blue part4(l).Intensity_Four - blue_back];
                end
            end
            blue(blue<0)=0;
            
            gfp = [];
            for m = 1:length(blue)
                gfp = [gfp cell2molecule(blue(m))];
            end
        end
        
        if violet_channel==1
            violet = [];
            samples = sample{s_l}{j};
            
            for k = samples
                s = strcat([d,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                for l = 1:length(part4)
                    violet = [violet part4(l).Intensity_One ];
                end
            end
            violet(violet<0)=0;
            
            
        elseif violet_channel == 2
            violet = [];
            
            samples = sample{s_l}{j};
            for k = samples
                s = strcat([d,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                for l = 1:length(part4)
                    violet = [violet part4(l).Intensity_Two ];
                end
            end
            violet(violet<0)=0;
            
            
        elseif violet_channel == 3
            violet = [];
            samples = sample{s_l}{j};
            
            for k = samples
                s = strcat([d,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                for l = 1:length(part4)
                    violet = [violet part4(l).Intensity_Three ];
                end
            end
            violet(violet<0)=0;
           
            
        elseif violet_channel == 4
            violet = [];
            
            samples = sample{s_l}{j};
            for k = samples
                s = strcat([d,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                for l = 1:length(part4)
                    violet = [violet part4(l).Intensity_Four];
                end
            end
            violet(violet<0)=0;
            
            
        end
        
        if blue_channel~=0
            
            l_blue = length(blue);
            [blue,blue_indices] = sort(blue);
            blue(l_blue-10:l_blue) = [];
            blue(1:10) = [];
            gfp(l_blue-10:l_blue) = [];
            gfp(1:10) = [];
            
            mean_BLUE(j) = mean(blue);
            mean_GFP(j) = mean(gfp);
            std_BLUE(j) = std(blue);
            std_GFP(j) = std(gfp);
            
            if figures == 1
                figure(1)
                w = blue;
                histogram(w);hold on
                mn = mean(w);    %%% Calculate the mean
                stdv = std(w);     %%% Calculate the standard deviation
                mnlabel = sprintf('Mean = %.2f', mn);
                stdlabel = sprintf('Std Deviation = %.2f', stdv);
                %%% Create the textbox
                h = annotation('textbox',[0.68 0.75 0.1 0.1]);
                set(h,'String',{mnlabel,stdlabel},'FontSize',16);
                title(strcat([strain_i, ' GFP Intensity, t=',num2str(time(j))]),'FontSize',28,'FontWeight','bold')
                xlabel('Whole Cell Intensity','FontSize',24)
                xt = get(gca, 'XTick');
                set(gca, 'FontSize', 18)
                set(gcf,'position',[835,883,868,667])
                file1 = strcat([file,'blue']);
                set(gcf,'PaperPositionMode','auto')
                print(file1,'-dpng','-r0')
                print(file1,'-painters','-depsc','-r0')
                file1_fig = strcat([file,'blue.fig']);
                savefig(gcf,file1_fig)
                
                figure(2)
                w = gfp;
                histogram(w);hold on
                mn = mean(w);    %%% Calculate the mean
                stdv = std(w);     %%% Calculate the standard deviation
                mnlabel = sprintf('Mean = %.2f', mn);
                stdlabel = sprintf('Std Deviation = %.2f', stdv);
                %%% Create the textbox
                h = annotation('textbox',[0.68 0.75 0.1 0.1]);
                set(h,'String',{mnlabel,stdlabel},'FontSize',16);
                title(strcat([strain_i, ' GFP Copy Number, t=',num2str(time(j))]),'FontSize',28,'FontWeight','bold')
                xlabel('GFP Copy Number','FontSize',24)
                xt = get(gca, 'XTick');
                set(gca, 'FontSize', 18)
                set(gcf,'PaperPositionMode','auto')
                print(file1,'-dpng','-r0')
                set(gcf,'position',[835,883,868,667])
                file1 = strcat([file,'gfp']);
                set(gcf,'PaperPositionMode','auto')
                print(file1,'-dpng','-r0')
                print(file1,'-painters','-depsc','-r0')
                file1_fig = strcat([file,'gfp.fig']);
                savefig(gcf,file1_fig)
            end
            
        end
        close all
        
        if green_channel ~=0
            
            l_green = length(green);
            if blue_channel~=0
                green = green(blue_indices);
                green(l_green-10:l_green) = [];
                green(1:10) = [];
                ptsg(l_green-10:l_green) = [];
                ptsg(1:10) = [];
            else
                [green,green_indices] = sort(green);
                green(l_green-10:l_green) = [];
                green(1:10) = [];
                ptsg(l_green-10:l_green) = [];
                ptsg(1:10) = [];
            end
            
            mean_GREEN(j) = mean(green);
            mean_PTSG(j) = mean(ptsg);
            std_GREEN(j) = std(green);
            std_PTSG(j) = std(ptsg);
            
            if figures == 1
                figure(1)
                w = green;
                histogram(w);hold on
                mn = mean(w);    %%% Calculate the mean
                stdv = std(w);     %%% Calculate the standard deviation
                mnlabel = sprintf('Mean = %.2f', mn);
                stdlabel = sprintf('Std Deviation = %.2f', stdv);
                %%% Create the textbox
                h = annotation('textbox',[0.68 0.75 0.1 0.1]);
                set(h,'String',{mnlabel,stdlabel},'FontSize',16);
                title(strcat([strain_i, ' ptsG Intensity, t=',num2str(time(j))]),'FontSize',28,'FontWeight','bold')
                xlabel('Whole Cell Intensity','FontSize',24)
                xt = get(gca, 'XTick');
                set(gca, 'FontSize', 18)
                set(gcf,'position',[835,883,868,667])
                file1 = strcat([file,'green']);
                set(gcf,'PaperPositionMode','auto')
                print(file1,'-dpng','-r0')
                print(file1,'-painters','-depsc','-r0')
                file1_fig = strcat([file,'green.fig']);
                savefig(gcf,file1_fig)
                
                figure(2)
                w = ptsg;
                histogram(w);hold on
                mn = mean(w);    %%% Calculate the mean
                stdv = std(w);     %%% Calculate the standard deviation
                mnlabel = sprintf('Mean = %.2f', mn);
                stdlabel = sprintf('Std Deviation = %.2f', stdv);
                %%% Create the textbox
                h = annotation('textbox',[0.68 0.75 0.1 0.1]);
                set(h,'String',{mnlabel,stdlabel},'FontSize',16);
                title(strcat([strain_i, ' ptsG Copy Number, t=',num2str(time(j))]),'FontSize',28,'FontWeight','bold')
                xlabel('ptsG Copy Number','FontSize',24)
                xt = get(gca, 'XTick');
                set(gca, 'FontSize', 18)
                set(gcf,'position',[835,883,868,667])
                file1 = strcat([file,'ptsg']);
                set(gcf,'PaperPositionMode','auto')
                print(file1,'-dpng','-r0')
                print(file1,'-painters','-depsc','-r0')
                file1_fig = strcat([file,'ptsg.fig']);
                savefig(gcf,file1_fig)
            end
            
        end
        close all
        
        if red_channel ~=0
            
            l_red = length(red);
            if blue_channel~=0
                red = red(blue_indices);
                red(l_red-10:l_red) = [];
                red(1:10) = [];
                sgrs = sgrs(blue_indices);
                sgrs(l_red-10:l_red) = [];
                sgrs(1:10) = [];
            elseif blue_channel == 0 && green_channel ~=0
                red = red(green_indices);
                red(l_red-10:l_red) = [];
                red(1:10) = [];
                sgrs = red(green_indices);
                sgrs(l_red-10:l_red) = [];
                sgrs(1:10) = [];
            else
                [red,red_indices] = sort(red);
                red(l_red-10:l_red) = [];
                red(1:10) = [];
                sgrs(l_red-10:l_red) = [];
                sgrs(1:10) = [];
            end
            
            mean_RED(j) = mean(red);
            mean_SGRS(j) = mean(sgrs);
            std_RED(j) = std(red);
            std_SGRS(j) = std(sgrs);
            
            if figures == 1
                figure(1)
                w = red;
                histogram(w);hold on
                mn = mean(w);    %%% Calculate the mean
                stdv = std(w);     %%% Calculate the standard deviation
                mnlabel = sprintf('Mean = %.2f', mn);
                stdlabel = sprintf('Std Deviation = %.2f', stdv);
                %%% Create the textbox
                h = annotation('textbox',[0.68 0.75 0.1 0.1]);
                set(h,'String',{mnlabel,stdlabel},'FontSize',16);
                title(strcat([strain_i, ' SgrS Intensity, t=',num2str(time(j))]),'FontSize',28,'FontWeight','bold')
                xlabel('Whole Cell Intensity','FontSize',24)
                xt = get(gca, 'XTick');
                set(gca, 'FontSize', 18)
                set(gcf,'position',[835,883,868,667])
                file1 = strcat([file,'red']);
                set(gcf,'PaperPositionMode','auto')
                print(file1,'-dpng','-r0')
                print(file1,'-painters','-depsc','-r0')
                file1_fig = strcat([file,'red.fig']);
                savefig(gcf,file1_fig)
                
                figure(2)
                w = sgrs;
                histogram(w);hold on
                mn = mean(w);    %%% Calculate the mean
                stdv = std(w);     %%% Calculate the standard deviation
                mnlabel = sprintf('Mean = %.2f', mn);
                stdlabel = sprintf('Std Deviation = %.2f', stdv);
                %%% Create the textbox
                h = annotation('textbox',[0.68 0.75 0.1 0.1]);
                set(h,'String',{mnlabel,stdlabel},'FontSize',16);
                title(strcat([strain_i, ' SgrS Copy Number, t=',num2str(time(j))]),'FontSize',28,'FontWeight','bold')
                xlabel('SgrS Copy Number','FontSize',24)
                xt = get(gca, 'XTick');
                set(gca, 'FontSize', 18)
                set(gcf,'position',[835,883,868,667])
                file1 = strcat([file,'sgrs']);
                set(gcf,'PaperPositionMode','auto')
                print(file1,'-dpng','-r0')
                print(file1,'-painters','-depsc','-r0')
                file1_fig = strcat([file,'sgrs.fig']);
                savefig(gcf,file1_fig)
            end
            
        end
        
        if violet_channel ~=0
            
            l_violet = length(violet);
            if blue_channel~=0
                violet = violet(blue_indices);
                violet(l_violet-10:l_violet) = [];
                violet(1:10) = [];
                
            elseif blue_channel == 0 && green_channel ~=0
                violet = violet(green_indices);
                violet(l_violet-10:l_violet) = [];
                violet(1:10) = [];
            elseif blue_channel == 0 && green_channel == 0 && red_channel ~=0
                violet = violet(red_indices); 
                violet(l_violet-10:l_violet) = [];
                violet(1:10) = [];
            else
                [violet,violet_indices] = sort(violet);
                violet(l_violet-10:l_violet) = [];
                violet(1:10) = [];
                
            end
            
            mean_VIOLET(j) = mean(violet);
            std_VIOLET(j) = std(violet);
            
            
            if figures == 1
                figure(1)
                w = violet;
                histogram(w);hold on
                mn = mean(w);    %%% Calculate the mean
                stdv = std(w);     %%% Calculate the standard deviation
                mnlabel = sprintf('Mean = %.2f', mn);
                stdlabel = sprintf('Std Deviation = %.2f', stdv);
                %%% Create the textbox
                h = annotation('textbox',[0.68 0.75 0.1 0.1]);
                set(h,'String',{mnlabel,stdlabel},'FontSize',16);
                title(strcat([strain_i, ' rRNA Intensity, t=',num2str(time(j))]),'FontSize',28,'FontWeight','bold')
                xlabel('Whole Cell Intensity','FontSize',24)
                xt = get(gca, 'XTick');
                set(gca, 'FontSize', 18)
                set(gcf,'position',[835,883,868,667])
                file1 = strcat([file,'violet']);
                set(gcf,'PaperPositionMode','auto')
                print(file1,'-dpng','-r0')
                print(file1,'-painters','-depsc','-r0')
                file1_fig = strcat([file,'violet.fig']);
                savefig(gcf,file1_fig)
            end
            
        end
        close all
        
        
        
        
        
    end
    
    copy_array = [time',mean_GFP',std_GFP',mean_PTSG',std_PTSG',mean_SGRS',std_SGRS'];
    int_array = [time',mean_BLUE',std_BLUE',mean_GREEN',std_GREEN',mean_RED',std_RED',mean_VIOLET',std_VIOLET'];
    
    copy_table = array2table(copy_array);copy_table.Properties.VariableNames = {'Time' 'Mean_GFP' 'GFP_Sigma' 'Mean_ptsG' 'ptsG_Sigma' 'Mean_SgrS' 'SgrS_Sigma'};
    int_table = array2table(int_array);int_table.Properties.VariableNames = {'Time' 'Mean_Blue' 'Blue_Sigma' 'Mean_Green' 'Green_Sigma' 'Mean_Red' 'Red_Sigma','Mean_Violet','Violet_Sigma'};
    
    copy_table_file = strcat([dateDir,'copy_numbers.csv']);
    int_table_file = strcat([dateDir,'intensity.csv']);
    
    writetable(copy_table,copy_table_file);
    writetable(int_table,int_table_file);
    
    dataSaveFile = strcat([dateDir,'summaryData.mat']);
    save(dataSaveFile,'copy_table','int_table');
        
    
end




