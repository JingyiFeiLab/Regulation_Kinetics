close all
clear all

blue_back_cell = 690.1429; %per pixel
figures = 0;

sub_red = 0;
sub_green = 0;

green_diff = 0;
blue_diff = 0;
red_diff = 0;


green_laser = [20];

time = [0,1,3,4,6,12,18,24];
strain1 = [156]; %Also, MR style
sample = {{[1,2,3],[1,2,3],[1,2,3],[1,2,3],[1,2,3],[1,2,3],[1,2,3],[1,2,3,4]}};
Date1 = {'May_30_2020'};
red_channel = 1;
green_channel = 2;
blue_channel = 3;
violet_channel = 4;



q = 0;
green_back = [];

mean_green = {};
mean_sgrs = {};
mean_red = {};
mean_violet = {};

std_green = {};
std_sgrs = {};
std_red = {};
std_violet = {};

for s_l = 1:length(strain1)
    strain_i = strcat(['MR',num2str(strain1(s_l))]);
    dateDir = strcat(['/Users/reyer/Data/SingleCellEpi/',strain_i,'/',Date1{s_l},'/']);
    
    green_back_cell = 11.6903*green_laser(s_l)+478.97;
    
    for j = 1:length(time)
        t = strcat(['t',num2str(time(j))]);
        close all
        
        timeDir = strcat([dateDir,t]);
        d2 = dir([timeDir, '/*.mat']);
        %samples = length(d2);
        samples = sample{s_l}{j};
        file = strcat([timeDir,'/']);
        if red_channel==1
            
            for k = samples
                red = [];
                red_back_pixel = 0;
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                red_back_image = (total_cells_one==0).*stack_one;
                red_back_image(red_back_image ==0) = [];
                red_back_pixel_intensity = median(red_back_image);
                for l = 1:length(part4)
                    red = [red (((part4(l).Intensity_One/(5*part4(l).Volume))-red_back_pixel_intensity - 1000)*5*part4(l).Volume)];
                end
                red(red == Inf) = [];
                mean_red{j,k}=mean(red);
                std_red{j,k} = std(red);
            end
            
        elseif red_channel == 2
            
            
            for k = samples
                red = [];
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                red_back_image = (total_cells_one==0).*stack_one;
                red_back_image(red_back_image ==0) = [];
                red_back_pixel_intensity = median(red_back_image);
                for l = 1:length(part4)
                    red = [red (((part4(l).Intensity_Two/(5*part4(l).Volume))-red_back_pixel_intensity - 1000)*5*part4(l).Volume)];
                end
                red(red == Inf) = [];
                mean_red{j,k}=mean(red);
                std_red{j,k} = std(red);
            end
            
            
           
        elseif red_channel == 3
            
            
            for k = samples
                red = [];
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                red_back_image = (total_cells_one==0).*stack_one;
                red_back_image(red_back_image ==0) = [];
                red_back_pixel_intensity = median(red_back_image);
                for l = 1:length(part4)
                    red = [red (((part4(l).Intensity_Three/(5*part4(l).Volume))-red_back_pixel_intensity - 1000)*5*part4(l).Volume)];
                end
                red(red == Inf) = [];
                mean_red{j,k}= mean(red);
                std_red{j,k} = std(red);
            end
            
            
        elseif red_channel == 4
            
            
            for k = samples
                red = [];
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                red_back_image = (total_cells_one==0).*stack_one;
                red_back_image(red_back_image ==0) = [];
                red_back_pixel_intensity = median(red_back_image);
                for l = 1:length(part4)
                    red = [red (((part4(l).Intensity_Four/(5*part4(l).Volume))-red_back_pixel_intensity - 1000)*5*part4(l).Volume)];
                end
                red(red == Inf) = [];
                mean_red{j,k}= mean(red);
                std_red{j,k} = std(red);
            end
            
        end
        
        if green_channel==1
            
            
            for k = samples
                green = [];
                q = q+1;
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                green_back_image = (total_cells_one==0).*stack_one;
                green_back_image(green_back_image ==0) = [];
                green_back_pixel_intensity = median(green_back_image);
                green_back(q) = green_back_pixel_intensity;
                for l = 1:length(part4)
                    green = [green (((part4(l).Intensity_One/(5*part4(l).Volume))-green_back_pixel_intensity - green_back_cell)*5*part4(l).Volume) ];
                end
                green(green == Inf) = [];
                mean_green{j,k}= mean(green);
                std_green{j,k} = std(green);
            end
            
        elseif green_channel == 2
            
            
            
            for k = samples
                green = [];
                q = q+1;
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                green_back_image = (total_cells_two==0).*stack_two;
                green_back_image(green_back_image ==0) = [];
                green_back_pixel_intensity = median(green_back_image);
                green_back(q) = green_back_pixel_intensity;
                for l = 1:length(part4)
                    green = [green (((part4(l).Intensity_Two/(5*part4(l).Volume))-green_back_pixel_intensity - green_back_cell)*5*part4(l).Volume) ];
                end
                green(green == Inf) = [];
                mean_green{j,k}= mean(green);
                std_green{j,k} = std(green);
            end
            
        elseif green_channel == 3
            
            
            
            for k = samples
                green = [];
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                green_back_image = (total_cells_three==0).*stack_three;
                green_back_image(green_back_image ==0) = [];
                green_back_pixel_intensity = median(green_back_image);
                green_back(q) = green_back_pixel_intensity;
                for l = 1:length(part4)
                    green = [green (((part4(l).Intensity_Three/(5*part4(l).Volume))-green_back_pixel_intensity - green_back_cell)*5*part4(l).Volume)];
                end
                green(green == Inf) = [];
                mean_green{j,k}= mean(green);
                std_green{j,k} = std(green);
            end
            
        elseif green_channel == 4
            
            for k = samples
                green = [];
                q = q+1;
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                green_back_image = (total_cells_one==0).*stack_one;
                green_back_image(green_back_image ==0) = [];
                green_back_pixel_intensity = median(green_back_image);
                green_back(q) = green_back_pixel_intensity;
                for l = 1:length(part4)
                    green = [green (((part4(l).Intensity_Four/(5*part4(l).Volume))-green_back_pixel_intensity - green_back_cell)*5*part4(l).Volume) ];
                end
                green(green == Inf) = [];
                mean_green{j,k}= mean(green);
                std_green{j,k} = std(green);
                
            end
            
        end
        
        if blue_channel==1
           
            for k = samples
                blue = [];
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                blue_back_image = (total_cells_one==0).*stack_one;
                blue_back_image(blue_back_image ==0) = [];
                blue_back_pixel_intensity = median(blue_back_image);
                
                for l = 1:length(part4)
                    blue = [blue (((part4(l).Intensity_One/(5*part4(l).Volume))-blue_back_pixel_intensity - blue_back_cell)*5*part4(l).Volume)];
                end
                blue(blue == Inf) = [];
                mean_blue{j,k}= mean(blue);
                std_blue{j,k} = std(blue);
            end
            
        elseif blue_channel == 2
            
            
            for k = samples
                blue = [];
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                blue_back_image = (total_cells_two==0).*stack_two;
                blue_back_image(blue_back_image ==0) = [];
                blue_back_pixel_intensity = median(blue_back_image);
                
                for l = 1:length(part4)
                    blue = [blue (((part4(l).Intensity_Two/(5*part4(l).Volume))-blue_back_pixel_intensity - blue_back_cell)*5*part4(l).Volume)];
                end
                blue(blue == Inf) = [];
                mean_blue{j,k}= mean(blue);
                std_blue{j,k} = std(blue);
            end
            
        elseif blue_channel == 3
            
            
            for k = samples
                blue = [];
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                blue_back_image = (total_cells_three==0).*stack_three;
                blue_back_image(blue_back_image ==0) = [];
                blue_back_pixel_intensity = median(blue_back_image);
                
                for l = 1:length(part4)
                    blue = [blue (((part4(l).Intensity_Three/(5*part4(l).Volume))-blue_back_pixel_intensity - blue_back_cell)*5*part4(l).Volume)];
                end
                blue(blue == Inf) = [];
                mean_blue{j,k}= mean(blue);
                std_blue{j,k} = std(blue);
            end
            
        elseif blue_channel == 4
            
            for k = samples
                blue = [];
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                blue_back_image = (total_cells_four==0).*stack_four;
                blue_back_image(blue_back_image ==0) = [];
                blue_back_pixel_intensity = median(blue_back_image);
                
                for l = 1:length(part4)
                    blue = [blue (((part4.Intensity_Four/(5*part4(l).Volume))-blue_back_pixel_intensity - blue_back_cell)*5*part4(l).Volume)];
                end
                blue(blue == Inf) = [];
                mean_blue{j,k}= mean(blue);
                std_blue{j,k} = std(blue);
            end
            
        end
        
        if violet_channel==1
           
            for k = samples
                violet = [];
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                violet_back_image = (total_cells_one==0).*stack_one;
                violet_back_image(violet_back_image ==0) = [];
                violet_back_pixel_intensity = median(violet_back_image);
                
                for l = 1:length(part4)
                    violet = [violet (((part4(l).Intensity_One/(5*part4(l).Volume))-violet_back_pixel_intensity)*5*part4(l).Volume)];
                end
                violet(violet == Inf) = [];
                mean_violet{j,k}= mean(violet);
                std_violet{j,k} = std(violet);
            end
            
        elseif violet_channel == 2
            
            
            for k = samples
                violet = [];
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                violet_back_image = (total_cells_two==0).*stack_two;
                violet_back_image(violet_back_image ==0) = [];
                violet_back_pixel_intensity = median(violet_back_image);
                
                for l = 1:length(part4)
                    violet = [violet (((part4(l).Intensity_Two/(5*part4(l).Volume))-violet_back_pixel_intensity)*5*part4(l).Volume)];
                end
                violet(violet == Inf) = [];
                mean_violet{j,k}= mean(violet);
                std_violet{j,k} = std(violet);
            end
            
        elseif violet_channel == 3
            
            
            for k = samples
                violet = [];
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                violet_back_image = (total_cells_three==0).*stack_three;
                violet_back_image(violet_back_image ==0) = [];
                violet_back_pixel_intensity = median(violet_back_image);
                
                for l = 1:length(part4)
                    violet = [violet (((part4(l).Intensity_Three/(5*part4(l).Volume))-violet_back_pixel_intensity)*5*part4(l).Volume)];
                end
                violet(violet == Inf) = [];
                mean_violet{j,k}= mean(violet);
                std_violet{j,k} = std(violet);
            end
            
        elseif violet_channel == 4
            
            for k = samples
                violet = [];
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                violet_back_image = (total_cells_four==0).*stack_four;
                violet_back_image(violet_back_image ==0) = [];
                violet_back_pixel_intensity = median(violet_back_image);
                
                for l = 1:length(part4)
                    violet = [violet (((part4(l).Intensity_Four/(5*part4(l).Volume))-violet_back_pixel_intensity)*5*part4(l).Volume)];
                end
                violet(violet == Inf) = [];
                mean_violet{j,k}= mean(violet);
                std_violet{j,k} = std(violet);
            end
            
        end
        
    end
end


for t = 1:length(time)
    samples = sample{s_l}{t};
    jt_mark = -1;
    offset = [];
    size_b = [];
    for jt = 0:length(samples)-1
        offset(jt+1) = time(t)+(.5*ceil(jt/2)*jt_mark);
        size_b(jt+1) = 100*(jt+1);
        jt_mark = jt_mark*-1;
    end
    
    blue_avs = [];
    green_avs = [];
    red_avs = [];
    violet_avs = [];
    
    blue_std = [];
    green_std = [];
    red_std = [];
    violet_std = [];
    
    kt_int = 0;
    for kt = samples
        kt_int = kt_int+1;
        blue_avs(kt_int) = mean_blue{t,kt};
        blue_std(kt_int) = std_blue{t,kt};
        green_avs(kt_int) = mean_green{t,kt};
        green_std(kt_int) = std_green{t,kt};
        if red_channel ~= 0
            red_avs(kt_int) = mean_red{t,kt};
            red_std(kt_int) = std_red{t,kt};
        end
        violet_avs(kt_int) = mean_violet{t,kt};
        violet_std(kt_int) = std_violet{t,kt};
        
    end
    
    if blue_channel == 0 && green_channel == 0
        for kt = samples
            kt_int = kt_int+1;
            red_avs(kt_int) = mean_red{t,kt};
            red_std(kt_int) = std_red{t,kt};
            violet_avs(kt_int) = mean_violet{t,kt};
            violet_std(kt_int) = std_violet{t,kt};
            
        end
    end
    
    if blue_channel~=0
        figure(1)
        scatter(offset,blue_avs,size_b,'b','filled')
        hold on
        errorbar(offset,blue_avs,blue_std,'b','LineStyle','none');
        hold on
    end
    
    if green_channel~=0
        figure(2)
        scatter(offset,green_avs,size_b,'g','filled')
        hold on
        errorbar(offset,green_avs,green_std,'g','LineStyle','none');
        hold on
    end
    
    if red_channel ~= 0
        figure(3)
        scatter(offset,red_avs,size_b,'r','filled')
        hold on
        errorbar(offset,red_avs,red_std,'r','LineStyle','none');
        hold on
    end
    
    figure(4)
    scatter(offset,violet_avs,size_b,'m','filled')
    hold on
    errorbar(offset,violet_avs,violet_std,'m','LineStyle','none');
    hold on
    
end

if blue_channel~=0
    figure(1)
    title(strcat([strain_i,' ',Date1,' GFP']),'FontSize',32)
    xlabel('Time after Induction (min)','FontSize',24)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 18)
    set(gca,'XLim',[-1 43])
    set(gca,'YLim',[0 20E6])
    ylabel('GFP Whole Cell Intensity','FontSize',24)
    set(gcf,'position',[835,883,868,667])
    file1 = strcat([dateDir,'/nvBlue_Samples']);
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([dateDir,'nvBlue_Samples.fig']);
    savefig(gcf,file1_fig)
end


if green_channel ~=0
    figure(2)
    title(strcat([strain_i,' ',Date1,' mRNA']),'FontSize',32)
    xlabel('Time after Induction (min)','FontSize',24)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 18)
    set(gca,'XLim',[-1 43])
    set(gca,'YLim',[0 3.5E6])
    ylabel('mRNA Whole Cell Intensity','FontSize',24)
    set(gcf,'position',[835,883,868,667])
    file1 = strcat([dateDir,'/nvGreen_Samples']);
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([dateDir,'nvGreen_Samples.fig']);
    savefig(gcf,file1_fig)
end

if red_channel ~= 0
    figure(3)
    title(strcat([strain_i,' ',Date1,' sRNA']),'FontSize',32)
    xlabel('Time after Induction (min)','FontSize',24)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 18)
    set(gca,'XLim',[-1 43])
    set(gca,'YLim',[0 3.5E6])
    ylabel('sRNA Whole Cell Intensity','FontSize',24)
    set(gcf,'position',[835,883,868,667])
    file1 = strcat([dateDir,'/nvRed_Samples']);
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([dateDir,'nvRed_Samples.fig']);
    savefig(gcf,file1_fig)
end

figure(4)
title(strcat([strain_i,' ',Date1,' rRNA']),'FontSize',32)
xlabel('Time after Induction (min)','FontSize',24)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 18)
set(gca,'XLim',[-1 43])
set(gca,'YLim',[0 1.8E6])
ylabel('rRNA Whole Cell Intensity','FontSize',24)
set(gcf,'position',[835,883,868,667])
file1 = strcat([dateDir,'/nvViolet_Samples']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([dateDir,'nvViolet_Samples.fig']);
savefig(gcf,file1_fig)

