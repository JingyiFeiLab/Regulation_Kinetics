close all
clear all

blue_back_cell = 690.1429; %per pixel
figures = 1;

sub_red = 1;
sub_green = 0;

green_diff = 0;
blue_diff = 0;
red_diff = 0;


green_laser = [20];

time = [0,6,12,18,24];
strain1 = [192]; %Also, MR style
sample = {{[1,2],[1,2],[1,2,3],[1,2,3],[1,2,3]}};
Date1 = {'July_18_2020_1'};
red_channel = 1;
green_channel = 0;
blue_channel = 0;
violet_channel = 2;

q = 0;
green_back = [];

mean_green = {};
mean_sgrs = {};
mean_red = {};

std_green = {};
std_sgrs = {};
std_red = {};

for s_l = 1:length(strain1)
    strain_i = strcat(['MR',num2str(strain1(s_l))]);
    dateDir = strcat(['/Users/reyer/Data/SingleCellEpi/',strain_i,'/',Date1{s_l},'/']);
    
    green_back_cell = 11.6903*green_laser(s_l)+478.97;
    red_back_cell = green_back_cell;
    
    green_zero = [];
    blue_zero = [];
    red_zero = [];
    
    for j = 1
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
                load(s,'part4','part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                red_back_image = (total_cells_one==0).*stack_one;
                red_back_image(red_back_image ==0) = [];
                red_back_pixel_intensity = median(red_back_image);
                for l = 1:length(part4_V_Normalized)
                    red_zero = [red_zero (((part4_V_Normalized(l).Original_One/(5*part4_V_Normalized(l).Volume))-red_back_pixel_intensity-red_back_cell)*5*part4_V_Normalized(l).Volume) ];
                end
                red_zero(red_zero == Inf) = [];
            end
            
        elseif red_channel == 2
            
            
            for k = samples
                red = [];
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                red_back_image = (total_cells_one==0).*stack_one;
                red_back_image(red_back_image ==0) = [];
                red_back_pixel_intensity = median(red_back_image);
                for l = 1:length(part4_V_Normalized)
                    red_zero = [red_zero (((part4_V_Normalized(l).Original_Two/(5*part4_V_Normalized(l).Volume))-red_back_pixel_intensity-red_back_cell)*5*part4_V_Normalized(l).Volume) ];
                end
                red_zero(red_zero == Inf) = [];
            end
            
           
        elseif red_channel == 3
            
            
            for k = samples
                red = [];
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                red_back_image = (total_cells_one==0).*stack_one;
                red_back_image(red_back_image ==0) = [];
                red_back_pixel_intensity = median(red_back_image);
                for l = 1:length(part4_V_Normalized)
                    red_zero = [red_zero (((part4_V_Normalized(l).Original_Three/(5*part4_V_Normalized(l).Volume))-red_back_pixel_intensity-red_back_cell)*5*part4_V_Normalized(l).Volume) ];
                end
                red_zero(red_zero == Inf) = [];
            end
            
            
        elseif red_channel == 4
            
            
            for k = samples
                red = [];
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                red_back_image = (total_cells_one==0).*stack_one;
                red_back_image(red_back_image ==0) = [];
                red_back_pixel_intensity = median(red_back_image);
                for l = 1:length(part4_V_Normalized)
                    red_zero = [red_zero (((part4_V_Normalized(l).Original_Four/(5*part4_V_Normalized(l).Volume))-red_back_pixel_intensity-red_back_cell)*5*part4_V_Normalized(l).Volume) ];
                end
                red_zero(red_zero == Inf) = [];
            end
            
        end
        
        
        if green_channel==1
            
            
            for k = samples
                
                q = q+1;
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                green_back_image = (total_cells_one==0).*stack_one;
                green_back_image(green_back_image ==0) = [];
                green_back_pixel_intensity = median(green_back_image);
                
                green_back(q) = green_back_pixel_intensity;
                for l = 1:length(part4_V_Normalized)
                    green_zero = [green_zero (((part4_V_Normalized(l).Original_One/(5*part4_V_Normalized(l).Volume))-green_back_pixel_intensity - green_back_cell)*5*part4_V_Normalized(l).Volume) ];
                end
                green_zero(green_zero == Inf) = [];
                
            end
            
        elseif green_channel == 2
            
            
            
            for k = samples
                if j == 5 && k == 3
                    k
                end
                
                q = q+1;
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                green_back_image = (total_cells_two==0).*stack_two;
                green_back_image(green_back_image ==0) = [];
                green_back_pixel_intensity = median(green_back_image);
                green_back(q) = green_back_pixel_intensity;
                for l = 1:length(part4_V_Normalized)
                    green_zero = [green_zero (((part4_V_Normalized(l).Original_Two/(5*part4_V_Normalized(l).Volume))-green_back_pixel_intensity - green_back_cell)*5*part4_V_Normalized(l).Volume) ];
                end
                green_zero(green_zero == Inf) = [];
                
            end
            
        elseif green_channel == 3
            
            
            
            for k = samples
                
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                green_back_image = (total_cells_three==0).*stack_three;
                green_back_image(green_back_image ==0) = [];
                green_back_pixel_intensity = median(green_back_image);
                green_back(q) = green_back_pixel_intensity;
                for l = 1:length(part4_V_Normalized)
                    green_zero = [green_zero (((part4_V_Normalized(l).Original_Three/(5*part4_V_Normalized(l).Volume))-green_back_pixel_intensity - green_back_cell)*5*part4_V_Normalized(l).Volume)];
                end
                green_zero(green_zero == Inf) = [];
                
            end
            
        elseif green_channel == 4
            
            for k = samples
                
                q = q+1;
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                green_back_image = (total_cells_one==0).*stack_one;
                green_back_image(green_back_image ==0) = [];
                green_back_pixel_intensity = median(green_back_image);
                green_back(q) = green_back_pixel_intensity;
                for l = 1:length(part4_V_Normalized)
                    green_zero = [green_zero (((part4_V_Normalized(l).Original_Four/(5*part4_V_Normalized(l).Volume))-green_back_pixel_intensity - green_back_cell)*5*part4_V_Normalized(l).Volume) ];
                end
                green_zero(green_zero == Inf) = [];
                
                
            end
            
        end
        
        if blue_channel==1
           
            for k = samples
                
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                blue_back_image = (total_cells_one==0).*stack_one;
                blue_back_image(blue_back_image ==0) = [];
                blue_back_pixel_intensity = median(blue_back_image);
                
                for l = 1:length(part4_V_Normalized)
                    blue_zero = [blue_zero (((part4_V_Normalized(l).Original_One/(5*part4_V_Normalized(l).Volume))-blue_back_pixel_intensity - blue_back_cell)*5*part4_V_Normalized(l).Volume)];
                end
                blue_zero(blue_zero == Inf) = [];
                
            end
            
        elseif blue_channel == 2
            
            
            for k = samples
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                blue_back_image = (total_cells_two==0).*stack_two;
                blue_back_image(blue_back_image ==0) = [];
                blue_back_pixel_intensity = median(blue_back_image);
                
                for l = 1:length(part4_V_Normalized)
                    blue_zero = [blue_zero (((part4_V_Normalized(l).Original_Two/(5*part4_V_Normalized(l).Volume))-blue_back_pixel_intensity - blue_back_cell)*5*part4_V_Normalized(l).Volume)];
                end
                blue_zero(blue_zero == Inf) = [];
               
            end
            
        elseif blue_channel == 3
            
            
            for k = samples
                
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                blue_back_image = (total_cells_three==0).*stack_three;
                blue_back_image(blue_back_image ==0) = [];
                blue_back_pixel_intensity = median(blue_back_image);
                
                for l = 1:length(part4_V_Normalized)
                    blue_zero = [blue_zero (((part4_V_Normalized(l).Original_Three/(5*part4_V_Normalized(l).Volume))-blue_back_pixel_intensity - blue_back_cell)*5*part4_V_Normalized(l).Volume)];
                end
                blue_zero(blue_zero == Inf) = [];
                
            end
            
        elseif blue_channel == 4
            
            for k = samples
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                blue_back_image = (total_cells_four==0).*stack_four;
                blue_back_image(blue_back_image ==0) = [];
                blue_back_pixel_intensity = median(blue_back_image);
                
                for l = 1:length(part4_V_Normalized)
                    blue_zero = [blue_zero (((part4_V_Normalized.Original_Four/(5*part4_V_Normalized(l).Volume))-blue_back_pixel_intensity - blue_back_cell)*5*part4_V_Normalized(l).Volume)];
                end
                blue_zero(blue_zero == Inf) = [];
               
            end
            
        end
        
    end
end


if sub_green == 1
    green_diff = 1000-mean(green_zero);
    blue_diff = 1000-mean(blue_zero);
    red_diff = 0;
end

if sub_red == 1
    green_diff = 0;
    blue_diff = 0;
    red_diff = 1000-mean(red_zero);
end

if sub_red == 1 && sub_green == 1
    green_diff = 1000-mean(green_zero);
    blue_diff = 1000-mean(blue_zero);
    red_diff = 1000-mean(red_zero);
end
    

for s_l = 1:length(strain1)
    strain_i = strcat(['MR',num2str(strain1(s_l))]);
    dateDir = strcat(['/Users/reyer/Data/SingleCellEpi/',strain_i,'/',Date1{s_l},'/']);
    
    green_back_cell = 11.6903*green_laser(s_l)+478.97;
    red_back_cell = green_back_cell;
    
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
                load(s,'part4','part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                red_back_image = (total_cells_one==0).*stack_one;
                red_back_image(red_back_image ==0) = [];
                red_back_pixel_intensity = median(red_back_image);
                for l = 1:length(part4_V_Normalized)
                    red = [red (((part4_V_Normalized(l).Original_One/(5*part4_V_Normalized(l).Volume))-red_back_pixel_intensity-red_back_cell)*5*part4_V_Normalized(l).Volume) + red_diff];
                end
                red(red == Inf) = [];
                mean_red{j,k}=mean(red);
                std_red{j,k} = std(red);
            end
            
        elseif red_channel == 2
            
            
            for k = samples
                red = [];
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                red_back_image = (total_cells_one==0).*stack_one;
                red_back_image(red_back_image ==0) = [];
                red_back_pixel_intensity = median(red_back_image);
                for l = 1:length(part4_V_Normalized)
                    red = [red (((part4_V_Normalized(l).Original_Two/(5*part4_V_Normalized(l).Volume))-red_back_pixel_intensity-red_back_cell)*5*part4_V_Normalized(l).Volume) + red_diff];
                end
                red(red == Inf) = [];
                mean_red{j,k}=mean(red);
                std_red{j,k} = std(red);
            end
            
           
        elseif red_channel == 3
            
            
            for k = samples
                red = [];
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                red_back_image = (total_cells_one==0).*stack_one;
                red_back_image(red_back_image ==0) = [];
                red_back_pixel_intensity = median(red_back_image);
                for l = 1:length(part4_V_Normalized)
                    red = [red (((part4_V_Normalized(l).Original_Three/(5*part4_V_Normalized(l).Volume))-red_back_pixel_intensity-red_back_cell)*5*part4_V_Normalized(l).Volume) + red_diff];
                end
                red(red == Inf) = [];
                mean_red{j,k}= mean(red);
                std_red{j,k} = std(red);
            end
            
            
        elseif red_channel == 4
            
            
            for k = samples
                red = [];
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                red_back_image = (total_cells_one==0).*stack_one;
                red_back_image(red_back_image ==0) = [];
                red_back_pixel_intensity = median(red_back_image);
                for l = 1:length(part4_V_Normalized)
                    red = [red (((part4_V_Normalized(l).Original_Four/(5*part4_V_Normalized(l).Volume))-red_back_pixel_intensity-red_back_cell)*5*part4_V_Normalized(l).Volume) + red_diff];
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
                load(s,'part4','part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                green_back_image = (total_cells_one==0).*stack_one;
                green_back_image(green_back_image ==0) = [];
                green_back_pixel_intensity = median(green_back_image);
                
                green_back(q) = green_back_pixel_intensity;
                for l = 1:length(part4_V_Normalized)
                    green = [green (((part4_V_Normalized(l).Original_One/(5*part4_V_Normalized(l).Volume))-green_back_pixel_intensity - green_back_cell)*5*part4_V_Normalized(l).Volume) + green_diff ];
                end
                green(green == Inf) = [];
                mean_green{j,k}= mean(green);
                std_green{j,k} = std(green);
            end
            
        elseif green_channel == 2
            
            
            
            for k = samples
                if j == 5 && k == 3
                    k
                end
                green = [];
                q = q+1;
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                green_back_image = (total_cells_two==0).*stack_two;
                green_back_image(green_back_image ==0) = [];
                green_back_pixel_intensity = median(green_back_image);
                green_back(q) = green_back_pixel_intensity;
                for l = 1:length(part4_V_Normalized)
                    green = [green (((part4_V_Normalized(l).Original_Two/(5*part4_V_Normalized(l).Volume))-green_back_pixel_intensity - green_back_cell)*5*part4_V_Normalized(l).Volume) + green_diff ];
                end
                green(green == Inf) = [];
                mean_green{j,k}= mean(green);
                std_green{j,k} = std(green);
            end
            
        elseif green_channel == 3
            
            
            
            for k = samples
                green = [];
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                green_back_image = (total_cells_three==0).*stack_three;
                green_back_image(green_back_image ==0) = [];
                green_back_pixel_intensity = median(green_back_image);
                green_back(q) = green_back_pixel_intensity;
                for l = 1:length(part4_V_Normalized)
                    green = [green (((part4_V_Normalized(l).Original_Three/(5*part4_V_Normalized(l).Volume))-green_back_pixel_intensity - green_back_cell)*5*part4_V_Normalized(l).Volume) + green_diff];
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
                load(s,'part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                green_back_image = (total_cells_one==0).*stack_one;
                green_back_image(green_back_image ==0) = [];
                green_back_pixel_intensity = median(green_back_image);
                green_back(q) = green_back_pixel_intensity;
                for l = 1:length(part4_V_Normalized)
                    green = [green (((part4_V_Normalized(l).Original_Four/(5*part4_V_Normalized(l).Volume))-green_back_pixel_intensity - green_back_cell)*5*part4_V_Normalized(l).Volume) + green_diff ];
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
                
                for l = 1:length(part4_V_Normalized)
                    blue = [blue (((part4_V_Normalized(l).Original_One/(5*part4_V_Normalized(l).Volume))-blue_back_pixel_intensity - blue_back_cell)*5*part4_V_Normalized(l).Volume) + blue_diff];
                end
                blue(blue == Inf) = [];
                mean_blue{j,k}= mean(blue);
                std_blue{j,k} = std(blue);
            end
            
        elseif blue_channel == 2
            
            
            for k = samples
                blue = [];
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                blue_back_image = (total_cells_two==0).*stack_two;
                blue_back_image(blue_back_image ==0) = [];
                blue_back_pixel_intensity = median(blue_back_image);
                
                for l = 1:length(part4_V_Normalized)
                    blue = [blue (((part4_V_Normalized(l).Original_Two/(5*part4_V_Normalized(l).Volume))-blue_back_pixel_intensity - blue_back_cell)*5*part4_V_Normalized(l).Volume) + blue_diff];
                end
                blue(blue == Inf) = [];
                mean_blue{j,k}= mean(blue);
                std_blue{j,k} = std(blue);
            end
            
        elseif blue_channel == 3
            
            
            for k = samples
                blue = [];
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4','part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                blue_back_image = (total_cells_three==0).*stack_three;
                blue_back_image(blue_back_image ==0) = [];
                blue_back_pixel_intensity = median(blue_back_image);
                
                for l = 1:length(part4_V_Normalized)
                    blue = [blue (((part4_V_Normalized(l).Original_Three/(5*part4_V_Normalized(l).Volume))-blue_back_pixel_intensity - blue_back_cell)*5*part4_V_Normalized(l).Volume) + blue_diff];
                end
                blue(blue == Inf) = [];
                mean_blue{j,k}= mean(blue);
                std_blue{j,k} = std(blue);
            end
            
        elseif blue_channel == 4
            
            for k = samples
                blue = [];
                s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                load(s,'part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                blue_back_image = (total_cells_four==0).*stack_four;
                blue_back_image(blue_back_image ==0) = [];
                blue_back_pixel_intensity = median(blue_back_image);
                
                for l = 1:length(part4_V_Normalized)
                    blue = [blue (((part4_V_Normalized.Original_Four/(5*part4_V_Normalized(l).Volume))-blue_back_pixel_intensity - blue_back_cell)*5*part4_V_Normalized(l).Volume) + blue_diff];
                end
                blue(blue == Inf) = [];
                mean_blue{j,k}= mean(blue);
                std_blue{j,k} = std(blue);
            end
            
        end
        
    end
end


for t = 1:length(time)
    if t == 5
        t
    end
    samples = sample{s_l}{t};
    jt_mark = -1;
    offset = [];
    size_b = [];
    for jt = 0:length(samples)-1
        offset(jt+1) = time(t)+(.34*ceil(jt/2)*jt_mark);
        size_b(jt+1) = 100*(jt+1);
        jt_mark = jt_mark*-1;
    end
    
    blue_avs = [];
    green_avs = [];
    red_avs = [];
    
    blue_std = [];
    green_std = [];
    red_std = [];
    
    kt_int = 0;
    
    if blue_channel ~=0 && green_channel~=0
        for kt = samples
            kt_int = kt_int+1;
            blue_avs(kt_int) = mean_blue{t,kt};
            blue_std(kt_int) = std_blue{t,kt};
            green_avs(kt_int) = mean_green{t,kt};
            green_std(kt_int) = std_green{t,kt};
            if red_channel ~=0
                red_avs(kt_int) = mean_red{t,kt};
                red_std(kt_int) = std_red{t,kt};
            end
            
        end
    else
        
        for kt = samples
            kt_int = kt_int+1;
            if red_channel ~=0
                red_avs(kt_int) = mean_red{t,kt};
                red_std(kt_int) = std_red{t,kt};
            end
            
        end
    end
    
    if blue_channel ~=0 && green_channel~=0
        figure(1)
        scatter(offset,blue_avs,size_b,'b','filled')
        hold on
        errorbar(offset,blue_avs,blue_std/2,'b','LineStyle','none');
        hold on
        
        figure(2)
        scatter(offset,green_avs,size_b,'g','filled')
        hold on
        errorbar(offset,green_avs,green_std/2,'g','LineStyle','none');
        hold on
    end
    
    if red_channel ~= 0
        figure(3)
        scatter(offset,red_avs,size_b,'r','filled')
        hold on
        errorbar(offset,red_avs,red_std/2,'r','LineStyle','none');
        hold on
    end
    
    
end

if blue_channel ~=0 && green_channel~=0
    figure(1)
    title(strcat([strain_i,' ',Date1,' GFP']),'FontSize',32)
    xlabel('Time after Induction (min)','FontSize',24)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 18)
    ylabel('GFP Whole Cell Intensity','FontSize',24)
    set(gca,'XLim',[-1 26])
    set(gca,'YLim',[0 15E6])
    %set(gca,'FontSize',18)
    set(gcf,'position',[835,883,868,667])
    file1 = strcat([dateDir,'/Original_Blue_Samples']);
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([dateDir,'Original_Blue_Samples.fig']);
    savefig(gcf,file1_fig)
    
    figure(2)
    title(strcat([strain_i,' ',Date1,' mRNA']),'FontSize',32)
    xlabel('Time after Induction (min)','FontSize',24)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 18)
    ylabel('mRNA Whole Cell Intensity','FontSize',24)
    set(gca,'XLim',[-1 26])
    %set(gca,'YLim',[0 4.5E6])
    set(gca,'YLim',[0 55E5])
    set(gca,'FontSize',18)
    set(gcf,'position',[835,883,868,667])
    file1 = strcat([dateDir,'/Original_Green_Samples']);
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([dateDir,'Original_Green_Samples.fig']);
    savefig(gcf,file1_fig)
end

if red_channel ~=0
    figure(3)
    title(strcat([strain_i,' ',Date1,' sRNA']),'FontSize',32)
    xlabel('Time after Induction (min)','FontSize',24)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 18)
    ylabel('sRNA Whole Cell Intensity','FontSize',24)
    set(gca,'XLim',[-1 26])
    set(gca,'YLim',[0 5.5E6])
    set(gca,'FontSize',18)
    set(gcf,'position',[835,883,868,667])
    file1 = strcat([dateDir,'/Original_Red_Samples']);
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([dateDir,'Original_Red_Samples.fig']);
    savefig(gcf,file1_fig)
end
