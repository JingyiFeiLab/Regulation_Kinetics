close all
clear all

blue_back_cell = 690.1429; %per pixel
figures = 0;
sub_red = 1;
sub_green = 0;

green_diff = 0;
blue_diff = 0;
red_diff = 0;

violet_norm = 1; % 1 if violet normalization


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

for s_l = 1:length(strain1)
    strain_i = strcat(['MR',num2str(strain1(s_l))]);
    dateDir = strcat(['/Users/reyer/Data/SingleCellEpi/',strain_i,'/',Date1{s_l},'/']);
    
    green_back_cell = 11.6903*green_laser(s_l)+478.97;
    red_back_cell = green_back_cell;
    
    
    
    
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
    
    se_SGRS = zeros(1,length(time));
    se_PTSG = zeros(1,length(time));
    se_GFP = zeros(1,length(time));
    se_RIB = zeros(1,length(time));
    se_RED = zeros(1,length(time));
    se_GREEN = zeros(1,length(time));
    se_BLUE = zeros(1,length(time));
    se_VIOLET = zeros(1,length(time));
    
    if violet_norm == 1
        for j = 1:length(time)
            t = strcat(['t',num2str(time(j))]);
            close all
            
            timeDir = strcat([dateDir,t]);
            d2 = dir([timeDir, '/*.mat']);
            %samples = length(d2);
            samples = sample{s_l}{j};
            file = strcat([timeDir,'/']);
            if red_channel==1
                red = [];
                red_back_pixel = 0;
                for k = samples
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'part4','part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    red_back_image = (total_cells_one==0).*stack_one;
                    red_back_image(red_back_image ==0) = [];
                    red_back_pixel_intensity = median(red_back_image);
                    for l = 1:length(part4_V_Normalized)
                        red = [red (((part4_V_Normalized(l).Original_One/(5*part4_V_Normalized(l).Volume))-red_back_pixel_intensity-red_back_cell)*5*part4_V_Normalized(l).Volume) + red_diff];
                    end
                end
                red(red<0) = 0;
                red(red == Inf) = [];
                
                sgrs = [];
                for m = 1:length(red)
                    sgrs = [sgrs cell2RNA(red(m),'sgrs')];
                end
            elseif red_channel == 2
                red = [];
                %samples = sample{s_l}{j};
                for k = samples
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    red_back_image = (total_cells_one==0).*stack_one;
                    red_back_image(red_back_image ==0) = [];
                    red_back_pixel_intensity = median(red_back_image);
                    for l = 1:length(part4_V_Normalized)
                        red = [red (((part4_V_Normalized(l).Original_Two/(5*part4_V_Normalized(l).Volume))-red_back_pixel_intensity-red_back_cell)*5*part4_V_Normalized(l).Volume) + red_diff];
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
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    red_back_image = (total_cells_one==0).*stack_one;
                    red_back_image(red_back_image ==0) = [];
                    red_back_pixel_intensity = median(red_back_image);
                    for l = 1:length(part4_V_Normalized)
                        red = [red (((part4_V_Normalized(l).Original_Three/(5*part4_V_Normalized(l).Volume))-red_back_pixel_intensity-red_back_cell)*5*part4_V_Normalized(l).Volume) + red_diff];
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
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    red_back_image = (total_cells_one==0).*stack_one;
                    red_back_image(red_back_image ==0) = [];
                    red_back_pixel_intensity = median(red_back_image);
                    for l = 1:length(part4_V_Normalized)
                        red = [red (((part4_V_Normalized(l).Original_Four/(5*part4_V_Normalized(l).Volume))-red_back_pixel_intensity-red_back_cell)*5*part4_V_Normalized(l).Volume) + red_diff];
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
                %samples = sample{s_l}{j};
                for k = samples
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
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    green_back_image = (total_cells_three==0).*stack_three;
                    green_back_image(green_back_image ==0) = [];
                    green_back_pixel_intensity = median(green_back_image);
                    green_back(q) = green_back_pixel_intensity;
                    for l = 1:length(part4_V_Normalized)
                        green = [green (((part4_V_Normalized(l).Original_Three/(5*part4_V_Normalized(l).Volume))-green_back_pixel_intensity - green_back_cell)*5*part4_V_Normalized(l).Volume) + green_diff];
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
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'part4','part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    blue_back_image = (total_cells_one==0).*stack_one;
                    blue_back_image(blue_back_image ==0) = [];
                    blue_back_pixel_intensity = median(blue_back_image);
                    
                    for l = 1:length(part4_V_Normalized)
                        blue = [blue (((part4_V_Normalized(l).Original_One/(5*part4_V_Normalized(l).Volume))-blue_back_pixel_intensity - blue_back_cell)*5*part4_V_Normalized(l).Volume) + blue_diff];
                    end
                end
                blue(blue<0)=0;
                
                gfp = [];
                for m = 1:length(blue)
                    gfp = [gfp cell2molecule(blue(m))];
                end
            elseif blue_channel == 2
                blue = [];
                
                %samples = sample{s_l}{j};
                for k = samples
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'part4','part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    blue_back_image = (total_cells_two==0).*stack_two;
                    blue_back_image(blue_back_image ==0) = [];
                    blue_back_pixel_intensity = median(blue_back_image);
                    
                    for l = 1:length(part4_V_Normalized)
                        blue = [blue (((part4_V_Normalized(l).Original_Two/(5*part4_V_Normalized(l).Volume))-blue_back_pixel_intensity - blue_back_cell)*5*part4_V_Normalized(l).Volume) + blue_diff];
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
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'part4','part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    blue_back_image = (total_cells_three==0).*stack_three;
                    blue_back_image(blue_back_image ==0) = [];
                    blue_back_pixel_intensity = median(blue_back_image);
                    
                    for l = 1:length(part4_V_Normalized)
                        blue = [blue (((part4_V_Normalized(l).Original_Three/(5*part4_V_Normalized(l).Volume))-blue_back_pixel_intensity - blue_back_cell)*5*part4_V_Normalized(l).Volume) + blue_diff];
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
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    blue_back_image = (total_cells_four==0).*stack_four;
                    blue_back_image(blue_back_image ==0) = [];
                    blue_back_pixel_intensity = median(blue_back_image);
                    
                    for l = 1:length(part4_V_Normalized)
                        blue = [blue (((part4_V_Normalized.Original_Four/(5*part4_V_Normalized(l).Volume))-blue_back_pixel_intensity - blue_back_cell)*5*part4_V_Normalized(l).Volume) + blue_diff];
                    end
                end
                blue(blue<0)=0;
                
                gfp = [];
                for m = 1:length(blue)
                    gfp = [gfp cell2molecule(blue(m))];
                end
            end
            
            if violet_channel==1
                violet_norm = [];
                samples = sample{s_l}{j};
                
                for k = samples
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    violet_back_image = (total_cells_one==0).*stack_one;
                    violet_back_image(violet_back_image ==0) = [];
                    violet_back_pixel_intensity = median(violet_back_image);
                    
                    for l = 1:length(part4_V_Normalized)
                        violet_norm = [violet_norm part4_V_Normalized(l).Original_One];
                    end
                end
                violet_norm(violet_norm<0)=0;
                
                
            elseif violet_channel == 2
                violet_norm = [];
                samples = sample{s_l}{j};
                
                for k = samples
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    violet_back_image = (total_cells_two==0).*stack_two;
                    violet_back_image(violet_back_image ==0) = [];
                    violet_back_pixel_intensity = median(violet_back_image);
                    
                    for l = 1:length(part4_V_Normalized)
                        violet_norm = [violet_norm part4_V_Normalized(l).Original_Two];
                    end
                end
                violet_norm(violet_norm<0)=0;
                
            elseif violet_channel == 3
                violet_norm = [];
                samples = sample{s_l}{j};
                
                for k = samples
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'part4','part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    violet_back_image = (total_cells_three==0).*stack_three;
                    violet_back_image(violet_back_image ==0) = [];
                    violet_back_pixel_intensity = median(violet_back_image);
                    
                    for l = 1:length(part4_V_Normalized)
                        violet_norm = [violet_norm part4_V_Normalized(l).Original_Three];
                        %violet = [violet (((part4_V_Normalized(l).Intensity_Three/(5*part4(l).Volume))-violet_back_pixel_intensity)*5*part4(l).Volume)];
                    end
                end
                violet_norm(violet_norm<0)=0;
                
            elseif violet_channel == 4
                violet_norm = [];
                samples = sample{s_l}{j};
                
                for k = samples
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'part4','part4_V_Normalized','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    violet_back_image = (total_cells_four==0).*stack_four;
                    violet_back_image(violet_back_image ==0) = [];
                    violet_back_pixel_intensity = median(violet_back_image);
                    
                    for l = 1:length(part4_V_Normalized)
                        violet_norm = [violet_norm part4_V_Normalized(l).Original_Four];
                    end
                end
                violet_norm(violet_norm<0)=0;
                
            end
            
            if isempty(sample{s_l}{j})
                continue
            end
            
            if blue_channel~=0
                
                l_blue = length(blue);
                [blue,blue_indices] = sort(blue);
%                 blue(l_blue-10:l_blue) = [];
%                 blue(1:10) = [];
%                 gfp(l_blue-10:l_blue) = [];
%                 gfp(1:10) = [];
                
                mean_BLUE(j) = mean(blue);
                mean_GFP(j) = mean(gfp);
                std_BLUE(j) = std(blue);
                std_GFP(j) = std(gfp);
                se_BLUE(j) = std(blue)/(sqrt(length(blue)));
                se_GFP(j) = std(gfp)/(sqrt(length(gfp)));
                
                if figures == 1
                    figure(1)
                    
                    w = blue;
                    histogram(w,20);hold on
                    mn = mean(w);    %%% Calculate the mean
                    stdv = std(w);     %%% Calculate the standard deviation
                    mnlabel = sprintf('Mean = %.2f', mn);
                    stdlabel = sprintf('Std Deviation = %.2f', stdv);
                    %%% Create the textbox
                    h = annotation('textbox',[0.68 0.75 0.1 0.1]);
                    set(h,'String',{mnlabel,stdlabel},'FontSize',20, 'FontWeight', 'bold');
                    title(strcat([Date1{1}, ' GFP Intensity, t=',num2str(time(j))]),'FontSize',28,'FontWeight','bold','Interpreter', 'none')
                    xlabel('Whole Cell Intensity','FontSize',24)
                    xt = get(gca, 'XTick');
                    set(gca, 'FontSize', 18)
                    set(gca,'XLim',[0,1E7])
                    set(gcf,'position',[835,883,868,667])
                    file1 = strcat([file,'ThreshBlue']);
                    set(gcf,'PaperPositionMode','auto')
                    print(file1,'-dpng','-r0')
                    print(file1,'-painters','-depsc','-r0')
                    file1_fig = strcat([file,'ThreshBlue.fig']);
                    savefig(gcf,file1_fig)
                    
                end
                
            end
            close all
            
            if green_channel ~=0
                
                
                
                l_green = length(green);
                if blue_channel~=0
                    green = green(blue_indices);
%                     green(l_green-10:l_green) = [];
%                     green(1:10) = [];
%                     ptsg(l_green-10:l_green) = [];
%                     ptsg(1:10) = [];
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
                se_GREEN(j) = std(green)/(sqrt(length(green)));
                se_PTSG(j) = std(ptsg)/(sqrt(length(ptsg)));
                
                
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
                    set(h,'String',{mnlabel,stdlabel},'FontSize',20, 'FontWeight', 'bold');
                    title(strcat([Date1{1}, ' ptsG Intensity, t=',num2str(time(j))]),'FontSize',28,'FontWeight','bold','Interpreter', 'none')
                    xlabel('Whole Cell Intensity','FontSize',24)
                    xt = get(gca, 'XTick');
                    set(gca, 'FontSize', 18)
                    set(gca,'XLim',[0,5E6])
                    set(gca, 'FontSize', 18)
                    set(gcf,'position',[835,883,868,667])
                    file1 = strcat([file,'ThreshGreen']);
                    set(gcf,'PaperPositionMode','auto')
                    print(file1,'-dpng','-r0')
                    print(file1,'-painters','-depsc','-r0')
                    file1_fig = strcat([file,'ThreshGreen.fig']);
                    savefig(gcf,file1_fig)
                    
                end
                
            end
            close all
            
            if red_channel ~=0
                
                l_red = length(red);
                if blue_channel~=0
%                     red = red(blue_indices);
%                     red(l_red-10:l_red) = [];
%                     red(1:10) = [];
%                     sgrs = sgrs(blue_indices);
%                     sgrs(l_red-10:l_red) = [];
%                     sgrs(1:10) = [];
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
                se_RED(j) = std(red)/(sqrt(length(red)));
                se_SGRS(j) = std(sgrs)/(sqrt(length(sgrs)));
                
                
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
                    set(h,'String',{mnlabel,stdlabel},'FontSize',20, 'FontWeight', 'bold');
                    title(strcat([Date1{1}, ' SgrS Intensity, t=',num2str(time(j))]),'FontSize',28,'FontWeight','bold','Interpreter', 'none')
                    xlabel('Whole Cell Intensity','FontSize',24)
                    xt = get(gca, 'XTick');
                    set(gca,'XLim',[0,5E6])
                    set(gca, 'FontSize', 18)
                    set(gcf,'position',[835,883,868,667])
                    file1 = strcat([file,'ThreshRed']);
                    set(gcf,'PaperPositionMode','auto')
                    print(file1,'-dpng','-r0')
                    print(file1,'-painters','-depsc','-r0')
                    file1_fig = strcat([file,'ThreshRed.fig']);
                    savefig(gcf,file1_fig)
                    
                    
                end
                
            end
            close all
            
        end
        
    elseif violet_norm == 0
        for j = 1:length(time)
            t = strcat(['t',num2str(time(j))]);
            close all
            
            timeDir = strcat([dateDir,t]);
            d2 = dir([timeDir, '/*.mat']);
            %samples = length(d2);
            samples = sample{s_l}{j};
            file = strcat([timeDir,'/']);
            if red_channel==1
                red = [];
                red_back_pixel = 0;
                for k = samples
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    red_back_image = (total_cells_one==0).*stack_one;
                    red_back_image(red_back_image ==0) = [];
                    red_back_pixel_intensity = median(red_back_image);
                    for l = 1:length(part4)
                        red = [red (((part4(l).Intensity_One/(5*part4(l).Volume))-red_back_pixel_intensity - 1000)*5*part4(l).Volume)];
                    end
                end
                red(red<0) = 0;
                
                
                sgrs = [];
                for m = 1:length(red)
                    sgrs = [sgrs cell2RNA(red(m),'sgrs')];
                end
            elseif red_channel == 2
                red = [];
                %samples = sample{s_l}{j};
                for k = samples
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    red_back_image = (total_cells_one==0).*stack_one;
                    red_back_image(red_back_image ==0) = [];
                    red_back_pixel_intensity = median(red_back_image);
                    for l = 1:length(part4)
                        red = [red (((part4(l).Intensity_Two/(5*part4(l).Volume))-red_back_pixel_intensity - 1000)*5*part4(l).Volume)];
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
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    red_back_image = (total_cells_one==0).*stack_one;
                    red_back_image(red_back_image ==0) = [];
                    red_back_pixel_intensity = median(red_back_image);
                    for l = 1:length(part4)
                        red = [red (((part4(l).Intensity_Three/(5*part4(l).Volume))-red_back_pixel_intensity - 1000)*5*part4(l).Volume)];
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
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    red_back_image = (total_cells_one==0).*stack_one;
                    red_back_image(red_back_image ==0) = [];
                    red_back_pixel_intensity = median(red_back_image);
                    for l = 1:length(part4)
                        red = [red (((part4(l).Intensity_Four/(5*part4(l).Volume))-red_back_pixel_intensity - 1000)*5*part4(l).Volume)];
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
                %samples = sample{s_l}{j};
                for k = samples
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
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    green_back_image = (total_cells_three==0).*stack_three;
                    green_back_image(green_back_image ==0) = [];
                    green_back_pixel_intensity = median(green_back_image);
                    green_back(q) = green_back_pixel_intensity;
                    for l = 1:length(part4)
                        green = [green (((part4(l).Intensity_Three/(5*part4(l).Volume))-green_back_pixel_intensity - green_back_cell)*5*part4(l).Volume)];
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
                    q = q+1;
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    green_back_image = (total_cells_one==0).*stack_one;
                    green_back_image(green_back_image ==0) = [];
                    green_back_pixel_intensity = median(green_back_image);
                    green_back(q) = green_back_pixel_intensity;
                    for l = 1:length(part4)
                        green = [green (((part4(l).Intensity_Four/(5*part4(l).Volume))-green_back_pixel_intensity - green_back_cell)*5*part4(l).Volume) ];
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
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    blue_back_image = (total_cells_one==0).*stack_one;
                    blue_back_image(blue_back_image ==0) = [];
                    blue_back_pixel_intensity = median(blue_back_image);
                    
                    for l = 1:length(part4)
                        blue = [blue (((part4(l).Intensity_One/(5*part4(l).Volume))-blue_back_pixel_intensity - blue_back_cell)*5*part4(l).Volume)];
                    end
                end
                blue(blue<0)=0;
                
                gfp = [];
                for m = 1:length(blue)
                    gfp = [gfp cell2molecule(blue(m))];
                end
            elseif blue_channel == 2
                blue = [];
                
                %samples = sample{s_l}{j};
                for k = samples
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    blue_back_image = (total_cells_two==0).*stack_two;
                    blue_back_image(blue_back_image ==0) = [];
                    blue_back_pixel_intensity = median(blue_back_image);
                    
                    for l = 1:length(part4)
                        blue = [blue (((part4(l).Intensity_Two/(5*part4(l).Volume))-blue_back_pixel_intensity - blue_back_cell)*5*part4(l).Volume)];
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
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    blue_back_image = (total_cells_three==0).*stack_three;
                    blue_back_image(blue_back_image ==0) = [];
                    blue_back_pixel_intensity = median(blue_back_image);
                    
                    for l = 1:length(part4)
                        blue = [blue (((part4(l).Intensity_Three/(5*part4(l).Volume))-blue_back_pixel_intensity - blue_back_cell)*5*part4(l).Volume)];
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
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    blue_back_image = (total_cells_four==0).*stack_four;
                    blue_back_image(blue_back_image ==0) = [];
                    blue_back_pixel_intensity = median(blue_back_image);
                    
                    for l = 1:length(part4)
                        blue = [blue (((part4.Intensity_Four/(5*part4(l).Volume))-blue_back_pixel_intensity - blue_back_cell)*5*part4(l).Volume)];
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
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    violet_back_image = (total_cells_one==0).*stack_one;
                    violet_back_image(violet_back_image ==0) = [];
                    violet_back_pixel_intensity = median(violet_back_image);
                    
                    for l = 1:length(part4)
                        violet = [violet part4(l).Intensity_One];
                    end
                end
                violet(violet<0)=0;
                
                
            elseif violet_channel == 2
                violet = [];
                samples = sample{s_l}{j};
                
                for k = samples
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    violet_back_image = (total_cells_two==0).*stack_two;
                    violet_back_image(violet_back_image ==0) = [];
                    violet_back_pixel_intensity = median(violet_back_image);
                    
                    for l = 1:length(part4)
                        violet = [violet part4(l).Intensity_Two];
                    end
                end
                violet(violet<0)=0;
                
            elseif violet_channel == 3
                violet = [];
                samples = sample{s_l}{j};
                
                for k = samples
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    violet_back_image = (total_cells_three==0).*stack_three;
                    violet_back_image(violet_back_image ==0) = [];
                    violet_back_pixel_intensity = median(violet_back_image);
                    
                    for l = 1:length(part4)
                        violet = [violet part4(l).Intensity_Three];
                        %violet = [violet (((part4(l).Intensity_Three/(5*part4(l).Volume))-violet_back_pixel_intensity)*5*part4(l).Volume)];
                    end
                end
                violet(violet<0)=0;
                
            elseif violet_channel == 4
                violet = [];
                samples = sample{s_l}{j};
                
                for k = samples
                    s = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
                    load(s,'part4','total_cells_one','stack_one','total_cells_two','stack_two','total_cells_three','stack_three','total_cells_four','stack_four')
                    violet_back_image = (total_cells_four==0).*stack_four;
                    violet_back_image(violet_back_image ==0) = [];
                    violet_back_pixel_intensity = median(violet_back_image);
                    
                    for l = 1:length(part4)
                        violet = [violet part4(l).Intensity_Four];
                    end
                end
                violet(violet<0)=0;
                
            end
            
            if isempty(sample{s_l}{j})
                continue
            end
            
            if blue_channel~=0
                
                l_blue = length(blue);
                [blue,blue_indices] = sort(blue);
%                 blue(l_blue-10:l_blue) = [];
%                 blue(1:10) = [];
%                 gfp(l_blue-10:l_blue) = [];
%                 gfp(1:10) = [];
                
                mean_BLUE(j) = mean(blue);
                mean_GFP(j) = mean(gfp);
                std_BLUE(j) = std(blue);
                std_GFP(j) = std(gfp);
                se_BLUE(j) = std(blue)/(sqrt(length(blue)));
                se_GFP(j) = std(gfp)/(sqrt(length(gfp)));
                
                if figures == 1
                    figure(1)
                    
                    w = blue;
                    histogram(w,20);hold on
                    mn = mean(w);    %%% Calculate the mean
                    stdv = std(w);     %%% Calculate the standard deviation
                    mnlabel = sprintf('Mean = %.2f', mn);
                    stdlabel = sprintf('Std Deviation = %.2f', stdv);
                    %%% Create the textbox
                    h = annotation('textbox',[0.68 0.75 0.1 0.1]);
                    set(h,'String',{mnlabel,stdlabel},'FontSize',20, 'FontWeight', 'bold');
                    title(strcat([Date1{1}, ' GFP Intensity, t=',num2str(time(j))]),'FontSize',28,'FontWeight','bold','Interpreter', 'none')
                    xlabel('Whole Cell Intensity','FontSize',24)
                    xt = get(gca, 'XTick');
                    set(gca, 'FontSize', 18)
                    set(gca,'XLim',[0,1E7])
                    set(gcf,'position',[835,883,868,667])
                    file1 = strcat([file,'ThreshBlue']);
                    set(gcf,'PaperPositionMode','auto')
                    print(file1,'-dpng','-r0')
                    print(file1,'-painters','-depsc','-r0')
                    file1_fig = strcat([file,'ThreshBlue.fig']);
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
                se_GREEN(j) = std(green)/(sqrt(length(green)));
                se_PTSG(j) = std(ptsg)/(sqrt(length(ptsg)));
                
                
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
                    set(h,'String',{mnlabel,stdlabel},'FontSize',20, 'FontWeight', 'bold');
                    title(strcat([Date1{1}, ' mRNA Intensity, t=',num2str(time(j))]),'FontSize',28,'FontWeight','bold','Interpreter', 'none')
                    xlabel('Whole Cell Intensity','FontSize',24)
                    xt = get(gca, 'XTick');
                    set(gca, 'FontSize', 18)
                    set(gca,'XLim',[0,5E6])
                    set(gca, 'FontSize', 18)
                    set(gcf,'position',[835,883,868,667])
                    file1 = strcat([file,'ThreshGreen']);
                    set(gcf,'PaperPositionMode','auto')
                    print(file1,'-dpng','-r0')
                    print(file1,'-painters','-depsc','-r0')
                    file1_fig = strcat([file,'ThreshGreen.fig']);
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
                se_RED(j) = std(red)/(sqrt(length(red)));
                se_SGRS(j) = std(sgrs)/(sqrt(length(sgrs)));
                
                
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
                    set(h,'String',{mnlabel,stdlabel},'FontSize',20, 'FontWeight', 'bold');
                    title(strcat([Date1{1}, ' SgrS Intensity, t=',num2str(time(j))]),'FontSize',28,'FontWeight','bold','Interpreter', 'none')
                    xlabel('Whole Cell Intensity','FontSize',24)
                    xt = get(gca, 'XTick');
                    set(gca,'XLim',[0,5E6])
                    set(gca, 'FontSize', 18)
                    set(gcf,'position',[835,883,868,667])
                    file1 = strcat([file,'ThreshRed']);
                    set(gcf,'PaperPositionMode','auto')
                    print(file1,'-dpng','-r0')
                    print(file1,'-painters','-depsc','-r0')
                    file1_fig = strcat([file,'ThreshRed.fig']);
                    savefig(gcf,file1_fig)
                    
                    
                end
                
            end
            close all
            
        end
    end
    
    copy_array = [time',mean_GFP',std_GFP',mean_PTSG',std_PTSG',mean_SGRS',std_SGRS'];
    int_array = [time',mean_BLUE',std_BLUE',se_BLUE',mean_GREEN',std_GREEN',se_GREEN',mean_RED',std_RED',se_RED',mean_VIOLET',std_VIOLET'];
    
    copy_table = array2table(copy_array);copy_table.Properties.VariableNames = {'Time' 'Mean_GFP' 'GFP_Sigma' 'Mean_ptsG' 'ptsG_Sigma' 'Mean_SgrS' 'SgrS_Sigma'};
    int_table = array2table(int_array);int_table.Properties.VariableNames = {'Time' 'Mean_Blue' 'Blue_Sigma' 'Blue_SE' 'Mean_Green' 'Green_Sigma' 'Green_SE' 'Mean_Red' 'Red_Sigma' 'Red_SE' 'Mean_Violet' 'Violet_Sigma'};
    
    copy_table_file = strcat([dateDir,'VThresh_BackSubCopy_numbers.csv']);
    int_table_file = strcat([dateDir,'VThresh_BackSubIntensity.csv']);
    
    writetable(copy_table,copy_table_file);
    writetable(int_table,int_table_file);
    
    dataSaveFile = strcat([dateDir,'VThresh_BackgroundSummaryData.mat']);
    save(dataSaveFile,'copy_table','int_table');
        
    
end




