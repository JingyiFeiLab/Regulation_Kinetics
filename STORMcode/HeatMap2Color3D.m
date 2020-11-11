clear cell_center_x cell_center_z cell_color cell_region
close all

%% End length of each bin in microns. User-Determined
% You can also determine how many bins there should be, e.g.
% micron_bins = [0,3,5,7,10] -> gives 4 bins, binning cells from 0 to 3
% microns in length, 3-5 microns in length, 5-7, and 7-10

parentFile = '/Users/reyer/Data/STORM';

%dataSet = 'SP98_NT';
dataSet = 'SP_2018_11_29';
dataFile = strcat([parentFile,'/',dataSet]);
dates = {'11_29'}; 
micron_bins = [0,3,5,7,10,12];
samples_per_date = [1];
spot_cutoff = .85;

end_lengths = micron_bins/.130; % Convert lengths to pixels

for bin = 1:length(micron_bins)-1
    cells_in_bin = [];
    widths_in_bin = [];
    lengths_in_bin = [];
    
    
    red_cell_center_x = [];
    red_cell_center_z = [];
    green_cell_center_x = [];
    green_cell_center_z = [];
    
    red_cell_center_x_dapi = [];
    red_cell_center_z_dapi = [];
    green_cell_center_x_dapi = [];
    green_cell_center_z_dapi = [];
    
    red_cell_center_x_membrane = [];
    red_cell_center_z_membrane = [];
    green_cell_center_x_membrane = [];
    green_cell_center_z_membrane = [];
    
    red_cell_center_x_pole = [];
    red_cell_center_z_pole = [];
    green_cell_center_x_pole = [];
    green_cell_center_z_pole = [];
    
    red_cell_center_x_cyto = [];
    red_cell_center_z_cyto = [];
    green_cell_center_x_cyto = [];
    green_cell_center_z_cyto = [];
    
    cell_diffusion = [];
    cell_color = [];
    cell_region = [];
    
    for id = 1:length(dates)
        for ip = 0:samples_per_date(id)-1
            %spaceFile = strcat([dataFile,'/workspace_2018_',dates{id},'_ptsG_mMaple3_CMptsG_mEOS3Tr',num2str(ip),'.mat']);
            %spaceFile = strcat([dataFile,'/workspace_2018_',dates{id},'_SP98_CM_RifSP98_NT_TR',num2str(ip),'.mat']);
            %load(spaceFile)
            for n = 1:length(cell_struct)
                if cell_struct(n).Cell_Y_Axis <= end_lengths(bin+1) && cell_struct(n).Cell_Y_Axis > end_lengths(bin) && ~isempty(cell_struct(n).Red_Spots) && ~isempty(cell_struct(n).Green_Spots)
                    widths_in_bin = [widths_in_bin cell_struct(n).Cell_X_Axis];
                end
            end
        end
    end
    
    for id = 1:length(dates)
        for ip = 1
        %for ip = 0:samples_per_date(id)-1
            %spaceFile = strcat([dataFile,'/workspace_2018_',dates{id},'_ptsG_mMaple3_CMptsG_mEOS3Tr',num2str(ip),'.mat']);
            %spaceFile = strcat([dataFile,'/workspace_2018_',dates{id},'_SP98_CM_RifSP98_NT_TR',num2str(ip),'.mat']);
            %load(spaceFile)
            for n = 1:length(cell_struct)
                if cell_struct(n).Cell_Y_Axis <= end_lengths(bin+1) && cell_struct(n).Cell_Y_Axis > end_lengths(bin) && ~isempty(cell_struct(n).Red_Spots) && ~isempty(cell_struct(n).Green_Spots)
                    cells_in_bin = [cells_in_bin n];
                    lengths_in_bin = [lengths_in_bin cell_struct(n).Cell_Y_Axis];
                    
                    if isempty(cell_struct(n).Red_Spots) || isempty(cell_struct(n).Green_Spots)
                        continue
                    end
                    
                    length_normalization = end_lengths(bin+1)/cell_struct(n).Cell_Y_Axis;
                    width_normalization = mean(widths_in_bin)/cell_struct(n).Cell_X_Axis;
                    
                    for si = 1:length(cell_struct(n).Red_Spots(:,1))
                        
                        red_cell_center_x = [red_cell_center_x cell_struct(n).Red_Spots(si,2)*width_normalization];
                        red_cell_center_z = [red_cell_center_z cell_struct(n).Red_Spots(si,4)*width_normalization];
                        
                        if cell_struct(n).Red_Spots(si,5) == 1
                            red_cell_center_x_membrane = [red_cell_center_x_membrane cell_struct(n).Red_Spots(si,2)*width_normalization];
                            red_cell_center_z_membrane = [red_cell_center_z_membrane cell_struct(n).Red_Spots(si,4)*width_normalization];
                        elseif cell_struct(n).Red_Spots(si,5) == 2
                            red_cell_center_x_cyto = [red_cell_center_x_cyto cell_struct(n).Red_Spots(si,2)*width_normalization];
                            red_cell_center_z_cyto = [red_cell_center_z_cyto cell_struct(n).Red_Spots(si,4)*width_normalization];
                        elseif cell_struct(n).Red_Spots(si,5) == 3
                            red_cell_center_x_pole = [red_cell_center_x_pole cell_struct(n).Red_Spots(si,2)*width_normalization];
                            red_cell_center_z_pole = [red_cell_center_z_pole cell_struct(n).Red_Spots(si,4)*width_normalization];
                        elseif cell_struct(n).Red_Spots(si,5) == 4
                            red_cell_center_x_dapi = [red_cell_center_x_dapi cell_struct(n).Red_Spots(si,2)*width_normalization];
                            red_cell_center_z_dapi = [red_cell_center_z_dapi cell_struct(n).Red_Spots(si,4)*width_normalization];
                        end
                        
                    end
                    
                    for sgi = 1:length(cell_struct(n).Green_Spots(:,1))
                        
                        green_cell_center_x = [green_cell_center_x cell_struct(n).Green_Spots(sgi,2)*width_normalization];
                        green_cell_center_z = [green_cell_center_z cell_struct(n).Green_Spots(sgi,4)*width_normalization];
                        
                        if cell_struct(n).Green_Spots(sgi,5) == 1
                            green_cell_center_x_membrane = [green_cell_center_x_membrane cell_struct(n).Green_Spots(sgi,2)*width_normalization];
                            green_cell_center_z_membrane = [green_cell_center_z_membrane cell_struct(n).Green_Spots(sgi,4)*width_normalization];
                        elseif cell_struct(n).Green_Spots(sgi,5) == 2
                            green_cell_center_x_cyto = [green_cell_center_x_cyto cell_struct(n).Green_Spots(sgi,2)*width_normalization];
                            green_cell_center_z_cyto = [green_cell_center_z_cyto cell_struct(n).Green_Spots(sgi,4)*width_normalization];
                        elseif cell_struct(n).Green_Spots(sgi,5) == 3
                            green_cell_center_x_pole = [green_cell_center_x_pole cell_struct(n).Green_Spots(sgi,2)*width_normalization];
                            green_cell_center_z_pole = [green_cell_center_z_pole cell_struct(n).Green_Spots(sgi,4)*width_normalization];
                        elseif cell_struct(n).Green_Spots(sgi,5) == 4
                            green_cell_center_x_dapi = [green_cell_center_x_dapi cell_struct(n).Green_Spots(sgi,2)*width_normalization];
                            green_cell_center_z_dapi = [green_cell_center_z_dapi cell_struct(n).Green_Spots(sgi,4)*width_normalization];
                        end
                        
                    end
                    
                end
                
            end
        end
    end
    
    if isempty(cells_in_bin)
       
        continue
    end
    
    red_x_range = int32(min(red_cell_center_x):max(red_cell_center_x));
    red_z_range = int32(min(red_cell_center_z):max(red_cell_center_z));
    
    green_x_range = int32(min(green_cell_center_x):max(green_cell_center_x));
    green_z_range = int32(min(green_cell_center_z):max(green_cell_center_z));
    
    %binned_diffusion = zeros(length(red_y_range),length(red_x_range));
    binned_red_spots = zeros(length(red_z_range),length(red_x_range));
    binned_red_spots_membrane = zeros(length(red_z_range),length(red_x_range));
    binned_red_spots_dapi = zeros(length(red_z_range),length(red_x_range));
    binned_red_spots_cyto = zeros(length(red_z_range),length(red_x_range));
    binned_red_spots_pole = zeros(length(red_z_range),length(red_x_range));
    
    for ix = 1:length(red_x_range)
        for iy = 1:length(red_z_range)
            temp_red_spots = [];
            for isp = 1:length(red_cell_center_x)
                if ix == length(red_x_range) && iy < length(red_z_range)
                    if red_cell_center_x(isp) > red_x_range(ix) && red_cell_center_z(isp) > red_z_range(iy) && red_cell_center_z(isp) <= red_z_range(iy+1)
                        temp_red_spots = [temp_red_spots isp];
                    end
                elseif ix < length(red_x_range) && iy == length(red_z_range)
                    if red_cell_center_x(isp) > red_x_range(ix) && red_cell_center_x(isp) <= red_x_range(ix+1) && red_cell_center_z(isp) > red_z_range(iy) 
                        temp_red_spots = [temp_red_spots isp];
                    end
                elseif ix == length(red_x_range) && iy == length(red_z_range)
                    if red_cell_center_x(isp) > red_x_range(ix) && red_cell_center_z(isp) > red_z_range(iy)
                        temp_red_spots = [temp_red_spots isp];
                    end
                elseif red_cell_center_x(isp) > red_x_range(ix) && red_cell_center_x(isp) <= red_x_range(ix+1) && red_cell_center_z(isp) > red_z_range(iy) && red_cell_center_z(isp) <= red_z_range(iy+1)
                    temp_red_spots = [temp_red_spots isp];
                end
            end
            binned_red_spots(iy,ix) = length(temp_red_spots);
        end
    end
    
    for ix = 1:length(red_x_range)
        for iy = 1:length(red_z_range)
            temp_red_spots_membrane = [];
            for isp = 1:length(red_cell_center_x_membrane)
                if ix == length(red_x_range) && iy < length(red_z_range)
                    if red_cell_center_x_membrane(isp) > red_x_range(ix) && red_cell_center_z_membrane(isp) > red_z_range(iy) && red_cell_center_z_membrane(isp) <= red_z_range(iy+1)
                        temp_red_spots_membrane = [temp_red_spots_membrane isp];
                    end
                elseif ix < length(red_x_range) && iy == length(red_z_range)
                    if red_cell_center_x_membrane(isp) > red_x_range(ix) && red_cell_center_x_membrane(isp) <= red_x_range(ix+1) && red_cell_center_z_membrane(isp) > red_z_range(iy) 
                        temp_red_spots_membrane = [temp_red_spots_membrane isp];
                    end
                elseif ix == length(red_x_range) && iy == length(red_z_range)
                    if red_cell_center_x_membrane(isp) > red_x_range(ix) && red_cell_center_z_membrane(isp) > red_z_range(iy)
                        temp_red_spots_membrane = [temp_red_spots_membrane isp];
                    end
                elseif red_cell_center_x_membrane(isp) > red_x_range(ix) && red_cell_center_x_membrane(isp) <= red_x_range(ix+1) && red_cell_center_z_membrane(isp) > red_z_range(iy) && red_cell_center_z_membrane(isp) <= red_z_range(iy+1)
                    temp_red_spots_membrane = [temp_red_spots_membrane isp];
                end
            end
            binned_red_spots_membrane(iy,ix) = length(temp_red_spots_membrane);
        end
    end
    
    for ix = 1:length(red_x_range)
        for iy = 1:length(red_z_range)
            temp_red_spots_pole = [];
            for isp = 1:length(red_cell_center_x_pole)
                if ix == length(red_x_range) && iy < length(red_z_range)
                    if red_cell_center_x_pole(isp) > red_x_range(ix) && red_cell_center_z_pole(isp) > red_z_range(iy) && red_cell_center_z_pole(isp) <= red_z_range(iy+1)
                        temp_red_spots_pole = [temp_red_spots_pole isp];
                    end
                elseif ix < length(red_x_range) && iy == length(red_z_range)
                    if red_cell_center_x_pole(isp) > red_x_range(ix) && red_cell_center_x_pole(isp) <= red_x_range(ix+1) && red_cell_center_z_pole(isp) > red_z_range(iy) 
                        temp_red_spots_pole = [temp_red_spots_pole isp];
                    end
                elseif ix == length(red_x_range) && iy == length(red_z_range)
                    if red_cell_center_x_pole(isp) > red_x_range(ix) && red_cell_center_z_pole(isp) > red_z_range(iy)
                        temp_red_spots_pole = [temp_red_spots_pole isp];
                    end
                elseif red_cell_center_x_pole(isp) > red_x_range(ix) && red_cell_center_x_pole(isp) <= red_x_range(ix+1) && red_cell_center_z_pole(isp) > red_z_range(iy) && red_cell_center_z_pole(isp) <= red_z_range(iy+1)
                    temp_red_spots_pole = [temp_red_spots_pole isp];
                end
            end
            binned_red_spots_pole(iy,ix) = length(temp_red_spots_pole);
        end
    end
    
    for ix = 1:length(red_x_range)
        for iy = 1:length(red_z_range)
            temp_red_spots_cyto = [];
            for isp = 1:length(red_cell_center_x_cyto)
                if ix == length(red_x_range) && iy < length(red_z_range)
                    if red_cell_center_x_cyto(isp) > red_x_range(ix) && red_cell_center_z_cyto(isp) > red_z_range(iy) && red_cell_center_z_cyto(isp) <= red_z_range(iy+1)
                        temp_red_spots_cyto = [temp_red_spots_cyto isp];
                    end
                elseif ix < length(red_x_range) && iy == length(red_z_range)
                    if red_cell_center_x_cyto(isp) > red_x_range(ix) && red_cell_center_x_cyto(isp) <= red_x_range(ix+1) && red_cell_center_z_cyto(isp) > red_z_range(iy) 
                        temp_red_spots_cyto = [temp_red_spots_cyto isp];
                    end
                elseif ix == length(red_x_range) && iy == length(red_z_range)
                    if red_cell_center_x_cyto(isp) > red_x_range(ix) && red_cell_center_z_cyto(isp) > red_z_range(iy)
                        temp_red_spots_cyto = [temp_red_spots_cyto isp];
                    end
                elseif red_cell_center_x_cyto(isp) > red_x_range(ix) && red_cell_center_x_cyto(isp) <= red_x_range(ix+1) && red_cell_center_z_cyto(isp) > red_z_range(iy) && red_cell_center_z_cyto(isp) <= red_z_range(iy+1)
                    temp_red_spots_cyto = [temp_red_spots_cyto isp];
                end
            end
            binned_red_spots_cyto(iy,ix) = length(temp_red_spots_cyto);
        end
    end
    
    for ix = 1:length(red_x_range)
        for iy = 1:length(red_z_range)
            temp_red_spots_dapi = [];
            for isp = 1:length(red_cell_center_x_dapi)
                if ix == length(red_x_range) && iy < length(red_z_range)
                    if red_cell_center_x_dapi(isp) > red_x_range(ix) && red_cell_center_z_dapi(isp) > red_z_range(iy) && red_cell_center_z_dapi(isp) <= red_z_range(iy+1)
                        temp_red_spots_dapi = [temp_red_spots_dapi isp];
                    end
                elseif ix < length(red_x_range) && iy == length(red_z_range)
                    if red_cell_center_x_dapi(isp) > red_x_range(ix) && red_cell_center_x_dapi(isp) <= red_x_range(ix+1) && red_cell_center_z_dapi(isp) > red_z_range(iy) 
                        temp_red_spots_dapi = [temp_red_spots_dapi isp];
                    end
                elseif ix == length(red_x_range) && iy == length(red_z_range)
                    if red_cell_center_x_dapi(isp) > red_x_range(ix) && red_cell_center_z_dapi(isp) > red_z_range(iy)
                        temp_red_spots_dapi = [temp_red_spots_dapi isp];
                    end
                elseif red_cell_center_x_dapi(isp) > red_x_range(ix) && red_cell_center_x_dapi(isp) <= red_x_range(ix+1) && red_cell_center_z_dapi(isp) > red_z_range(iy) && red_cell_center_z_dapi(isp) <= red_z_range(iy+1)
                    temp_red_spots_dapi = [temp_red_spots_dapi isp];
                end
            end
            binned_red_spots_dapi(iy,ix) = length(temp_red_spots_dapi);
        end
    end
    
    binned_green_spots = zeros(length(green_z_range),length(green_x_range));
    binned_green_spots_membrane = zeros(length(green_z_range),length(green_x_range));
    binned_green_spots_dapi = zeros(length(green_z_range),length(green_x_range));
    binned_green_spots_cyto = zeros(length(green_z_range),length(green_x_range));
    binned_green_spots_pole = zeros(length(green_z_range),length(green_x_range));
    
    for ix = 1:length(green_x_range)
        for iy = 1:length(green_z_range)
            temp_green_spots = [];
            for isp = 1:length(green_cell_center_x)
                if ix == length(green_x_range) && iy < length(green_z_range)
                    if green_cell_center_x(isp) > green_x_range(ix) && green_cell_center_z(isp) > green_z_range(iy) && green_cell_center_z(isp) <= green_z_range(iy+1)
                        temp_green_spots = [temp_green_spots isp];
                    end
                elseif ix < length(green_x_range) && iy == length(green_z_range)
                    if green_cell_center_x(isp) > green_x_range(ix) && green_cell_center_x(isp) <= green_x_range(ix+1) && green_cell_center_z(isp) > green_z_range(iy) 
                        temp_green_spots = [temp_green_spots isp];
                    end
                elseif ix == length(green_x_range) && iy == length(green_z_range)
                    if green_cell_center_x(isp) > green_x_range(ix) && green_cell_center_z(isp) > green_z_range(iy)
                        temp_green_spots = [temp_green_spots isp];
                    end
                elseif green_cell_center_x(isp) > green_x_range(ix) && green_cell_center_x(isp) <= green_x_range(ix+1) && green_cell_center_z(isp) > green_z_range(iy) && green_cell_center_z(isp) <= green_z_range(iy+1)
                    temp_green_spots = [temp_green_spots isp];
                end
            end
            binned_green_spots(iy,ix) = length(temp_green_spots);
        end
    end
    
    for ix = 1:length(green_x_range)
        for iy = 1:length(green_z_range)
            temp_green_spots_membrane = [];
            for isp = 1:length(green_cell_center_x_membrane)
                if ix == length(green_x_range) && iy < length(green_z_range)
                    if green_cell_center_x_membrane(isp) > green_x_range(ix) && green_cell_center_z_membrane(isp) > green_z_range(iy) && green_cell_center_z_membrane(isp) <= green_z_range(iy+1)
                        temp_green_spots_membrane = [temp_green_spots_membrane isp];
                    end
                elseif ix < length(green_x_range) && iy == length(green_z_range)
                    if green_cell_center_x_membrane(isp) > green_x_range(ix) && green_cell_center_x_membrane(isp) <= green_x_range(ix+1) && green_cell_center_z_membrane(isp) > green_z_range(iy) 
                        temp_green_spots_membrane = [temp_green_spots_membrane isp];
                    end
                elseif ix == length(green_x_range) && iy == length(green_z_range)
                    if green_cell_center_x_membrane(isp) > green_x_range(ix) && green_cell_center_z_membrane(isp) > green_z_range(iy)
                        temp_green_spots_membrane = [temp_green_spots_membrane isp];
                    end
                elseif green_cell_center_x_membrane(isp) > green_x_range(ix) && green_cell_center_x_membrane(isp) <= green_x_range(ix+1) && green_cell_center_z_membrane(isp) > green_z_range(iy) && green_cell_center_z_membrane(isp) <= green_z_range(iy+1)
                    temp_green_spots_membrane = [temp_green_spots_membrane isp];
                end
            end
            binned_green_spots_membrane(iy,ix) = length(temp_green_spots_membrane);
        end
    end
    
    for ix = 1:length(green_x_range)
        for iy = 1:length(green_z_range)
            temp_green_spots_pole = [];
            for isp = 1:length(green_cell_center_x_pole)
                if ix == length(green_x_range) && iy < length(green_z_range)
                    if green_cell_center_x_pole(isp) > green_x_range(ix) && green_cell_center_z_pole(isp) > green_z_range(iy) && green_cell_center_z_pole(isp) <= green_z_range(iy+1)
                        temp_green_spots_pole = [temp_green_spots_pole isp];
                    end
                elseif ix < length(green_x_range) && iy == length(green_z_range)
                    if green_cell_center_x_pole(isp) > green_x_range(ix) && green_cell_center_x_pole(isp) <= green_x_range(ix+1) && green_cell_center_z_pole(isp) > green_z_range(iy) 
                        temp_green_spots_pole = [temp_green_spots_pole isp];
                    end
                elseif ix == length(green_x_range) && iy == length(green_z_range)
                    if green_cell_center_x_pole(isp) > green_x_range(ix) && green_cell_center_z_pole(isp) > green_z_range(iy)
                        temp_green_spots_pole = [temp_green_spots_pole isp];
                    end
                elseif green_cell_center_x_pole(isp) > green_x_range(ix) && green_cell_center_x_pole(isp) <= green_x_range(ix+1) && green_cell_center_z_pole(isp) > green_z_range(iy) && green_cell_center_z_pole(isp) <= green_z_range(iy+1)
                    temp_green_spots_pole = [temp_green_spots_pole isp];
                end
            end
            binned_green_spots_pole(iy,ix) = length(temp_green_spots_pole);
        end
    end
    
    for ix = 1:length(green_x_range)
        for iy = 1:length(green_z_range)
            temp_green_spots_cyto = [];
            for isp = 1:length(green_cell_center_x_cyto)
                if ix == length(green_x_range) && iy < length(green_z_range)
                    if green_cell_center_x_cyto(isp) > green_x_range(ix) && green_cell_center_z_cyto(isp) > green_z_range(iy) && green_cell_center_z_cyto(isp) <= green_z_range(iy+1)
                        temp_green_spots_cyto = [temp_green_spots_cyto isp];
                    end
                elseif ix < length(green_x_range) && iy == length(green_z_range)
                    if green_cell_center_x_cyto(isp) > green_x_range(ix) && green_cell_center_x_cyto(isp) <= green_x_range(ix+1) && green_cell_center_z_cyto(isp) > green_z_range(iy) 
                        temp_green_spots_cyto = [temp_green_spots_cyto isp];
                    end
                elseif ix == length(green_x_range) && iy == length(green_z_range)
                    if green_cell_center_x_cyto(isp) > green_x_range(ix) && green_cell_center_z_cyto(isp) > green_z_range(iy)
                        temp_green_spots_cyto = [temp_green_spots_cyto isp];
                    end
                elseif green_cell_center_x_cyto(isp) > green_x_range(ix) && green_cell_center_x_cyto(isp) <= green_x_range(ix+1) && green_cell_center_z_cyto(isp) > green_z_range(iy) && green_cell_center_z_cyto(isp) <= green_z_range(iy+1)
                    temp_green_spots_cyto = [temp_green_spots_cyto isp];
                end
            end
            binned_green_spots_cyto(iy,ix) = length(temp_green_spots_cyto);
        end
    end
    
    for ix = 1:length(green_x_range)
        for iy = 1:length(green_z_range)
            temp_green_spots_dapi = [];
            for isp = 1:length(green_cell_center_x_dapi)
                if ix == length(green_x_range) && iy < length(green_z_range)
                    if green_cell_center_x_dapi(isp) > green_x_range(ix) && green_cell_center_z_dapi(isp) > green_z_range(iy) && green_cell_center_z_dapi(isp) <= green_z_range(iy+1)
                        temp_green_spots_dapi = [temp_green_spots_dapi isp];
                    end
                elseif ix < length(green_x_range) && iy == length(green_z_range)
                    if green_cell_center_x_dapi(isp) > green_x_range(ix) && green_cell_center_x_dapi(isp) <= green_x_range(ix+1) && green_cell_center_z_dapi(isp) > green_z_range(iy) 
                        temp_green_spots_dapi = [temp_green_spots_dapi isp];
                    end
                elseif ix == length(green_x_range) && iy == length(green_z_range)
                    if green_cell_center_x_dapi(isp) > green_x_range(ix) && green_cell_center_z_dapi(isp) > green_z_range(iy)
                        temp_green_spots_dapi = [temp_green_spots_dapi isp];
                    end
                elseif green_cell_center_x_dapi(isp) > green_x_range(ix) && green_cell_center_x_dapi(isp) <= green_x_range(ix+1) && green_cell_center_z_dapi(isp) > green_z_range(iy) && green_cell_center_z_dapi(isp) <= green_z_range(iy+1)
                    temp_green_spots_dapi = [temp_green_spots_dapi isp];
                end
            end
            binned_green_spots_dapi(iy,ix) = length(temp_green_spots_dapi);
        end
    end
    
    
    figure(10*bin-9)
    if median(binned_red_spots(:)) < 1
        heatmap(binned_red_spots)
    else
        heatmap(binned_red_spots)%,'','',[],'MaxColorValue',(1+spot_cutoff)*median(binned_red_spots(:)),'MinColorValue',(1-spot_cutoff)*median(binned_red_spots(:)))
    end
    graph_title = strcat(['Total Red Spots for Cells between Lengths ', num2str(micron_bins(bin)), ' and ', num2str(micron_bins(bin+1)),' micron, n= ', num2str(length(cells_in_bin)),' Cells']);
    title(graph_title)
    xlabel('Cell X Axis')
    ylabel('Cell Z Axis')
    colormap autumn
    colorbar
    %pbaspect([mean(widths_in_bin)/2 mean(lengths_in_bin)/2 mean(widths_in_bin)/2])
    set(gcf,'position',[835,883,868,667])
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([dataFile,'/cells',num2str(micron_bins(bin)), '_', num2str(micron_bins(bin+1)),'_TotalRedSpots']);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([dataFile,'/cells',num2str(micron_bins(bin)), '_', num2str(micron_bins(bin+1)),'_TotalRedSpots.fig']);
    savefig(gcf,file1_fig)

    figure(10*bin-8)
    if median(binned_red_spots_membrane(:)) < 1
        heatmap(binned_red_spots_membrane)
    else
        heatmap(binned_red_spots_membrane)%,'','',[],'MaxColorValue',(1+spot_cutoff)*median(binned_red_spots(:)),'MinColorValue',(1-spot_cutoff)*median(binned_red_spots_membrane(:)))
    end
    graph_title = strcat(['Red Membrane Spots for Cells between Lengths ', num2str(micron_bins(bin)), ' and ', num2str(micron_bins(bin+1)),' micron, n= ', num2str(length(cells_in_bin)),' Cells']);
    title(graph_title)
    xlabel('Cell X Axis')
    ylabel('Cell Z Axis')
    colormap autumn
    colorbar
    %pbaspect([mean(widths_in_bin)/2 mean(lengths_in_bin)/2 mean(widths_in_bin)/2])
    set(gcf,'position',[835,883,868,667])
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([dataFile,'/cells',num2str(micron_bins(bin)), '_', num2str(micron_bins(bin+1)),'_RedMembraneSpots']);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([dataFile,'/cells',num2str(micron_bins(bin)), '_', num2str(micron_bins(bin+1)),'_RedMembraneSpots.fig']);
    savefig(gcf,file1_fig)
    
    figure(10*bin-7)
    if median(binned_red_spots_dapi(:)) < 1
        heatmap(binned_red_spots_dapi)
    else
        heatmap(binned_red_spots_dapi)%,'','',[],'MaxColorValue',(1+spot_cutoff)*median(binned_red_spots(:)),'MinColorValue',(1-spot_cutoff)*median(binned_red_spots_dapi(:)))
    end
    graph_title = strcat(['Red Dapi Spots for Cells between Lengths ', num2str(micron_bins(bin)), ' and ', num2str(micron_bins(bin+1)),' micron, n= ', num2str(length(cells_in_bin)),' Cells']);
    title(graph_title)
    xlabel('Cell X Axis')
    ylabel('Cell Z Axis')
    colormap autumn
    colorbar
    %pbaspect([mean(widths_in_bin)/2 mean(lengths_in_bin)/2 mean(widths_in_bin)/2])
    set(gcf,'position',[835,883,868,667])
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([dataFile,'/cells',num2str(micron_bins(bin)), '_', num2str(micron_bins(bin+1)),'_RedDapiSpots']);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([dataFile,'/cells',num2str(micron_bins(bin)), '_', num2str(micron_bins(bin+1)),'_RedDapiSpots.fig']);
    savefig(gcf,file1_fig)
    
    figure(10*bin-6)
    if median(binned_red_spots_cyto(:)) < 1
        heatmap(binned_red_spots_cyto)
    else
        heatmap(binned_red_spots_cyto)%,'','',[],'MaxColorValue',(1+spot_cutoff)*median(binned_red_spots(:)),'MinColorValue',(1-spot_cutoff)*median(binned_red_spots_cyto(:)))
    end
    graph_title = strcat(['Red Cytoplasm Spots for Cells between Lengths ', num2str(micron_bins(bin)), ' and ', num2str(micron_bins(bin+1)),' micron, n= ', num2str(length(cells_in_bin)),' Cells']);
    title(graph_title)
    xlabel('Cell X Axis')
    ylabel('Cell Z Axis')
    colormap autumn
    colorbar
    %pbaspect([mean(widths_in_bin)/2 mean(lengths_in_bin)/2 mean(widths_in_bin)/2])
    set(gcf,'position',[835,883,868,667])
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([dataFile,'/cells',num2str(micron_bins(bin)), '_', num2str(micron_bins(bin+1)),'_RedCytoplasmSpots']);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([dataFile,'/cells',num2str(micron_bins(bin)), '_', num2str(micron_bins(bin+1)),'_RedCytoplasmSpots.fig']);
    savefig(gcf,file1_fig)
    
    figure(10*bin-5)
    if median(binned_red_spots_pole(:)) < 1
        heatmap(binned_red_spots_pole)
    else
        heatmap(binned_red_spots_pole)%,'','',[],'MaxColorValue',(1+spot_cutoff)*median(binned_red_spots(:)),'MinColorValue',(1-spot_cutoff)*median(binned_red_spots_pole(:)))
    end
    graph_title = strcat(['Red Pole Spots for Cells between Lengths ', num2str(micron_bins(bin)), ' and ', num2str(micron_bins(bin+1)),' micron, n= ', num2str(length(cells_in_bin)),' Cells']);
    title(graph_title)
    xlabel('Cell X Axis')
    ylabel('Cell Z Axis')
    colormap autumn
    colorbar
    %pbaspect([mean(widths_in_bin)/2 mean(lengths_in_bin)/2 mean(widths_in_bin)/2])
    set(gcf,'position',[835,883,868,667])
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([dataFile,'/cells',num2str(micron_bins(bin)), '_', num2str(micron_bins(bin+1)),'_RedPoleSpots']);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([dataFile,'/cells',num2str(micron_bins(bin)), '_', num2str(micron_bins(bin+1)),'_RedPoleSpots.fig']);
    savefig(gcf,file1_fig)
    
    figure(10*bin-4)
    if median(binned_green_spots(:)) < 1
        heatmap(binned_green_spots)
    else
        heatmap(binned_green_spots)%,'','',[],'MaxColorValue',(1+spot_cutoff)*median(binned_green_spots(:)),'MinColorValue',(1-spot_cutoff)*median(binned_green_spots(:)))
    end
    graph_title = strcat(['Total Green Spots for Cells between Lengths ', num2str(micron_bins(bin)), ' and ', num2str(micron_bins(bin+1)),' micron, n= ', num2str(length(cells_in_bin)),' Cells']);
    title(graph_title)
    xlabel('Cell X Axis')
    ylabel('Cell Z Axis')
    colormap winter
    colorbar
    %pbaspect([mean(widths_in_bin)/2 mean(lengths_in_bin)/2 mean(widths_in_bin)/2])
    set(gcf,'position',[835,883,868,667])
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([dataFile,'/cells',num2str(micron_bins(bin)), '_', num2str(micron_bins(bin+1)),'_TotalGreenSpots']);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([dataFile,'/cells',num2str(micron_bins(bin)), '_', num2str(micron_bins(bin+1)),'_TotalGreenSpots.fig']);
    savefig(gcf,file1_fig)

    figure(10*bin-3)
    if median(binned_green_spots_membrane(:)) < 1
        heatmap(binned_green_spots_membrane)
    else
        heatmap(binned_green_spots_membrane)%,'','',[],'MaxColorValue',(1+spot_cutoff)*median(binned_green_spots(:)),'MinColorValue',(1-spot_cutoff)*median(binned_green_spots_membrane(:)))
    end
    graph_title = strcat(['Green Membrane Spots for Cells between Lengths ', num2str(micron_bins(bin)), ' and ', num2str(micron_bins(bin+1)),' micron, n= ', num2str(length(cells_in_bin)),' Cells']);
    title(graph_title)
    xlabel('Cell X Axis')
    ylabel('Cell Z Axis')
    colormap winter
    colorbar
    %pbaspect([mean(widths_in_bin)/2 mean(lengths_in_bin)/2 mean(widths_in_bin)/2])
    set(gcf,'position',[835,883,868,667])
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([dataFile,'/cells',num2str(micron_bins(bin)), '_', num2str(micron_bins(bin+1)),'_GreenMembraneSpots']);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([dataFile,'/cells',num2str(micron_bins(bin)), '_', num2str(micron_bins(bin+1)),'_GreenMembraneSpots.fig']);
    savefig(gcf,file1_fig)
    
    figure(10*bin-2)
    if median(binned_green_spots_dapi(:)) < 1
        heatmap(binned_green_spots_dapi)
    else
        heatmap(binned_green_spots_dapi)%,'','',[],'MaxColorValue',(1+spot_cutoff)*median(binned_green_spots(:)),'MinColorValue',(1-spot_cutoff)*median(binned_green_spots_dapi(:)))
    end
    graph_title = strcat(['Green Dapi Spots for Cells between Lengths ', num2str(micron_bins(bin)), ' and ', num2str(micron_bins(bin+1)),' micron, n= ', num2str(length(cells_in_bin)),' Cells']);
    title(graph_title)
    xlabel('Cell X Axis')
    ylabel('Cell Z Axis')
    colormap winter
    colorbar
    %pbaspect([mean(widths_in_bin)/2 mean(lengths_in_bin)/2 mean(widths_in_bin)/2])
    set(gcf,'position',[835,883,868,667])
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([dataFile,'/cells',num2str(micron_bins(bin)), '_', num2str(micron_bins(bin+1)),'_GreenDapiSpots']);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([dataFile,'/cells',num2str(micron_bins(bin)), '_', num2str(micron_bins(bin+1)),'_GreenDapiSpots.fig']);
    savefig(gcf,file1_fig)
    
    figure(10*bin-1)
    if median(binned_green_spots_cyto(:)) < 1
        heatmap(binned_green_spots_cyto)
    else
        heatmap(binned_green_spots_cyto)%,'','',[],'MaxColorValue',(1+spot_cutoff)*median(binned_green_spots(:)),'MinColorValue',(1-spot_cutoff)*median(binned_green_spots_cyto(:)))
    end
    graph_title = strcat(['Green Cytoplasm Spots for Cells between Lengths ', num2str(micron_bins(bin)), ' and ', num2str(micron_bins(bin+1)),' micron, n= ', num2str(length(cells_in_bin)),' Cells']);
    title(graph_title)
    xlabel('Cell X Axis')
    ylabel('Cell Z Axis')
    colormap winter
    colorbar
    %pbaspect([mean(widths_in_bin)/2 mean(lengths_in_bin)/2 mean(widths_in_bin)/2])
    set(gcf,'position',[835,883,868,667])
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([dataFile,'/cells',num2str(micron_bins(bin)), '_', num2str(micron_bins(bin+1)),'_GreenCytoplasmSpots']);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([dataFile,'/cells',num2str(micron_bins(bin)), '_', num2str(micron_bins(bin+1)),'_GreenCytoplasmSpots.fig']);
    savefig(gcf,file1_fig)
    
    figure(10*bin)
    if median(binned_green_spots_pole(:)) < 1
        heatmap(binned_green_spots_pole)
    else
        heatmap(binned_green_spots_pole)%,'','',[],'MaxColorValue',(1+spot_cutoff)*median(binned_green_spots(:)),'MinColorValue',(1-spot_cutoff)*median(binned_green_spots_pole(:)))
    end
    graph_title = strcat(['Green Pole Spots for Cells between Lengths ', num2str(micron_bins(bin)), ' and ', num2str(micron_bins(bin+1)),' micron, n= ', num2str(length(cells_in_bin)),' Cells']);
    title(graph_title)
    xlabel('Cell X Axis')
    ylabel('Cell Z Axis')
    colormap winter
    colorbar
    %pbaspect([mean(widths_in_bin)/2 mean(lengths_in_bin)/2 mean(widths_in_bin)/2])
    set(gcf,'position',[835,883,868,667])
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([dataFile,'/cells',num2str(micron_bins(bin)), '_', num2str(micron_bins(bin+1)),'_GreenPoleSpots']);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([dataFile,'/cells',num2str(micron_bins(bin)), '_', num2str(micron_bins(bin+1)),'_GreenPoleSpots.fig']);
    savefig(gcf,file1_fig)
    
    
    
end

