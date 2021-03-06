%Output structure has one entry for each mRNA. Each entry contains the
%frames, location, and corresponding GFP intensity of that mRNA. After
%running this, you can run "suntag_scatter.m" without changing anything to
%generate images of mRNA spots over GFP images. You can also run
%"sunTag_Intensity.m" without any changes to generate the intensity graphs
%for each mRNA

%Set filepath to folder containing .tif stack
%cell_image_filepath = '/Users/reyer/Data/STORM/amine52919/sample/c5_sptPALM_m3.tif';


saveFile = '/Users/reyer/Data/STORM/heatmap72319/NT_Tr5'; %Where you want the images to save
treatment = 'NT';
% Direct to .txt file with mRNA tracking results


% These might change. Which columns represent x, y, z coordinates
x_col = 1;
y_col = 2;
z_col = 3;
spot_cutoff = .85;

% Don't Change any of these
channels = 1;
slices = 1;
pixel_scaling = .130;
window = 5; % Size of mRNA tracking window for GFP projection



%particles = textread(msd_file);

%Assumes diffusion coordinate is in column 4, ID is in Column 5. If not,
%change these

cell_center_x = [];
cell_center_y = [];
cell_diffusion = [];
cell_color = [];
cell_region = [];

for n = 1:length(cell_struct)
    
    if isempty(cell_struct(n).Spots)
        continue
    end
    
        
    for si = 1:length(cell_struct(n).Spots(:,1))
        
        cell_center_x = [cell_center_x cell_struct(n).Spots(si,2)];
        cell_center_y = [cell_center_y cell_struct(n).Spots(si,3)];
        cell_region = [cell_region cell_struct(n).Spots(si,5)];
        cell_diffusion = [cell_diffusion cell_struct(n).Spots(si,6)];
        
        
        [cell_region,region_sort] = sort(cell_region);
        cell_center_x = cell_center_x(region_sort);
        cell_center_y = cell_center_y(region_sort);
        
        
    end
    
    
    x_range = int32(min(cell_center_x):max(cell_center_x));
    y_range = int32(min(cell_center_y):max(cell_center_y));
    
    binned_diffusion = zeros(length(y_range),length(x_range));
    binned_spots = zeros(length(y_range),length(x_range));
    for ix = 1:length(x_range)
        for iy = 1:length(y_range)
            temp_diffusion = [];
            for isp = 1:length(cell_center_x)
                if ix == length(x_range) && iy < length(y_range)
                    if cell_center_x(isp) > x_range(ix) && cell_center_y(isp) > y_range(iy) && cell_center_y(isp) <= y_range(iy+1)
                        temp_diffusion = [temp_diffusion cell_diffusion(isp)];
                    end
                elseif ix < length(x_range) && iy == length(y_range)
                    if cell_center_x(isp) > x_range(ix) && cell_center_x(isp) <= x_range(ix+1) && cell_center_y(isp) > y_range(iy) 
                        temp_diffusion = [temp_diffusion cell_diffusion(isp)];
                    end
                elseif ix == length(x_range) && iy == length(y_range)
                    if cell_center_x(isp) > x_range(ix) && cell_center_y(isp) > y_range(iy)
                        temp_diffusion = [temp_diffusion cell_diffusion(isp)];
                    end
                elseif cell_center_x(isp) > x_range(ix) && cell_center_x(isp) <= x_range(ix+1) && cell_center_y(isp) > y_range(iy) && cell_center_y(isp) <= y_range(iy+1)
                    temp_diffusion = [temp_diffusion cell_diffusion(isp)];
                end
            end
            if isempty(temp_diffusion)
                binned_diffusion(iy,ix) = 0;
            else
                binned_diffusion(iy,ix) = mean(temp_diffusion);
                binned_spots(iy,ix) = length(temp_diffusion);
            end
        end
    end
    
    se = [1 1 1 ; 1 1 1 ; 1 1 1];
    new_spots = im2bw(imdilate(imdilate(imerode(binned_spots,se),se),se),.01);
    new_diffusion = new_spots.*binned_diffusion;
    
    figure(1)
    s = pcolor(new_diffusion);
    colormap([1 1 1; jet])
    shading flat
    axis square tight
    s.FaceColor = 'interp';
    caxis([0,.03])
    set(gca,'ydir','reverse')
    graph_title = strcat([treatment,', Diffusion Coefficients, Cell ', num2str(n)]);
    title(graph_title)
    xlabel('Cell X Axis')
    ylabel('Cell Y Axis')
    colorbar
    pbaspect([max(cell_center_x)/2 max(cell_center_y)/2 max(cell_center_x)/2])
    set(gcf,'position',[835,883,868,667])
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([saveFile,'/smoothDiffusion_cell',num2str(n)]);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([saveFile,'/smoothDiffusion_cell',num2str(n),'.fig']);
    savefig(gcf,file1_fig)
    
    figure(2)
    if median(binned_spots(:)) < 1
        heatmap(binned_spots,x_range,y_range)
    else
        heatmap(binned_spots,'','',[],'MaxColorValue',(1+spot_cutoff)*median(binned_spots(:)),'MinColorValue',(1-spot_cutoff)*median(binned_spots(:)))
    end
    graph_title = strcat([treatment,', Number of Spots for Cell ', num2str(n)]);
    title(graph_title)
    xlabel('Cell X Axis')
    ylabel('Cell Y Axis')
    colorbar
    pbaspect([max(cell_center_x)/2 max(cell_center_y)/2 max(cell_center_x)/2])
    set(gcf,'position',[835,883,868,667])
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([saveFile,'/cellSpots_cell',num2str(n)]);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([saveFile,'/cellSpots_cell',num2str(n),'.fig']);
    savefig(gcf,file1_fig)
    
    close all
end                     % mRNA Spot ID











