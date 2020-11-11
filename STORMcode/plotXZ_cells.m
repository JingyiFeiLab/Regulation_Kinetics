clear cell_center_x cell_center_z cell_color cell_region

one_cell_id = 10; % If you want to create a scatter plot for just one cell, put it's ID here
plot_all_cells = 0; % If you want to create a heatmap for all spots in all cells, set to 1
percent_from_pole = .1; % Cutoff spots within this percentage of cell axis of cell pole

cell_center_x = [];
cell_center_z = [];
cell_color = [];
cell_region = [];

from_pole_temp = [];

levels = 10; %How many heights for cell spot density map? More takes longer
             % 10 takes a little bit of time, 100 a lot

if plot_all_cells == 0
    for si = 1:length(cell_struct(one_cell_id).Spots(:,1))
        from_pole = abs(cell_struct(one_cell_id).Spots(si,3));
        if (.5*cell_struct(one_cell_id).Cell_Y_Axis-from_pole)/(.5*cell_struct(one_cell_id).Cell_Y_Axis) > percent_from_pole
            
            cell_center_x = [cell_center_x cell_struct(one_cell_id).Spots(si,2)];
            cell_center_z = [cell_center_z cell_struct(one_cell_id).Spots(si,4)];
            cell_region = [cell_region cell_struct(one_cell_id).Spots(si,5)];
            if cell_struct(one_cell_id).Spots(si,5) == 4
                cell_color = [cell_color; [1 0 1]];
            elseif cell_struct(one_cell_id).Spots(si,5) == 3
                cell_color = [cell_color; [0 0 1]];
            elseif cell_struct(one_cell_id).Spots(si,5) == 1
                cell_color = [cell_color; [0 0 0]];
            elseif cell_struct(one_cell_id).Spots(si,5) == 2
                cell_color = [cell_color; [0 1 0]];
            end
            
            [cell_region,region_sort] = sort(cell_region);
            cell_center_x = cell_center_x(region_sort);
            cell_center_z = cell_center_z(region_sort);
            
               
        end
    end
    
    if sum(cell_region==3) ~= 0
        figure(3);gscatter(cell_center_x,cell_center_z,cell_region,'kgbm','.',10,'filled');grid on
    else
        figure(3);gscatter(cell_center_x,cell_center_z,cell_region,'kgm','.',10,'filled');grid on
    end
    title(strcat(['Cell ' num2str(one_cell_id) ' XZ Spot Projection']),'FontSize',24)
    ylabel('Vertical Axis (Z)','FontSize',24)
    xlabel('Short Axis (X)','FontSize',24)

end


if plot_all_cells == 1
    for ci = 1:length(cell_struct)
        
        if size(cell_struct(ci).Spots) == [0,0];
            continue
        end
        
        for si = 1:length(cell_struct(ci).Spots(:,1))
            from_pole = abs(cell_struct(ci).Spots(si,3));
            from_pole_temp = [from_pole_temp from_pole];
            if (.5*cell_struct(ci).Cell_Y_Axis-from_pole)/cell_struct(ci).Cell_Y_Axis > percent_from_pole
                
                cell_center_x = [cell_center_x cell_struct(ci).Spots(si,2)];
                cell_center_z = [cell_center_z cell_struct(ci).Spots(si,4)];
                
            end
        end
    end
    
    figure(3);
    map = dataDensity(cell_center_x, cell_center_z, 256, 256);
    map = map - min(min(map));
    map = floor(map ./ max(max(map)) * (levels-1));
    
    image(map);
    colormap(jet(levels));
    set(gca, 'XTick', [1 256]);
    set(gca, 'XTickLabel', [min(cell_center_x) max(cell_center_x)]);
    set(gca, 'YTick', [1 256]);
    set(gca, 'YTickLabel', [min(cell_center_z) max(cell_center_z)]);
    
    colorbar
    grid on
    title(strcat(['All Cells XZ Spot Projection']),'FontSize',24)
    ylabel('Vertical Axis (Z)','FontSize',24)
    xlabel('Short Axis (X)','FontSize',24)

end

cell_coordinates = [cell_center_x ; cell_center_z];

