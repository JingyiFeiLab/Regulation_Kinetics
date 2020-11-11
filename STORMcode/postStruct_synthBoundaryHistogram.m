close all

folderTitle = '/Users/reyer/Data/STORM/05042019_Syn_data/Drug_Synthetic/DIP/';
stress = 'DIP';

bound_x = [];
bound_y = [];
expanded_bound_x = [];
expanded_bound_y = [];
constricted_bound_x = [];
constricted_bound_y = [];
spot_x = [];
spot_y = [];
membrane_correction = .1*boundary;

for i = 1:length(cell_struct)
    if isempty(cell_struct(i).Spots)
        continue
    end
    
    for j = 1:length(cell_struct(i).Transformed_Boundaries{1,1})
        bound_x = [bound_x; cell_struct(i).Transformed_Boundaries{1,1}(j,2)];
        bound_y = [bound_y; cell_struct(i).Transformed_Boundaries{1,1}(j,1)];
    end
    
%     for k = 1:length(cell_struct(i).Transformed_Expanded_Boundaries{1,1})
%         expanded_bound_x = [expanded_bound_x; cell_struct(i).Transformed_Expanded_Boundaries{1,1}(k,2)];
%         expanded_bound_y = [expanded_bound_y; cell_struct(i).Transformed_Expanded_Boundaries{1,1}(k,1)];
%     end
%     
%     for l = 1:length(cell_struct(i).Transformed_Constricted_Boundaries{1,1})
%         constricted_bound_x = [constricted_bound_x; cell_struct(i).Transformed_Constricted_Boundaries{1,1}(l,2)];
%         constricted_bound_y = [constricted_bound_y; cell_struct(i).Transformed_Constricted_Boundaries{1,1}(l,1)];
%     end
%     
    for j = 1:length(cell_struct(i).Spots(:,2))
        if spot_struct(cell_struct(i).Spots(j,1)).Distance2Membrane < membrane_correction
            continue
        end
        
        if cell_struct(i).Spots(j,5) == 0
            continue
        end
        
        
        spot_x = [spot_x; spot_struct(cell_struct(i).Spots(j,1)).Collapsed_2D_Coordinate(1)];
        spot_y = [spot_y; spot_struct(cell_struct(i).Spots(j,1)).Collapsed_2D_Coordinate(2)];
        
%         spot_x = [spot_x; spot_struct(cell_struct(i).Spots(j,1)).Transform_3D_Coordinate(1)];
%         spot_y = [spot_y; spot_struct(cell_struct(i).Spots(j,1)).Transform_3D_Coordinate(2)];
    end
    
    
end


figure(1);histogram(spot_x);histogram(spot_x);title(strcat([stress, ' Post-Cutoff Spot Distribution']),'FontSize',24)
set(gcf,'position',[835,883,868,667])
file1 = strcat([folderTitle,'postCut_spot_distribution']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'postCut_spot_distribution.fig']);
savefig(gcf,file1_fig)




