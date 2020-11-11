% Creates 2 output structures. 1) One entry for each cell containing its
% coordinates, center, boundaries, angle, and axis, and 2) One entry
% for each spot containing the cell ID it belongs to its 3D distance to
% center, and its distance to the membrane.

% This code will potentially take a while to run, depending on how many
% spots there are

% You should only have to change "pixelscaling", "dataset", "dic_file",
% "proper_set", and potentially the x and y columns for the STORM output
% (Seongjin's data is not consistent, so I have to check manually.
clear all
close all
graph_title = '100 GFP V A647 Correlation ';
parentDir = '/Users/reyer/Data/STORM/1xPBS_1ugPmL_lysozyme_no_DNAseI_2019_06_21/';

samples = [0,1,2];

for s = 1:length(samples)
    
    sampleDir = strcat([parentDir,'sample',num2str(samples(s))]);
    
    corrDir = strcat([sampleDir,'/gfp_correlations']);
    if exist(corrDir,'dir')~=7
        mkdir(sampleDir,'/gfp_correlations')
    end
    
    pixelscaling = 130; % Sometimes Seongjin's data is in nm, sometimes in um. If nm, set to 130. If um, set to .130
    % set "dataset" to STORM output .txt file
    
    
    dic_file = strcat([sampleDir,'/DIC_',num2str(samples(s)),'.tif']);
    gfp_file = strcat([sampleDir,'/GFP']);
    a647_file = strcat([sampleDir,'/A647']);
    stack_gfp = imFormat_wholeImage(gfp_file,1);
    stack_a647 = imFormat_wholeImage(a647_file,1);
    
    gfp_objects = dapiMask2(stack_gfp,.2);
    a647_objects = dapiMask2(stack_a647,.2);
    se = [1 1 1 ; 1 1 1 ; 1 1 1];
    
    
    %*********Scaling the X Y Coordinates as per the pixel of the input image
    filepath_dic = dic_file;
    dic_image = imread(filepath_dic);                  % Put the reference image based on which you want to select the points
    X = size(dic_image,1) ;                              % CapitalX is the number of X pixels in the reference image
    Y = size(dic_image,2) ;
    
    [mask, area,ellipticity,centers] = track_mask_SP(dic_file,pixelscaling);
    
    
%     m2 = im2bw(mask,.01);
%     d2 = im2bw(gfp_objects,.01);
%     C = imfuse(m2,d2,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
%     close all
%     C_lg =imresize(C,2);
%     imshow(C_lg)
%     continuetranslation=1;
%     while continuetranslation==1
%         prompt={'X_Translation(- = Left, + = Right):','Y_Translation(- = Down, + = Up):','X_Stretch (0-1 : Compress, >1 = Stretch','Y_Stretch (0-1 : Compress, >1 = Stretch','Continue the translation(1 == Yes && 0== NO) ?'};   % A box will take in the values for the X/Ytranslation
%         mask_title='GFP (green) Translation';                             % The title of the box
%         answer=inputdlg(prompt,mask_title);
%         Xtranslation = str2num(answer{1});
%         Ytranslation = str2num(answer{2});
%         YStretch = round(str2num(answer{4})*Y);
%         XStretch = round(str2num(answer{3})*X);
%         if isempty(answer{1}) || str2num(answer{1}) == 0
%             Xtranslation = 0;
%         end
%         if isempty(answer{2}) || str2num(answer{2}) == 0
%             Ytranslation = 0;
%         end
%         if isempty(answer{4}) || str2num(answer{4}) == 0
%             YStretch = size(d2,1);
%         end
%         if isempty(answer{3}) || str2num(answer{3}) == 0
%             XStretch = size(d2,2);
%         end
%         continuetranslation = str2num(answer{5});
%         if isempty(answer{5})
%             continuetranslation = 1;
%         end
%         d2 = imtranslate(d2,[Xtranslation,-Ytranslation]);
%         d2 = imresize(d2,[YStretch XStretch]);
%         stack_gfp = imtranslate(stack_gfp,[Xtranslation,-Ytranslation]);
%         stack_gfp = imresize(stack_gfp,[YStretch XStretch]);
%         
%         C = imfuse(m2,d2,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
%         
%         C_lg=imresize(C,2);
%         imshow(C_lg) %imshow(C)
%     end
%     
%     d3 = im2bw(a647_objects,.01);
%     C = imfuse(m2,d3,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
%     close all
%     C_lg =imresize(C,2);
%     imshow(C_lg)
%     continuetranslation=1;
%     while continuetranslation==1
%         prompt={'X_Translation(- = Left, + = Right):','Y_Translation(- = Down, + = Up):','X_Stretch (0-1 : Compress, >1 = Stretch','Y_Stretch (0-1 : Compress, >1 = Stretch','Continue the translation(1 == Yes && 0== NO) ?'};   % A box will take in the values for the X/Ytranslation
%         mask_title='Dapi (green) Translation';                             % The title of the box
%         answer=inputdlg(prompt,mask_title);
%         Xtranslation = str2num(answer{1});
%         Ytranslation = str2num(answer{2});
%         YStretch = round(str2num(answer{4})*Y);
%         XStretch = round(str2num(answer{3})*X);
%         if isempty(answer{1}) || str2num(answer{1}) == 0
%             Xtranslation = 0;
%         end
%         if isempty(answer{2}) || str2num(answer{2}) == 0
%             Ytranslation = 0;
%         end
%         if isempty(answer{4}) || str2num(answer{4}) == 0
%             YStretch = size(d3,1);
%         end
%         if isempty(answer{3}) || str2num(answer{3}) == 0
%             XStretch = size(d3,2);
%         end
%         continuetranslation = str2num(answer{5});
%         if isempty(answer{5})
%             continuetranslation = 1;
%         end
%         d3 = imtranslate(d3,[Xtranslation,-Ytranslation]);
%         d3 = imresize(d3,[YStretch XStretch]);
%         stack_a647 = imtranslate(stack_a647,[Xtranslation,-Ytranslation]);
%         stack_a647 = imresize(stack_a647,[YStretch XStretch]);
%         
%         C = imfuse(m2,d3,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
%         
%         C_lg=imresize(C,2);
%         imshow(C_lg) %imshow(C)
%     end
    
    
    field1 = 'Cell'; % All Objects, single and multi, labeled
    field2 = 'End_Points';
    field3 = 'Gfp_Intensities';
    field4 = 'a647_Intensities';
    field5 = 'Area';
    field6 = 'Cell_Angle';
    field7 = 'Cell_Y_Axis';
    field8 = 'Cell_X_Axis';
    cell_struct = struct(field1,[],field2, [], field3, [], field4, [],field5,[],field6, [], field7, [], field8, []);
    
    se = [1 1 1 ; 1 1 1 ; 1 1 1];
    counted = [];
    
    all_cells_good = 0;
    
    gv_over_mean = [];
    rv_over_mean = [];
    corr_of_coefs = [];
    cell_nums = [];
    cell_nums2 = [];
    
    for i = 1:max(mask(:))
        strcat(['Working on Cell ' num2str(i)])
        cell_struct(i).Cell = i;
        cell_struct(i).Center = [centers(i,2),centers(i,1)];
        cell_struct(i).Area = sum(mask(:)==i);
        cell_struct(i).Cell_Angle = -ellipticity(i,2)*(pi/180);
        cell_struct(i).Cell_Y_Axis = ellipticity(i,3);
        cell_struct(i).Cell_X_Axis = ellipticity(i,4);
        cell_struct(i).Gfp_Intensities = [];
        cell_struct(i).a647_Intensities = [];
        
        mask_i = imerode(mask==i,se);
        [cell_row, cell_col] = ind2sub(size(mask_i),find(mask_i));
        
        line = makeLine([cell_col(1),cell_col(length(cell_col))],[cell_row(1),cell_row(length(cell_row))]);
        if length(line(:,1)) < 3 
            continue
        end
        
        cell_struct(i).End_Points = [line(1,2),line(1,1);line(length(line),2),line(length(line),1)];
        
        
        
        for j = 3:length(line)-2
            cell_struct(i).Gfp_Intensities = [cell_struct(i).Gfp_Intensities; stack_gfp(line(j,1), line(j,2))];
            cell_struct(i).a647_Intensities = [cell_struct(i).a647_Intensities; stack_a647(line(j,1), line(j,2))];
        end
        
        if length(cell_struct(i).Gfp_Intensities) < 3 || length(cell_struct(i).a647_Intensities) < 3
            continue
        end
        
        cell_nums2 = [cell_nums2 i];
        pixels = 3:length(line(:,1))-2;
        gfp_i = cell_struct(i).Gfp_Intensities;
        a647_i = cell_struct(i).a647_Intensities;
        
        dline = polyfit(pixels,gfp_i',1);
        aline = polyfit(pixels,a647_i',1);
        
        d_data_line = dline(1)*pixels;
        d_data_line = d_data_line + (-1*gfp_i(1)-d_data_line(1));
        d_data_line = d_data_line - d_data_line(1);
        a_data_line = aline(1)*pixels;
        a_data_line = a_data_line + (-1*a647_i(1)-a_data_line(1));
        a_data_line = a_data_line - a_data_line(1);
        
        gfp_i_corrected = gfp_i' - d_data_line;
        a647_i_corrected = a647_i' - a_data_line;
        
        gfp_i = gfp_i_corrected;
        a647_i = a647_i_corrected;
        
        cell_struct(i).Gfp_Intensities = gfp_i;
        cell_struct(i).a647_Intensities = a647_i;
        
        
        
        
        gv_over_mean = [gv_over_mean std(gfp_i)/mean(gfp_i)];
        rv_over_mean = [rv_over_mean std(a647_i)/mean(a647_i)];
        [r,p] = corrcoef(gfp_i,a647_i);
        if p(1,2) < .9
            corr_of_coefs = [corr_of_coefs r(1,2)];
            cell_nums = [cell_nums i];
        end
        
        
        
        
        
        figure(1)
        ylabels{1}='Normalized GFP Intensity';
        ylabels{2}='Normalized A647 Intensity';
        [ax,hlines(1),hlines(2)] = plotyy(pixels,gfp_i/mean(gfp_i),pixels,a647_i/mean(a647_i));
        cfig = get(gcf,'color');
        pos = [0.1  0.1  0.7  0.8];
        offset = pos(3)/5.5;
        hlines(1).LineWidth = 6;
        hlines(2).LineWidth = 6;
        hlines(1).Color = 'b';
        hlines(2).Color = 'r';
        set(get(ax(1),'ylabel'),'string',ylabels{1})
        set(get(ax(2),'ylabel'),'string',ylabels{2})
        set(ax,{'ycolor'},{'b';'r'})
        
        
        
        
        title(strcat(['Cell ', num2str(i), ' Correlation Coefficient = ', num2str(r(1,2), '%.3f'),' P-value = ', num2str(p(1,2), '%.3f')]),'FontSize',26)
        xlabel('Pixel along Long Axis','FontSize',24)
        ax(1).FontSize = 18;
        ax(2).FontSize = 18;
        file1 = strcat([corrDir,'/gfp_a647_',num2str(i)]);
        set(gcf,'position',[835,883,868,667])
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-painters','-depsc','-r0')
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-dpng','-r0')
        file1_fig = strcat([corrDir,'/gfp_a647_',num2str(i),'.fig']);
        savefig(gcf,file1_fig)
        close all
        
        
    end
    
    cell_nums(isnan(corr_of_coefs)) = [];
    
    corr_of_coefs(isnan(corr_of_coefs)) = [];
    
    figure(1)
    hold on
    plot(cell_nums,zeros(length(cell_nums),1),'-k','LineWidth',3);
    hold on
    scatter(cell_nums,corr_of_coefs,100,'g','filled');
    corr2 = [corr2 corr_of_coefs];
    
    %title(strcat(['Correlation Coefficients, Mean CorrCoef = ', num2str(mean(corr_of_coefs), '%.4f')]),'FontSize',28)
    title(strcat(['Correlation Coefficients, Mean CorrCoef = ', num2str(mean(corr2), '%.4f')]),'FontSize',28)
    
    xlabel('Cell ID','FontSize',24)
    ylabel('Correlation Coefficient','FontSize',24)
    %file1 = strcat([corrDir,'/corrcoef']);
    set(gcf,'position',[835,883,868,667])
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([corrDir,'/corrcoef.fig']);
    savefig(gcf,file1_fig)
    
    
    figure(2)
    ylabels{1}='GFP C of Variation';
    ylabels{2}='A647 C of Variation';
    [ax,hlines(1),hlines(2)] = plotyy(cell_nums2,gv_over_mean,cell_nums2,rv_over_mean);
    cfig = get(gcf,'color');
    pos = [0.1  0.1  0.7  0.8];
    offset = pos(3)/5.5;
    hlines(1).LineWidth = 3;
    hlines(2).LineWidth = 3;
    hlines(1).Color = 'b';
    hlines(2).Color = 'r';
    set(get(ax(1),'ylabel'),'string',ylabels{1})
    set(get(ax(2),'ylabel'),'string',ylabels{2})
    set(ax,{'ycolor'},{'b';'r'})
    
    title(strcat(['Coefficients of Variation, GFP CoV = ', num2str(mean(gv_over_mean), '%.4f'),' A647 CoV = ', num2str(mean(rv_over_mean), '%.4f')]),'FontSize',26)
    xlabel('Cell ID','FontSize',24)
    ax(1).FontSize = 18;
    ax(2).FontSize = 18;
    file1 = strcat([corrDir,'/cov']);
    set(gcf,'position',[835,883,868,667])
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([corrDir,'/cov.fig']);
    savefig(gcf,file1_fig)
    close all
    
    saveFile = strcat([sampleDir,'/cellStructure_GFP.mat']);
    save(saveFile)
    
end



