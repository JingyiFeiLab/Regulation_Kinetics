[mask,area,ellipticity,centers] = track_mask(dic_file,pixelscaling);
dapi_objects = dapiMask2(stack_dapi,.5);
d2 = im2bw(dapi_objects,.01);
C = imfuse(m2,d2,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]); 
close all
imshow(C)
continuetranslation=1;
while continuetranslation==1
prompt={'X_Translation(- = Left, + = Right):','Y_Translation(- = Down, + = Up):','X_Stretch (0-1 : Compress, >1 = Stretch','Y_Stretch (0-1 : Compress, >1 = Stretch','Continue the translation(1 == Yes && 0== NO) ?'};   % A box will take in the values for the X/Ytranslation
mask_title='Dapi (green) Translation';                             % The title of the box
answer=inputdlg(prompt,mask_title);
Xtranslation = str2num(answer{1}); 
Ytranslation = str2num(answer{2});
YStretch = round(str2num(answer{4})*Y);
XStretch = round(str2num(answer{3})*X);
if isempty(answer{1}) || str2num(answer{1}) == 0
    Xtranslation = 0;
end
if isempty(answer{2}) || str2num(answer{2}) == 0
    Ytranslation = 0;
end
if isempty(answer{4}) || str2num(answer{4}) == 0
    YStretch = size(d2,1);
end
if isempty(answer{3}) || str2num(answer{3}) == 0
    XStretch = size(d2,2);
end
continuetranslation = str2num(answer{5});
if isempty(answer{5}) 
    continuetranslation = 1;
end
d2 = imtranslate(d2,[Xtranslation,-Ytranslation]);
d2 = imresize(d2,[YStretch XStretch]);
C = imfuse(m2,d2,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
imshow(C)
end

xmask = size(mask,1);
ymask = size(mask,2);

d2 = bwlabel(C(:,:,2));
xd2 = size(d2,1);
yd2 = size(d2,2);
if xd2 > xmask
    d2(xmask+1:xd2,:) = [];
end

if yd2 > xmask
    d2(:,ymask+1:yd2,:) = [];
end
   
d3 = d2;
num_dapi = max(d2(:));

for j = 1:num_dapi
    temp_dapi_mask = d2==j;
    dapi_ids = [];
    dapi_ids = unique(temp_dapi_mask.*mask);
    dapi_ids(dapi_ids==0) = [];
    
    if length(dapi_ids) > 1
        dapi_split = concave_split(temp_dapi_mask,1,.130,2.5);
        d2(d2 == j) = 0;
        d2 = d2 + dapi_split;
    end
        
end

d2 = bwlabel(imdilate(imerode(d2,se),se),4);
d3 = d2;

for j = 1:num_dapi
    temp_dapi_mask = d2==j;
    dapi_ids = [];
    dapi_ids = unique(temp_dapi_mask.*mask);
    dapi_ids(dapi_ids==0) = [];
    if isempty(dapi_ids)
        continue
    end
    
    if length(dapi_ids) > 1
        mask_temp_split = zeros(size(mask));
        for split_i = 1:length(dapi_ids)
            mask_temp_split = mask_temp_split + (mask==dapi_ids(split_i));
            mask(mask==dapi_ids(split_i)) = 0;
        end
        mask_temp_split = bwconvhull(mask_temp_split);
        mask = mask + dapi_ids(1)*mask_temp_split;
        dapi_ids = dapi_ids(1);
    end
   
    d3(d2 == j) = dapi_ids;
    
    
end

d3 = bwlabel(d3,4);
d4 = d3;
num_dapi = max(d3(:));
clear centers area ellipticity
[new_mask,area,ellipticity,centers] = post_track_mask(mask,d3,pixelscaling);

for j = 1:max(new_mask(:))
    temp_dapi_mask = new_mask ==j;
    dapi_ids = [];
    dapi_ids = unique(temp_dapi_mask.*d3);
    dapi_ids(dapi_ids==0) = [];
    if isempty(dapi_ids)
        %d4(d3 == j) = 0;
        continue
    end
    
    for ij = 1:length(dapi_ids)
        d4(d3 == dapi_ids(ij)) = j;
    end
    
%     if length(dapi_ids) > 1
%         mask_temp_split = zeros(size(mask));
%         for split_i = 1:length(dapi_ids)
%             mask_temp_split = mask_temp_split + (new_mask==dapi_ids(split_i));
%             new_mask(new_mask==dapi_ids(split_i)) = 0;
%         end
%         mask_temp_split = bwconvhull(mask_temp_split);
%         new_mask = new_mask + dapi_ids(1)*mask_temp_split;
%         dapi_ids = dapi_ids(1);
%     end
%     
%     d4(d3 == j) = dapi_ids;
    
end

m4 = im2bw(new_mask,.01);
d4 = m4.*d4;
d4l = im2bw(d4,.01);
C = imfuse(m4,d4l,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
close all

for i = 1:max(new_mask(:))
    strcat(['Working on Cell ' num2str(i)])
    
    if sum(new_mask(:)==i) == 0 || isempty(cell_struct(i).Center)
        continue
    end
    
    cell_struct(i).Dapi_Boundaries = bwboundaries(d4==i);
    cell_struct(i).Transformed_Dapi_Boundaries = cell_struct(i).Dapi_Boundaries;
    
    if ~isempty(cell_struct(i).Dapi_Boundaries) 
        
        for di = 1:length(cell_struct(i).Dapi_Boundaries)
            
            cell_struct(i).Transformed_Dapi_Boundaries{di,1}(:,1) = cell_struct(i).Transformed_Dapi_Boundaries{di,1}(:,1) - cell_struct(i).Center(2);
            cell_struct(i).Transformed_Dapi_Boundaries{di,1}(:,2) = cell_struct(i).Transformed_Dapi_Boundaries{di,1}(:,2) - cell_struct(i).Center(1);
            dapi_row_border = cell_struct(i).Transformed_Dapi_Boundaries{di,1}(:,1);
            dapi_col_border = cell_struct(i).Transformed_Dapi_Boundaries{di,1}(:,2);
            cell_struct(i).Transformed_Dapi_Boundaries{di,1}(:,1) = (dapi_col_border*sin(cell_struct(i).Cell_Angle)+dapi_row_border*cos(cell_struct(i).Cell_Angle));
            cell_struct(i).Transformed_Dapi_Boundaries{di,1}(:,2) = (dapi_col_border*cos(cell_struct(i).Cell_Angle)-dapi_row_border*sin(cell_struct(i).Cell_Angle));
            
        end
    end
    
end
    