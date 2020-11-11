function [cell_mask,area,ellipticity,center] = track_mask2(dic_path,pixelscaling)
%file_array = numbers of samples to analyze
%dic_path = path to dic images

pix_size = pixelscaling; %Microns

int_thresh = .00001; % Intensity Threshold
shape2D_thresh = 4;  % <--- Splitting threshold, you can change this

background_thresh = .2;
pix_neigh = 5; %floor((.08/pix_size)*12); <- If you have no idea, try this

% Path to main file (i.e. channel) that you will use for segmentation
%filepath_dic = strcat(['/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/August_24_17_convert/+SgrS/t20/dic',num2str(cell_num),'.tif']);
filepath_dic = dic_path;
stack_o = mat2gray(imread(filepath_dic));



se = [1 1 1; 1 1 1 ; 1 1 1]; % Structuring Element for basic Erosion and dilation


field2 = 'Objects'; % All Objects, single and multi, labeled
field3 = 'Area';
field4 = 'Centers';
field5 = 'Ellipticity'; % Single Cell Selections
part1 = struct(field2, [], field3, [], field4, [], field5, []);  %Table for Part 1


stack2 = zeros(size(stack_o));
xdim = size(stack_o,1);
ydim = size(stack_o,2);

edge_cut = 10;



stack2 = anisodiff2D(stack_o,15,1/7,30,1);
I=stack_o;
[r,c] = size(I);
I2 = stack2;

I_low_pass = low_pass(I2,.01); %.025

stack2 = stack2 - I_low_pass;
stack2 = stack2./max(max(stack2));
a = bwareaopen(imdilate(imerode(bradley(stack2,[pix_neigh,pix_neigh],int_thresh),se),se),20);
%a = a - bwareaopen(a,400);

b = bwlabel(a,4);
a_temp = a;

for i = 1:max(max(b))
    if  sum(sum(((b==i).*stack2)))/cellArea(b,i) < background_thresh
        a_temp(b==i) = 0;
        b(b == i) = 0;
    end
end

a = a_temp;
a(1:edge_cut,:) = 0; a(r-edge_cut+1:r,:) = 0; a(:,1:edge_cut) = 0; a(:,c-edge_cut+1:c) = 0;
BW = bwareaopen(a,10);
[~,BW] = edgeBreak(BW);
BW = imfill(imclearborder(smallID(bwlabel(BW,4))),'holes');



objects = bwlabel(BW,4);

num = max(objects(:));

ellipse_error = zeros(num,1);
test_ellipse = {};

for i = 1:num
    
    [ellipse1,test1] = ellipseError(objects,i);
    if isempty(ellipse1) == 1 || isempty(test1) == 1
        ellipse_error(i) = shape2D_thresh+1;
        continue
    else
        i
        test_ellipse(i) = test1;
        ellipse_error(i) = ellipseTest(ellipse1,test1,cellArea(objects,i,pix_size),pix_size);
        
    end
end

for i = 1:num
    
    if ellipse_error(i) < shape2D_thresh
        objects(objects==i) = 0;
        object_temp = zeros(size(objects));
        for l = 1:length(test_ellipse{i})
            object_temp(test_ellipse{i}(l,1),test_ellipse{i}(l,2)) = i;
        end
        object_temp = imfill(object_temp);
        objects(object_temp == i) = i;
        objects = smallID(imfill(objects));
        
    end
end

uni_obs = unique(objects);
for numb = 1:length(unique(objects))-1
    objects(objects == uni_obs(numb+1)) = numb;
end
num = max(objects(:));

clear centers area ellipticity

area=zeros(num,1);

for i=1:num
    area(i) = cellArea(objects,i,pix_size);
end

centers=zeros(max(objects(:)),2);

for i=1:num
    centers(i,:) = cellCenter(objects,i);
end

ellipticity = zeros(num,4);

% Re-done Ellipticity Calculation
for i=1:num
    ellipticity(i,:) = cellEllipseSPT(objects,i);
end

%
% Structured Array with our Data
part1.Objects = objects;
part1.Area = area;
part1.Centers = centers;
part1.Ellipticity = ellipticity;

cell_mask = part1.Objects;
area = part1.Area;
ellipticity = part1.Ellipticity;
center = part1.Centers;

end