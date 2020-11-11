%Run this with your dataset loaded in the workspace already

image_slice = 3; %Which slice of the image do you want to check
channel = 3; %Which color channel do you want to check
se = [1,1,1;1,1,1;1,1,1];

if channel == 1
    J = imadjust(mat2gray(stack_one(:,:,image_slice)));
    IR = J;IG = J;IB = J;
    IR(bwperim(imdilate(total_cells_one(:,:,image_slice),se))) = 0;IG(bwperim(imdilate(total_cells_one(:,:,3),se))) = 0;IB(bwperim(imdilate(total_cells_one(:,:,3),se))) = 0;
    IR(bwperim(imdilate(total_cells_one(:,:,3),se))) = 100;
    IRGB = cat(3,IR,IG,IB);
    figure(1);imshow(IRGB)
elseif channel == 2
    J = imadjust(mat2gray(stack_two(:,:,image_slice)));
    IR = J;IG = J;IB = J;
    IR(bwperim(imdilate(total_cells_two(:,:,image_slice),se))) = 0;IG(bwperim(imdilate(total_cells_two(:,:,3),se))) = 0;IB(bwperim(imdilate(total_cells_two(:,:,3),se))) = 0;
    IR(bwperim(imdilate(total_cells_two(:,:,3),se))) = 100;
    IRGB = cat(3,IR,IG,IB);
    figure(1);imshow(IRGB)
elseif channel == 3
    J = imadjust(mat2gray(stack_three(:,:,image_slice)));
    IR = J;IG = J;IB = J;
    IR(bwperim(imdilate(total_cells_three(:,:,image_slice),se))) = 0;IG(bwperim(imdilate(total_cells_three(:,:,3),se))) = 0;IB(bwperim(imdilate(total_cells_three(:,:,3),se))) = 0;
    IR(bwperim(imdilate(total_cells_three(:,:,3),se))) = 100;
    IRGB = cat(3,IR,IG,IB);
    figure(1);imshow(IRGB)
elseif channel == 4
    J = imadjust(mat2gray(stack_four(:,:,image_slice)));
    IR = J;IG = J;IB = J;
    IR(bwperim(imdilate(total_cells_four(:,:,image_slice),se))) = 0;IG(bwperim(imdilate(total_cells_four(:,:,3),se))) = 0;IB(bwperim(imdilate(total_cells_four(:,:,3),se))) = 0;
    IR(bwperim(imdilate(total_cells_four(:,:,3),se))) = 100;
    IRGB = cat(3,IR,IG,IB);
    figure(1);imshow(IRGB)
end
    
    
    