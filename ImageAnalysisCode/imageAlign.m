channel = 3;
slice = 3;
num_slices = 5;
new_mask = total_cells;


if channel == 1
    new_stack = stack_one(:,:,slice);
elseif channel == 2
    new_stack = stack_two(:,:,slice);
elseif channel == 3
    new_stack = stack_three(:,:,slice);    
elseif channel == 4
    new_stack = stack_four(:,:,slice);
end

close all
figure(1)
imshow(imadjust(mat2gray(new_stack))+bwperim(imdilate(new_mask,se)))
set(figure(1),'Position',[431 823 642 574])
continuetranslation=1;
while continuetranslation==1
prompt={'X_Translation(- = Down, + = Up):','Y_Translation(- = right, + = left):','Continue the translation(1 == Yes && 0== NO) ?'};   % A box will take in the values for the X/Ytranslation
title='Translation';                             % The title of the box
answer=inputdlg(prompt,title);
Xtranslation = str2num(answer{1}); 
Ytranslation = str2num(answer{2});
continuetranslation = str2num(answer{3});
new_mask = imtranslate(new_mask,[-1*Ytranslation,-1*Xtranslation]);

imshow(imadjust(mat2gray(new_stack))+bwperim(imdilate(new_mask,se)))
set(figure(1),'Position',[431 823 642 574])
end

for g = 1:num_slices
    new_total_cells(:,:,g) = new_mask;
end
    
if channel == 1
    for i = 1:index
        part4(i).Intensity_One = cellIntensity(new_total_cells,stack_one,i);
    end
elseif channel == 2
    for i = 1:index
        part4(i).Intensity_Two = cellIntensity(new_total_cells,stack_two,i);
    end
elseif channel == 3
    for i = 1:index
        part4(i).Intensity_Three = cellIntensity(new_total_cells,stack_three,i);
    end
elseif channel == 4
    for i = 1:index
        part4(i).Intensity_Four = cellIntensity(new_total_cells,stack_four,i);
    end
end

close all
save(filename)