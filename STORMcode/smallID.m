function BW = smallID(BW,small_size)

BWi = bwlabel(BW,4);
num = max(BWi(:));

for i = 1:num
    if sum(BWi(:) == i) < small_size
        BW(BWi == i) = 0;
    end
end

BW;

end