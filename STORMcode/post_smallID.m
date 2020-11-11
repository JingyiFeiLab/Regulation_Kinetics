function BW = post_smallID(BW,small_size)

num = max(BW(:));

for i = 1:num
    if sum(BW(:) == i) < small_size
        BW(BW == i) = 0;
    end
end

BW;

end