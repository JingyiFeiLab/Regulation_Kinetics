function RNA = calc_copy_num(copy_pdf,cell_spot_matrix)

RNA=zeros(length(cell_spot_matrix),1);  %this is where the calucated RNA # is saved
l = length(cell_spot_matrix);
p=zeros(1,1000); 
P=zeros(1,1000);
n=1:1000;

%for sgrS
for i=1:l %along at ith cell
    if cell_spot_matrix(i,3)==0  %no spot in the cell, then # of RNA=0
        RNA(i,1)=0;
    else
    p(:)= copy_pdf(int32(cell_spot_matrix(i,3)),n); % at N(i,1)th row, going through 
                               %every column
                               %sum(p(:)) is the denominator 
                               %of eq (5) in 1.6.3
    P(:)=p(:)/sum(p(:));
    RNA(i,1)=sum(P.*n);  %multiplication of two vectors, eq(5) in 1.6.3
    end    
end


    
       
        