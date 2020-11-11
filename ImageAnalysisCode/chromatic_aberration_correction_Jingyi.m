%the purpose of this code is to correct the chromatic aberration
%happening in Nikon microscope - Aug10,2016, Seongjin Park

path = 'E:\Windows\System32\Data\Aug8_2016_beads_two_color_correction\3D_calibration';


[x, y, z, c] =textread('output_two.txt', '%f %f %f %f'); %reading the file which contains (x,y)
leng=length([x, y, x, c]);

corrected=zeros(leng, 4);


for i=1:1:leng
    if c(i)==3
        corrected(i,1)=x(i)-(-15.253+0.00133*x(i));

        corrected(i,2)=y(i)-(-12.13189+(7.97552*10^-4)*y(i));
    else
        corrected(i,1)=x(i);

        corrected(i,2)=y(i);
        
    end
    
        corrected(i,3)=z(i);
        corrected(i,4)=c(i);
    
end

dlmwrite('corrected_output_A568.txt',corrected,' ');