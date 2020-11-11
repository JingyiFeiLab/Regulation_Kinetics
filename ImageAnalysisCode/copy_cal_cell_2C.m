clear;
pdf_1=load('sgrS_nbin_vector.dat');
pdf_2=load('ptsG_nbin_vector.dat');
N=load('N_of_spots_per_cell_0_BC.dat');
RNA=zeros(length(N),2);
p=zeros(1,1000);
P=zeros(1,1000);
n=1:1:1000;

%for sgrS
for i=1:length(N)
    if N(i,1)==0
        RNA(i,1)=0;
    else
    p(:)= pdf_1(N(i,1), 1:length(n));
    P(:)=p(:)/sum(p(:));
    RNA(i,1)=sum(P.*n);
    end    
end

%for ptsG
for i=1:length(N)
    if N(i,2)==0
        RNA(i,2)=0;
    else
    p(:)= pdf_2(N(i,2), 1:length(n));
    P(:)=p(:)/sum(p(:));
    RNA(i,2)=sum(P.*n);
    end    
end

    
       
        