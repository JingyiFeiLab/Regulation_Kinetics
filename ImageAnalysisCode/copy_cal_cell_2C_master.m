pdf_1=load('sgrS_nbin_vector.dat');
pdf_2=load('ptsG_nbin_vector.dat');
A1=zeros(6, 6);
A2=zeros(6, 6);
A=zeros(6, 6);

for i=1:6
    filename = ['Nofspot_per_cell_' int2str((i-1)*2) '_set1_BC.dat'];
    N = load(filename);
    RNA1 = copy_cal_cell_2C_2(N, pdf_1, pdf_2);
    filename=sprintf('RNA_per_cell_%d_set1_BC.dat',(i-1)*2); 
    dlmwrite(filename,RNA1,' ');
    clear N;
    A1(i,1)= mean(RNA1(:,1));
    A1(i,2)=std(RNA1(:,1))/sqrt(length(RNA1));
    A1(i,3)= var(RNA1(:,1));
    A1(i,4)= mean(RNA1(:,2));
    A1(i,5)=std(RNA1(:,2))/sqrt(length(RNA1));
    A1(i,6)= var(RNA1(:,2));
    
    
    filename = ['Nofspot_per_cell_' int2str((i-1)*2) '_set2_BC.dat'];
    N = load(filename);
    RNA2 = copy_cal_cell_2C_2(N, pdf_1, pdf_2);
    filename=sprintf('RNA_per_cell_%d_set2_BC.dat',(i-1)*2); 
    dlmwrite(filename,RNA2,' ');
    clear N;
    A2(i,1)= mean(RNA2(:,1));
    A2(i,2)=std(RNA2(:,1))/sqrt(length(RNA2));
    A2(i,3)= var(RNA2(:,1));
    A2(i,4)= mean(RNA2(:,2));
    A2(i,5)=std(RNA2(:,2))/sqrt(length(RNA2));
    A2(i,6)= var(RNA2(:,2));
    
    RNA=[RNA1; RNA2];
    filename=sprintf('RNA_per_cell_%d_BC.dat',(i-1)*2); 
    dlmwrite(filename,RNA,' ');   
    A(i,1)= mean(RNA(:,1));
    A(i,2)=std(RNA(:,1))/sqrt(length(RNA));
    A(i,3)= var(RNA(:,1));
    A(i,4)= mean(RNA(:,2));
    A(i,5)=std(RNA(:,2))/sqrt(length(RNA));
    A(i,6)= var(RNA(:,2));
    clear RNA;
    
    
    
end
 save ('copy_new_BC_mean_std_set1.dat', 'A1', '-ascii');
 save ('copy_new_BC_mean_std_set2.dat', 'A2', '-ascii');
 save ('copy_new_BC_mean_std.dat', 'A', '-ascii'); 
    
    
   


