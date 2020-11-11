for i=1:6
    filename1 = ['RNA_per_cell_' int2str((i-1)*2) '_set1_BC.dat'];
    N1 = load(filename1);
    filename2 = ['complex_' int2str((i-1)*2) '_set1.dat'];
    N2 = load(filename2); 
    N1(:,3)=N2-0.000137*N1(:,1); %baseline corrected percentage
    select = find(N1(:,3)<0);
    N1(select,3)=0;
    N1(:,4)=N1(:,2).*N1(:,3); %SP copy number
    N1(:,5)=N1(:,1)-N1(:,4); %sgrS copy - SP
    N1(:,6)=N1(:,2)-N1(:,4); %ptsG copy-SP
    N1(:,1)=round(N1(:,1));
    N1(:,2)=round(N1(:,2));
    N1(:,4)=round(N1(:,4));
    N1(:,5)=round(N1(:,5));
    N1(:,6)=round(N1(:,6));
    filename=sprintf('RNA_per_cell_SPadded_%d_set1.dat',(i-1)*2); 
    dlmwrite(filename,N1,' ');
    SET1=N1;
    clear N1 N2
    
    filename1 = ['RNA_per_cell_' int2str((i-1)*2) '_set2_BC.dat'];
    N1 = load(filename1);
    filename2 = ['complex_' int2str((i-1)*2) '_set2.dat'];
    N2 = load(filename2); 
    N1(:,3)=N2-0.000137*N1(:,1);
    select = find(N1(:,3)<0);
    N1(select,3)=0;
    N1(:,4)=N1(:,2).*N1(:,3);
    N1(:,5)=N1(:,1)-N1(:,4);
    N1(:,6)=N1(:,2)-N1(:,4);
    N1(:,1)=round(N1(:,1));
    N1(:,2)=round(N1(:,2));
    N1(:,4)=round(N1(:,4));
    N1(:,5)=round(N1(:,5));
    N1(:,6)=round(N1(:,6));
    filename=sprintf('RNA_per_cell_SPadded_%d_set2.dat',(i-1)*2); 
    dlmwrite(filename,N1,' ');
    SET2=N1;
    clear N1 N2
    
    Total=[SET1;SET2];
    filename=sprintf('RNA_per_cell_SPadded_%d_total.dat',(i-1)*2); 
    dlmwrite(filename,Total,' ');
    
    ave1=mean(SET1(:,3));
    ave2=mean(SET2(:,3));
    ave(i,1:2)=[ave1 ave2];
    ave(i,3)=mean(ave(i, 1:2));
    ave(i,4)=std(ave(i, 1:2));
    
end
save('SP_cell_ave.dat', 'ave', '-ascii');