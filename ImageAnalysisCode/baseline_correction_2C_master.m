for i=1:6
    filename = ['Nofspot_per_cell_' int2str((i-1)*2) '_set1.dat'];
    N = load(filename);
    cor=baseline_correction_2C(N);
    filename=sprintf('Nofspot_per_cell_%d_set1_BC.dat',(i-1)*2); 
    dlmwrite(filename,cor,' ');
    clear N cor
    
    filename = ['Nofspot_per_cell_' int2str((i-1)*2) '_set2.dat'];
    N = load(filename);
    cor=baseline_correction_2C(N);
    filename=sprintf('Nofspot_per_cell_%d_set2_BC.dat',(i-1)*2); 
    dlmwrite(filename,cor,' ');
    clear N cor
end