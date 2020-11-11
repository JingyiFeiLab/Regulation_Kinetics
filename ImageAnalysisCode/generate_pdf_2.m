%c1=load('nbin_1R.dat');
pdf=zeros(10001, 1000);
x=linspace(0,10000, 10001);
for i=1:1000
    pdf(:,i)= nbinpdf(x, 0.9959*i, 0.0217);
end