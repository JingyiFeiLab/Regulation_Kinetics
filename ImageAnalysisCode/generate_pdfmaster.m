function copy_pdf = generate_pdfmaster(params)
copy_pdf=zeros(10001, 1000);

x=linspace(0,10000, 10001);
for i=1:1000
    copy_pdf(:,i)= nbinpdf(x, params(1)*i, params(2));
end
copy_pdf;
end