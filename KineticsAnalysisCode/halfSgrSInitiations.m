folderTitle = '/Users/reyer/Data/SingleCellEpi/MR156/halfSgrS73020/';

b_e = 1.33E-03; %std = 2.23E-3
be_sd = 0.0010905149974;

a_s_wt = 0.242; %a_s = 2.29E+03;
a_s_wt_sd = 0.02771858582;
a_s_rne = 0.37649; %a_s = 2.29E+03; 
a_s_rne_sd = 0.05064298767;
b_s_wt = 1.46E-03; 
b_s_wt_sd = 0.0002935200249;
b_s_rne = 0.85E-03; 
b_s_rne_sd = 0.0001349089028;

b_m = 0.003229666667;
b_m_sd = 0.00007212489168;

k_on = 1.34E-03;
k_on_sd = 0.0006514260861;
k_off = 3.03E-01;
k_off_sd = 0.04203270155;

k_nuc = 1.34E-03;
k_nuc_sd = 0.0006514260861;
k_off_nuc = 3.03E-01;
k_off_nuc_sd = 0.04203270155;

b_ms = 0.004487666667; % b_ms = 4.88E-03;
b_ms_sd = 0.001467482561;
b_p = .00018;

b_nuc_wt = b_e + b_ms;
b_nuc_wt_sd = 0.001594820993;
b_nuc_rne = b_ms;
b_nuc_rne_sd = 0;

kx_wt_minus =  14.0582;
kx_wt_plus = 13.72; % kx_wt = 0.0040191;
kxw_sd = 0.7665063622;

kx_rne_minus = 6.12195;
kx_rne_plus = 5.91575; %kx_rne = 0.0023639;
kxr_sd = 0.3636650176;

kx_s = 0.1395; %std = .15
kxs_sd = 0.03322472724;

k_init_wt = 0.08996; %std = 310.41
k_init_wt_sd = 0.01332896283;
k_init_rne = 0.06525; %std = 1164.1
k_init_rne_sd = 0.005303300859;

k_elon = 0.0643;
k_elon_prime = 2.18E-02;
k_elon_prime_sd = .0001;

n_plasmids = 20;


close all
tv = linspace(0, 1800);
tvp = linspace(0, 1440);
tvm = linspace(0,1620);
%tvm = linspace(time(1), max(time)-time(3));

timex = [0,60,180,360,720,1080,1440,1800];

%Protein
%mRNA
%sRNA


molecules_wt_plus1 = [63299.71435,84141.76378,92402.08862,94059.45704,273371.2844,879098.9228,1331884.175,1607995.719;
                      16451.1872,408600.674,859039.2094,1042661.409,725810.8373,1077776.2,1948487.552,1571600.626;
                      915553.04,969278.3673,692863.6179,640696.8966,352026.9624,450922.4441,339613.4235,405109.419];
                  
molecules_wt_plus2 = [78020.12364,60847.61886,50041.19003,112402.7926,485958.8642,679304.8251,879693.1022,1430223.153;
                      17396.48268,452189.5981,802101.0655,1142454.893,1288669.698,1393308.53,1651538.696,1885986.362;
                      1016205.565,1028096.508,731813.9828,659914.1046,543955.6785,484966.5779,460182.1843,426409.5266];                  
                  
              
moleculeswt_plus = [mean((([molecules_wt_plus1(3,1),molecules_wt_plus2(3,1)])+5.2879e+04)/7.3302e+03) 0 0 mean(((([molecules_wt_plus1(2,1),molecules_wt_plus2(2,1)])+1.3257e+05)/6.0368e+03)+1)/20 0 mean([molecules_wt_plus1(1,4),molecules_wt_plus2(1,4)])];
std_wt_plus = [std((([molecules_wt_plus1(3,1),molecules_wt_plus2(3,1)])+5.2879e+04)/7.3302e+03) 0 0 std(((([molecules_wt_plus1(2,1),molecules_wt_plus2(2,1)])+1.3257e+05)/6.0368e+03)+1) 0 std([molecules_wt_plus1(1,4),molecules_wt_plus2(1,4)])];


wtplusp = zeros(100,120);
wtplusm = zeros(100,120);
wtpluss = zeros(100,120);

wtplus_pspread = zeros(100,1);
wtplus_mspread = zeros(100,1);
wtplus_sspread = zeros(100,1);



wtplusp(:,1:6) = bestInit_fit(tvp,a_s_wt,b_s_wt,k_init_wt,k_elon,k_elon_prime,b_m,k_on,k_off,kx_wt_plus,kx_s*kx_wt_plus,b_ms,b_e,b_p,k_nuc,b_nuc_wt,n_plasmids,k_off_nuc,moleculeswt_plus);
wtplusm(:,1:6) = bestInit_fit(tv,a_s_wt,b_s_wt,k_init_wt,k_elon,k_elon_prime,b_m,k_on,k_off,kx_wt_plus,kx_s*kx_wt_plus,b_ms,b_e,b_p,k_nuc,b_nuc_wt,n_plasmids,k_off_nuc,moleculeswt_plus);


i_rand = 7:12;
i_rand_check = [];
for i = 1:19
    
    
    
    a3 = bestInit_fit(tvp,normrnd(a_s_wt,a_s_wt_sd),normrnd(b_s_wt,b_s_wt_sd),normrnd(k_init_wt,k_init_wt_sd),k_elon,normrnd(k_elon_prime,k_elon_prime_sd),normrnd(b_m,b_m_sd),normrnd(k_on,k_on_sd),normrnd(k_off,k_off_sd),normrnd(kx_wt_plus,kxr_sd),normrnd(kx_s,kxs_sd)*kx_wt_plus,normrnd(b_ms,b_ms_sd),normrnd(b_e,be_sd),b_p,normrnd(k_nuc,k_nuc_sd),normrnd(b_nuc_wt,b_nuc_wt_sd),n_plasmids,normrnd(k_off_nuc,k_off_nuc_sd),normrnd(moleculeswt_plus,std_wt_plus));
    if length(a3(:,1)) == 100
        wtplusp(:,i_rand) = a3;
    else
        i_rand_check = [i_rand_check i_rand];
    end
    
    a4 = bestInit_fit(tv,normrnd(a_s_wt,a_s_wt_sd),normrnd(b_s_wt,b_s_wt_sd),normrnd(k_init_wt,k_init_wt_sd),k_elon,normrnd(k_elon_prime,k_elon_prime_sd),normrnd(b_m,b_m_sd),normrnd(k_on,k_on_sd),normrnd(k_off,k_off_sd),normrnd(kx_wt_plus,kxr_sd),normrnd(kx_s,kxs_sd)*kx_wt_plus,normrnd(b_ms,b_ms_sd),normrnd(b_e,be_sd),b_p,normrnd(k_nuc,k_nuc_sd),normrnd(b_nuc_wt,b_nuc_wt_sd),n_plasmids,normrnd(k_off_nuc,k_off_nuc_sd),normrnd(moleculeswt_plus,std_wt_plus));
    if length(a4(:,1)) == 100
        wtplusm(:,i_rand) = a4;
    else
        i_rand_check = [i_rand_check i_rand];
    end
   
    i_rand = i_rand + 6;
end
    

wt_plus_protein_error = zeros(100,20);

wt_plus_mRNA_error = zeros(100,20);

wt_plus_sRNA_error = zeros(100,20);

for i = 1:100
    ij = 1;
    for j = 4:6:120
        
        wt_plus_mRNA_error(i,ij) = wtplusm(i,j) + wtplusm(i,j+1);
        ij = ij + 1;
        
    end

    ik = 1;
    for k = 6:6:120
        
        wt_plus_protein_error(i,k) = wtplusp(i,k);
        ik = ik + 1;
    end
    
    il = 1;
    for j = 1:6:120
        
        
        wt_plus_sRNA_error(i,il) = wtplusm(i,j) + wtplusm(i,j+2) + wtplusm(i,j+4);
        il = il + 1;
        
    end
    
end

for i = 1:100
    
    
    percentiles = [];
    nonOutliers = [];
    percentiles = prctile(wt_plus_protein_error(i,:),[5,95]);
    outlierIndexes = wt_plus_protein_error(i,:) < percentiles(1) | wt_plus_protein_error(i,:) > percentiles(2);
    nonOutliers = wt_plus_protein_error(i,~outlierIndexes);
    
    wtplus_pspread(i) = std(nonOutliers);
    
    percentiles = [];
    nonOutliers = [];
    percentiles = prctile(wt_plus_mRNA_error(i,:),[5,95]);
    outlierIndexes = wt_plus_mRNA_error(i,:) < percentiles(1) | wt_plus_mRNA_error(i,:) > percentiles(2);
    nonOutliers = wt_plus_mRNA_error(i,~outlierIndexes);
    
    
    wtplus_mspread(i) = std(nonOutliers);
    
    percentiles = [];
    nonOutliers = [];
    percentiles = prctile(wt_plus_sRNA_error(i,:),[5,95]);
    outlierIndexes = wt_plus_sRNA_error(i,:) < percentiles(1) | wt_plus_sRNA_error(i,:) > percentiles(2);
    nonOutliers = wt_plus_sRNA_error(i,~outlierIndexes);
    
    
    wtplus_sspread(i) = std(nonOutliers);
    
    
end  

% wtplus_pspread = sgolayfilt(wtplus_pspread,3,11);
% wtplus_mspread = sgolayfilt(wtplus_mspread,3,11);
% wtplus_sspread = sgolayfilt(wtplus_sspread,3,11);


molecules_wt_plus = zeros(6,8);

for i = 1:8
     
    molecules_wt_plus(1,i) = mean([molecules_wt_plus1(1,i),molecules_wt_plus2(1,i)]);
     
    molecules_wt_plus(2,i) = abs(diff([molecules_wt_plus1(1,i),molecules_wt_plus2(1,i)]));
    
    molecules_wt_plus(3,i) = mean((([molecules_wt_plus1(2,i),molecules_wt_plus2(2,i)]+1.3257e+05)/6.0368e+03)+1);
     
    molecules_wt_plus(4,i) = abs(((diff([molecules_wt_plus1(2,i),molecules_wt_plus2(2,i)])+1.3257e+05)/6.0368e+03)+1);
     
    molecules_wt_plus(5,i) = mean((([molecules_wt_plus1(3,i),molecules_wt_plus2(3,i)]+5.2879e+04)/7.3302e+03));
     
    molecules_wt_plus(6,i) = abs(((diff([molecules_wt_plus1(3,i),molecules_wt_plus2(3,i)])+5.2879e+04)/7.3302e+03));
    
    
end


figure(2)
e3 = errorbar(timex,molecules_wt_plus(3,:),molecules_wt_plus(4,:),'s');
e3.MarkerFaceColor = [0,.3,0];
e3.MarkerSize = 5;
e3.Color = [0,.3,0];
e3.LineWidth = 2;
hold on

hold on
shadedErrorBar(tv,(wtplusm(:,4)+wtplusm(:,5)),wtplus_mspread,'lineProps', {'Color',[0,.40,0]})
hold on
plot(tv,(wtplusm(:,4)+wtplusm(:,5)),'LineWidth',1,'Color','k')
hold on

e3 = errorbar(timex,molecules_wt_plus(3,:),molecules_wt_plus(4,:),'s');
e3.MarkerFaceColor = [0,.3,0];
e3.MarkerSize = 5;
e3.Color = [0,.3,0];
e3.LineWidth = 2;
hold on

lgd = legend('WT SgrS');
lgd.FontSize = 10;
lgd.FontWeight = 'bold';
lgd.Location = 'northwest';
title('WT{\it ptsG-sfGFP}, {0.25% \alphaMG}','FontSize',14,'FontWeight','bold','Color',[0,.40,0],'FontName','Arial')
xlabel('Time after Induction (min)','FontSize',18,'FontName','Arial')
ylabel('Copy Number','FontSize',18,'FontName','Arial')
set(gca,'YLim',[-100 2.5E6/3025.7])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10,'FontWeight','bold')
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[626,281,248,201])
file1 = strcat([folderTitle,'WTptsGmRNA_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'WTptsGmRNA_fit.fig']);
savefig(gcf,file1_fig)

figure(4)
e7 = errorbar(timex,molecules_wt_plus(1,:),molecules_wt_plus(2,:),'s');
e7.MarkerFaceColor = [.05,.05,.53];
e7.MarkerSize = 5;
e7.Color = [.05,.05,.53];
e7.LineWidth = 2;
hold on

shadedErrorBar(tvp+360,wtplusp(:,6),wtplus_pspread,'lineProps', {'Color',[.05,.05,.38]})
hold on
plot(tvp+360,wtplusp(:,6),'LineWidth',1,'Color','k')

e7 = errorbar(timex,molecules_wt_plus(1,:),molecules_wt_plus(2,:),'s');
e7.MarkerFaceColor = [.05,.05,.53];
e7.MarkerSize = 5;
e7.Color = [.05,.05,.53];
e7.LineWidth = 2;

lgd = legend('WT SgrS');
lgd.FontSize = 10;
lgd.FontWeight = 'bold';
lgd.Location = 'northwest';
title('WT sfGFP, {0.25% \alphaMG}','FontSize',14,'FontWeight','bold','Color','b')
xlabel('Time after Induction (min)','FontSize',12,'FontWeight','bold')
ylabel('Fluorescence (A.U)','FontSize',12,'FontWeight','bold')
set(gca,'YLim',[-10000 8E6])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10,'FontWeight','bold')
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[626,281,248,201])
file1 = strcat([folderTitle,'WTptsGProtein_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'WTptsGProtein_fit.fig']);
savefig(gcf,file1_fig)

figure(6)
e11 = errorbar(timex,molecules_wt_plus(5,:),molecules_wt_plus(6,:),'s');
e11.MarkerFaceColor = [1,.75,.79];
e11.MarkerSize = 5;
e11.Color = [1,.75,.79];
e11.LineWidth = 2;
hold on

hold on
shadedErrorBar(tv,wtplusm(:,1)+wtplusm(:,3)+wtplusm(:,5),wtplus_sspread,'lineProps', {'Color',[1,.42,.71]})
hold on
plot(tv,wtplusm(:,1)+wtplusm(:,3)+wtplusm(:,5),'LineWidth',1,'Color','k')

e11 = errorbar(timex,molecules_wt_plus(5,:),molecules_wt_plus(6,:),'s');
e11.MarkerFaceColor = [1,.75,.79];
e11.MarkerSize = 5;
e11.Color = [1,.75,.79];
e3.LineWidth = 2;
hold on

lgd = legend('WT SgrS');
lgd.FontSize = 10;
lgd.FontWeight = 'bold';
lgd.Location = 'northwest';
title('WT SgrS, {0.25% \alphaMG} ','FontSize',18,'FontWeight','bold','Color','r')
xlabel('Time after Induction (min)','FontSize',18,'FontName','Arial')
ylabel('Copy Number','FontSize',18,'FontName','Arial')
set(gca,'YLim',[-10 3.0E6/7096.25])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10,'FontWeight','bold')
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[626,281,248,201])
file1 = strcat([folderTitle,'WTptsGsRNA_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'WTptsGsRNA_fit.fig']);
savefig(gcf,file1_fig)