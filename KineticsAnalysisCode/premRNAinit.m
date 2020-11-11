close all
folderTitle = '/Users/reyer/Data/SingleCellEpi/MR156/premRNA51820/';

b_e = 0.0157; %std = 2.23E-3
be_sd = 1E-7;

a_s_wt = 0.305475; %a_s = 2.29E+03;
a_s_wt_sd = 0.006187184335;
a_s_rne = 0.361325; %a_s = 2.29E+03; 
a_s_rne_sd = 0.01586747617;
b_s_wt = 1.46E-03; 
b_s_wt_sd = 0.0002935200249;
b_s_rne = 0.85E-03; 
b_s_rne_sd = 0.0001349089028;

b_m = 0.003229666667;
b_m_sd = 0.00007212489168;

k_on = 0.000254;
k_on_sd = 0;
k_off = 0.0292;
k_off_sd = 0;

k_nuc = 0.00624;
k_nuc_sd = 0;
k_off_nuc = 0.01081;
k_off_nuc_sd = 0;

b_ms = 0.009078; % b_ms = 4.88E-03;
b_ms_sd = 0;
b_p = .00018;

b_nuc_wt = b_e + b_ms;
b_nuc_wt_sd = 0;
b_nuc_rne = b_ms;
b_nuc_rne_sd = 0;

kx_wt_minus =  3.650233333;
kx_wt_plus = 3.65; % kx_wt = 0.0040191;
kxw_sd = 0.7665063622;

kx_rne_minus = 2.5366;
kx_rne_plus = 2.16075; %kx_rne = 0.0023639;
kxr_sd = 0.9100464274;

kx_s = 0.016; %std = .15
kxs_sd = 0.068;

k_init_wt = 0.2192666667; %std = 310.41
k_init_wt_sd = 0.009295339334;
k_init_rne = 0.2165; %std = 1164.1
k_init_rne_sd = 0.01909188309;

k_elon = 0.0743;
k_elon_prime = 0.0218;
k_elon_prime_sd = .0001;

n_plasmids = 20;



close all
tv = linspace(0, 2520);
tvp = linspace(0, 2160);
tvm = linspace(0,2520);

timex = [0,360,720,1080,1440,1800,2160,2520];






wpp1 = [4893286.68
8906893.953
10632693.03
11936906.68
15594274.71
16179536.96
16037966.92
14967956.67];

wpp2 = [6109745.756
8789932.785
10997674.36
12458214.45
13373057.81
16732476.06
16587509.59
15798323.27];

wpp3 = [5208002.736
7807893.587
9536727.262
12178635.4
12550445.12
14109669.44
14016284.97
14699703.79];


wmp1 = [3435632.493
3900085.168
4040177.559
3673281.347
1942951.852
1490869.54
2189226.711
2395907.563];

wmp2 = [1535724.849
1640236.135
1226712.194
999647.2196
1167564.493
1006417.939
625753.2134
714870.1859];

wmp3 = [2618705.985
2673116.879
2499702.136
2103895.669
1112001.644
1085069.313
593656.5789
491957.996];


wsp1 = [20409.87246
271104.5671
273452.7785
311335.4054
425947.5461
421265.5618
396374.4361
383874.6717];

wsp2 = [28193.16431
163545.9793
156515.6786
212156.765
275338.023
308021.621
197243.1084
221072.7948];

wsp3 = [33843.95136
240669.2611
285195.0278
404143.1438
246759.5766
342928.9136
284345.1964
274038.4447];



molecules_wt_plus = zeros(6,8);
for i = 1:8
    

    molecules_wt_plus(1,i) = mean([wpp1(i),wpp2(i),wpp3(i)]);
    
    
     
    molecules_wt_plus(2,i) = std([wpp1(i),wpp2(i),wpp3(i)])/1;
    
    
    
    molecules_wt_plus(3,i) = mean(((([wmp1(i),wmp2(i),wmp3(i)])+4.3472e+04)/1.6831e+03)+1);
    
     
    molecules_wt_plus(4,i) = std(((([wmp1(i),wmp2(i),wmp3(i)])+4.3472e+04)/1.6831e+03)+1)/2.25;
    
    
    molecules_wt_plus(5,i) = mean((([wsp1(i),wsp2(i),wsp3(i)])+1.7582e+04)/7.8067e+03);
    
    
     
    molecules_wt_plus(6,i) = std((([wsp1(i),wsp2(i),wsp3(i)])+1.7582e+04)/7.8067e+03)/1;
    
end

moleculeswt_plus = [0 0 molecules_wt_plus(5,1) molecules_wt_plus(3,1) 0 molecules_wt_plus(1,2)];
std_wt_plus = [ 0 0 molecules_wt_plus(6,1) molecules_wt_plus(4,1) 0 molecules_wt_plus(2,2)];


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
    
    wtplus_mspread(i) = std(nonOutliers)/1;
    
    percentiles = [];
    nonOutliers = [];
    percentiles = prctile(wt_plus_sRNA_error(i,:),[5,95]);
    outlierIndexes = wt_plus_sRNA_error(i,:) < percentiles(1) | wt_plus_sRNA_error(i,:) > percentiles(2);
    nonOutliers = wt_plus_sRNA_error(i,~outlierIndexes);
    
    wtplus_sspread(i) = std(nonOutliers)/1;
    
    
end  


figure(2)
e3 = errorbar(timex,molecules_wt_plus(3,:),molecules_wt_plus(4,:),'s');
e3.MarkerFaceColor = [0,.3,0];
e3.MarkerSize = 5;
e3.Color = [0,.3,0];
e3.LineWidth = 2;
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
title('WT{\it ptsG-sfGFP}','FontSize',18,'FontWeight','bold','Color','g','FontName','Arial')
xlabel('Time after Induction (min)','FontSize',18,'FontName','Arial')
ylabel('Copy Number','FontSize',18,'FontName','Arial')
set(gca,'YLim',[-100 5.5E6/1125.7])
set(gca,'XLim',[-1 2650])
set(gca,'XTick',[0,360,720,1080,1440,1800,2160,2520])
set(gca,'XTickLabel',[0,6,12,18,24,30,36,42],'FontSize',10,'FontWeight','bold')
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
hold on

e7 = errorbar(timex,molecules_wt_plus(1,:),molecules_wt_plus(2,:),'s');
e7.MarkerFaceColor = [.05,.05,.53];
e7.MarkerSize = 5;
e7.Color = [.05,.05,.53];
e7.LineWidth = 2;
hold on

lgd = legend('WT SgrS');
lgd.FontSize = 10;
lgd.FontWeight = 'bold';
lgd.Location = 'northwest';
title('WT sfGFP','FontSize',14,'FontWeight','bold','Color','b')
xlabel('Time after Induction (min)','FontSize',12,'FontWeight','bold')
ylabel('Fluorescence (A.U)','FontSize',12,'FontWeight','bold')
set(gca,'YLim',[-10000 20E6])
set(gca,'XLim',[-1 2650])
set(gca,'XTick',[0,360,720,1080,1440,1800,2160,2520])
set(gca,'XTickLabel',[0,6,12,18,24,30,36,42],'FontSize',10,'FontWeight','bold')
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

shadedErrorBar(tv,wtplusm(:,1)+wtplusm(:,3)+wtplusm(:,5),wtplus_sspread,'lineProps', {'Color',[1,.42,.71]})
hold on
plot(tv,wtplusm(:,1)+wtplusm(:,3)+wtplusm(:,5),'LineWidth',1,'Color','k')
hold on

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
title('WT SgrS','FontSize',18,'FontWeight','bold','Color','r')
xlabel('Time after Induction (min)','FontSize',18,'FontName','Arial')
ylabel('Copy Number','FontSize',18,'FontName','Arial')
set(gca,'YLim',[-10 1E6/8096.25])
set(gca,'XLim',[-1 2650])
set(gca,'XTick',[0,360,720,1080,1440,1800,2160,2520])
set(gca,'XTickLabel',[0,6,12,18,24,30,36,42],'FontSize',10,'FontWeight','bold')
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

figure(7)
plot(linspace(0,1440,100),wtplusm(:,1),'LineWidth',3); 
title('Free sRNA','FontSize',18);
file1 = strcat([folderTitle,'sRNA']);
set(gcf,'PaperPositionMode','auto');
print(file1,'-dpng','-r0')

figure(8)
plot(linspace(0,1440,100),wtplusm(:,2),'LineWidth',3); 
title('Elongating mRNA','FontSize',18);
file1 = strcat([folderTitle,'elong_mRNA']);
set(gcf,'PaperPositionMode','auto');
print(file1,'-dpng','-r0')

figure(9)
plot(linspace(0,1440,100),wtplusm(:,3),'LineWidth',3); 
title('Co-Transcriptionally Bound mRNA','FontSize',18);
file1 = strcat([folderTitle,'cobound_mRNA']);
set(gcf,'PaperPositionMode','auto');
print(file1,'-dpng','-r0')

figure(10)
plot(linspace(0,1440,100),wtplusm(:,4),'LineWidth',3); 
title('Full mRNA','FontSize',18);
file1 = strcat([folderTitle,'full_mRNA']);
set(gcf,'PaperPositionMode','auto');
print(file1,'-dpng','-r0')

figure(11)
plot(linspace(0,1440,100),wtplusm(:,5),'LineWidth',3); 
title('Post sRNA-mRNA Complex','FontSize',18);
file1 = strcat([folderTitle,'s_mRNA']);
set(gcf,'PaperPositionMode','auto');
print(file1,'-dpng','-r0')

figure(12)
plot(linspace(0,1440,100),wtplusm(:,6),'LineWidth',3); 
title('Protein','FontSize',18);
file1 = strcat([folderTitle,'protein']);
set(gcf,'PaperPositionMode','auto');
print(file1,'-dpng','-r0')

copy_array = [k_on,k_off,kx_s,b_e,b_ms,k_elon_prime,k_init_wt,k_init_rne,kx_wt_plus,kx_rne_plus,b_m,kx_wt_minus,kx_rne_minus,k_elon,k_nuc,b_nuc_wt,b_nuc_rne,k_off_nuc;
    k_on_sd,k_off_sd,kxs_sd,be_sd,b_ms_sd,k_elon_prime_sd,k_init_wt_sd,k_init_rne_sd,kxw_sd,kxr_sd,b_m_sd,kxw_sd,kxr_sd,0,k_nuc_sd,b_nuc_wt_sd,b_nuc_rne_sd,k_off_nuc_sd];

copy_table = array2table(copy_array);copy_table.Properties.VariableNames = {'k_on' 'k_off' 'kx_s' 'b_e' 'b_ms' 'k_elon_prime' 'k_init_wt' 'k_init_rne_' 'kx_wt_plus' 'kx_rne_plus' 'b_m' 'kx_wt_minus' 'kx_rne_minus' 'k_elon' 'k_nuc' 'b_nuc_wt' 'b_nuc_rne' 'k_off_nuc'};

copy_table_file = strcat([folderTitle,'parameters.csv']);

writetable(copy_table,copy_table_file);