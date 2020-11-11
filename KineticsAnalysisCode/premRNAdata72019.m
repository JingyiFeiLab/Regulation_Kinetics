close all
folderTitle = '/Users/reyer/Data/SingleCellEpi/MR156/premRNA62720/';

b_e = 1.90E-03; %std = 2.23E-3
be_sd = 0.0003917371568;
a_m_plus_wt = [6.8318,5.5318,5.5318,5.5318];% std =589.18; 
amp_wt_sd = 0.7506645589;
a_m_plus_rne = 5.8198; %std = 796.64
amp_rne_sd = 0.2046367025;
a_s_wt = 0.29845; %a_s = 2.29E+03;
a_s_wt_sd = 0.09369164851;
a_s_rne = 0.4121; %a_s = 2.29E+03; 
a_s_rne_sd = 0.04369919908;
b_s_wt = 1.44E-03; 
b_s_wt_sd = 0.00005939696962;
b_s_rne = 0.0000487903679; 
b_s_rne_sd = 0.0001447447581;
b_m = 0.003279;
b_m_sd = 0.00007212489168;
k_on = 9.63E-05;
k_on_sd = 0.0001287245469;
k_off = 6.97E-02;
k_off_sd = 0.01347745525;
kx_wt_minus = 2.5196;
kx_wt_plus = [3.0196,3.0196,3.0196,3.0196]; % kx_wt = 0.0040191;
kxw_sd = 0.1046518036;
% kx_wt_plus=mode(kx_fit); %kx_wt_plus=0.0033;
kx_rne_minus = 2.1407;
kx_rne_plus = 1.66; %kx_rne = 0.0023639;
kxr_sd = 0.03889087297;
kx_s = 0.5395; %std = .15
kxs_sd = 0.1463711037;
b_ms_rne = 0.00356; % b_ms = 4.88E-03;
b_ms_wt = [0.00356,0.00356,0.00356,0.00356]; % b_ms = 4.88E-03;
b_ms_sd = 0.0001414213562;
b_p = .00018;
a_m_minus_wt = 5.5318; %std = 310.41
amw_sd = 0.7506645589;
a_m_minus_rne = 5.8198; %std = 1164.1
amm_sd = 0.2046367025;


close all
% tv = linspace(0, 1800);
% tvp = linspace(0, 1440);
% tvm = linspace(0,1620);
tv1 = linspace(0, 360);
tvp1 = linspace(0, 360);
tv2 = linspace(0, 360);
tvp2 = linspace(0, 360);
tv3 = linspace(0, 480);
tvp3 = linspace(0, 480);
tv4 = linspace(0, 1320);
tvp4 = linspace(0, 960);

tvm = linspace(0,1620);
%tvm = linspace(time(1), max(time)-time(3));

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
    
    
    
    molecules_wt_plus(3,i) = mean(((([wmp1(i),wmp2(i),wmp3(i)])-20022)/1125.7)+1);
    
     
    molecules_wt_plus(4,i) = std(((([wmp1(i),wmp2(i),wmp3(i)])-20022)/1125.7)+1)/2.25;
    
    
    molecules_wt_plus(5,i) = mean((([wsp1(i),wsp2(i),wsp3(i)])-48376)/8096.5);
    
    
     
    molecules_wt_plus(6,i) = std((([wsp1(i),wsp2(i),wsp3(i)])-48376)/8096.5)/1;
    
end

moleculeswt_plus = [molecules_wt_plus(5,1) molecules_wt_plus(3,1) 0 molecules_wt_plus(1,2)];


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
plot([tv1';(tv2+360)';(tv3+720)';(tv4+1200)'],[wtplusm1(:,2)+wtplusm1(:,3);wtplusm2(:,2)+wtplusm2(:,3);wtplusm3(:,2)+wtplusm3(:,3);wtplusm4(:,2)+wtplusm4(:,3)],'LineWidth',1,'Color','k')

% 
% legend('WT +SgrS','WT -SgrS')
shadedErrorBar([tv1';(tv2+360)';(tv3+720)';(tv4+1200)'],[wtplusm1(:,2)+wtplusm1(:,3);wtplusm2(:,2)+wtplusm2(:,3);wtplusm3(:,2)+wtplusm3(:,3);wtplusm4(:,2)+wtplusm4(:,3)],wtplus_mspread,'lineProps', {'Color',[0,.40,0]})
hold on
plot([tv1';(tv2+360)';(tv3+720)';(tv4+1200)'],[wtplusm1(:,2)+wtplusm1(:,3);wtplusm2(:,2)+wtplusm2(:,3);wtplusm3(:,2)+wtplusm3(:,3);wtplusm4(:,2)+wtplusm4(:,3)],'LineWidth',1,'Color','k')
hold on

e3 = errorbar(timex,molecules_wt_plus(3,:),molecules_wt_plus(4,:),'s');
e3.MarkerFaceColor = [0,.3,0];
e3.MarkerSize = 5;
e3.Color = [0,.3,0];
e3.LineWidth = 2;
% 
title('Pre-Induced mRNA, {\it ptsG-sfGFP}','FontSize',14,'FontWeight','bold','Color','g')
xlabel('Time after Induction (min)','FontSize',12,'FontWeight','bold')
ylabel('Fluorescence (A.U)','FontSize',12,'FontWeight','bold')
set(gca,'YLim',[0 4.5E6/1125.7])
set(gca,'XLim',[0,2520])
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
plot([(tvp1+360)';(tvp2+360+360)';(tvp3+720+360)';(tvp4+1200+360)'],[wtplusp1(:,4);wtplusp2(:,4);wtplusp3(:,4);wtplusp4(:,4)],'LineWidth',1,'Color','k')

% 
% legend('WT +SgrS','WT -SgrS')
shadedErrorBar([(tvp1+360)';(tvp2+360+360)';(tvp3+720+360)';(tvp4+1200+360)'],[wtplusp1(:,4);wtplusp2(:,4);wtplusp3(:,4);wtplusp4(:,4)],wtplus_pspread,'lineProps', {'Color',[.05,.05,.38]})
hold on
plot([(tvp1+360)';(tvp2+360+360)';(tvp3+720+360)';(tvp4+1200+360)'],[wtplusp1(:,4);wtplusp2(:,4);wtplusp3(:,4);wtplusp4(:,4)],'LineWidth',1,'Color','k')
hold on

e7 = errorbar(timex,molecules_wt_plus(1,:),molecules_wt_plus(2,:),'s');
e7.MarkerFaceColor = [.05,.05,.53];
e7.MarkerSize = 5;
e7.Color = [.05,.05,.53];
e7.LineWidth = 2;
% 
title('Pre-Induced mRNA, sfGFP','FontSize',14,'FontWeight','bold','Color','b')
xlabel('Time after Induction (min)','FontSize',12,'FontWeight','bold')
ylabel('Fluorescence (A.U)','FontSize',12,'FontWeight','bold')
set(gca,'YLim',[0 20E6])
set(gca,'XLim',[0 2520])
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
plot([tv1';(tv2+360)';(tv3+720)';(tv4+1200)'],[wtplusm1(:,1)+wtplusm1(:,3);wtplusm2(:,1)+wtplusm2(:,3);wtplusm3(:,1)+wtplusm3(:,3);wtplusm4(:,1)+wtplusm4(:,3)],'LineWidth',1,'Color','k')

% 
hold on
legend('WT +SgrS')
shadedErrorBar([tv1';(tv2+360)';(tv3+720)';(tv4+1200)'],[wtplusm1(:,1)+wtplusm1(:,3);wtplusm2(:,1)+wtplusm2(:,3);wtplusm3(:,1)+wtplusm3(:,3);wtplusm4(:,1)+wtplusm4(:,3)],wtplus_sspread,'lineProps', {'Color',[1,.42,.71]})
hold on
plot([tv1';(tv2+360)';(tv3+720)';(tv4+1200)'],[wtplusm1(:,1)+wtplusm1(:,3);wtplusm2(:,1)+wtplusm2(:,3);wtplusm3(:,1)+wtplusm3(:,3);wtplusm4(:,1)+wtplusm4(:,3)],'LineWidth',1,'Color','k')

hold on
e11 = errorbar(timex,molecules_wt_plus(5,:),molecules_wt_plus(6,:),'s');
e11.MarkerFaceColor = [1,.75,.79];
e11.MarkerSize = 5;
e11.Color = [1,.75,.79];
e3.LineWidth = 2;
hold on

title('Pre-Induced mRNA, SgrS','FontSize',14,'FontWeight','bold','Color','r')
xlabel('Time after Induction (min)','FontSize',12,'FontWeight','bold')
ylabel('Fluorescence (A.U)','FontSize',12,'FontWeight','bold')
set(gca,'YLim',[0 3.5E6/16096.25])
set(gca,'XLim',[0,2520])
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