folderTitle = '/Users/reyer/Data/SingleCellEpi/params/';

b_e = logspace(-4,-1.9,11); %std = 2.23E-3
be_sd = 0.0009905149974;

a_s_wt = 0.3646; %a_s = 2.29E+03;
a_s_wt_sd = 0.02771858582;
a_s_rne = 0.37649; %a_s = 2.29E+03; 
a_s_rne_sd = 0.05064298767;
b_s_wt = 1.46E-03; 
b_s_wt_sd = 0.0002935200249;
b_s_rne = 0.85E-03; 
b_s_rne_sd = 0.0001349089028;

b_m = 0.003229666667;
b_m_sd = 0.00007212489168;

k_on = logspace(-5,-2.5,11);
k_on_sd = 0.000651426;
k_off = 0.303;
k_off_sd = 0.042032702;

k_nuc = logspace(-4,-2.5,11);
k_nuc_sd = 0.000651426;
k_off_nuc = 0.303;
k_off_nuc_sd = 0.042032702;

b_ms = 0.004487667; % b_ms = 4.88E-03;
b_ms_sd = 0.001467483;
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

kx_s = 0.5; %std = .15
kxs_sd = 0.03322472724;

k_init_wt = 0.08996; %std = 310.41
k_init_wt_sd = 0.01332896283;
k_init_rne = 0.06525; %std = 1164.1
k_init_rne_sd = 0.005303300859;

k_elon = 0.0643;
k_elon_prime = 2.18E-02;
k_elon_prime_sd = .0001;

Pt = linspace(0,1,11);

n_plasmids = 20;

close all
tv = linspace(0, 18000);
tvp = linspace(0, 200000);
tvm = linspace(0,1620);
%tvm = linspace(time(1), max(time)-time(3));

timex = [0,60,180,360,720,1080,1440,1800];

%Protein
%mRNA
%sRNA

molecules_rne_plus1= [364800.2635,336436.0658,341534.1991,380389.8331,390146.4319,691803.3185,1185416.509,1382649.746;
                      67765.94014,81801.67861,178800.4914,958651.5462,1841993.122,2606387.04,2085851.336,2046367.339;
                      917753.0684, 1928070.802, 2287291.285, 1866986.14, 1789976.436, 2234806.001, 2115092.332, 2253681.743];

molecules_rne_plus2= [323450.4837,271414.5415,255897.5024,289121.3282,552624.045,1179711.879,1256733.641,1940635.146;
                      120987.0961,102501.1293,218328.8313,1211313.97,1527640.739,1949001.673,1161004.273,2012571.706;
                      1539983.611,1814226.119,1640451.051,1700309.258,1157810.216,1291511.045,1139946.553,1250870.797];

molecules_rne_plus3= [257171.8156,242054.4163,248115.1836,246027.831,383643.0025,802367.1203,1131430.97,1382649.746;
                      16225.77785,7646.443305,194166.3494,714018.5826,1418654.816,1533096.423,1617435.188,2046367.339;
                      2188990.237,2050903.278,2390768.473,1788277.574,1296136.946,1018016.032,973140.0654, 2253681.743];

molecules_rne_plus4= [64338.27187,36900.77477,47182.78534,52802.52599,298986.6166,731084.9614,1055270.669,1940635.146;
                      28424.46802,39375.74472,237406.821,1027314.934,2003857.809,2042997.469,1758470.575,2012571.706;
                      2547665.889,2079076.002,2324980.48,2178616.805,1944435.717,956261.5957,1099876.149,1250870.797];
                  
                  
molecules_wt_plus1 = [134917.1029, 184452.1771, 180579.3257,216245.1063, 432250.9255, 839879.6454, 1170873.494,1903248.354;
                      22443.71086, 339380.1687, 774276.8132, 989754.3786, 1178364.294, 1168843.171, 1039977.762,1152377.24;
                      1299487.304, 2158748.53, 1461797.642, 1479650.818, 1213052.738, 1132760.79, 930608.6544,743058.2862];
                  
molecules_wt_plus2 = [144202.1426, 333919.142, 149161.7576, 304584.7737, 738063.1426, 808832.9779, 1392651.226,1681171.266;
                      33318.96697, 262752.4716, 627298.9573, 851325.6287, 1159651.072, 1151829.933, 1016580.917,1118273.902;
                      1788327.436, 1787030.495, 1406293.534, 1188055.269, 1088649.355, 825880.4722, 785990.9093,693407.8427];

molecules_rne_minus1 = [84255.81349,95803.25762,116352.407,191341.6814,516093.6138,1269280.982,1693173.504,3620573.777;
                        22559.45533,32680.31201,202125.5704,1548604.718,1870091.512,2429377.867,2163771.481,1579967.962;
                        0,0,0,0,0,0,0,0];
                    
molecules_rne_minus2 = [83829.39758,97323.40231,82115.75258,79705.12663,600125.1929,1447703.631,1907118.6,4021270.604;
                        24304.45749,21222.77473,283979.2055,1260949.779,2468968.323,2861939.516,2436779.396,1460208.711;
                        0,0,0,0,0,0,0,0]; 

molecules_wt_minus1 = [148497.7287,245486.3533,121653.7209,515704.765,1928256.077,3850694.509,5850005.209,6647973.429;
                       20457.03281,138813.3189,1862038.607,3535042.352,3611851.875,2688972.649,3095618.788,2093476.027;
                       0,0,0,0,0,0,0,0];
                   
molecules_wt_minus2 = [139946.6233,156878.962,166032.6582,685644.4883,2295568.63,4960237.252,5549831.123,4276264.083;
                       29649.85822,155400.9189,2041858.881,3431922.71,3002638.619,3547152.203,2794706.764,2667938.955;
                       0,0,0,0,0,0,0,0];  
                   
% molecules_wt_minus3 = [147008.5757,118225.8312,129874.9853,119087.0998,806120.481,1973703.918,3249202.757,4276264.083;
%                        19725.86219,139812.3909,1366620.452,1361329.174,2415583.235,1479284.733,2253799.432,2253799.432;
%                        0,0,0,0,0,0,0,0];
                   

                   
moleculesrne_plus = [mean((([molecules_rne_plus1(3,1),molecules_rne_plus2(3,1),molecules_rne_plus3(3,1),molecules_rne_plus4(3,1)])+5.2879e+04)/7.3302e+03) 0 0 mean(((([molecules_rne_plus1(2,1),molecules_rne_plus2(2,1),molecules_rne_plus3(2,1),molecules_rne_plus4(2,1)])+1.3257e+05)/6.0368e+03)+1) 0 mean([molecules_rne_plus1(1,4),molecules_rne_plus2(1,4),molecules_rne_plus3(1,4),molecules_rne_plus4(1,4)])];
moleculeswt_plus = [mean((([molecules_wt_plus1(3,1),molecules_wt_plus2(3,1)])+5.2879e+04)/7.3302e+03) 0 0 mean(((([molecules_wt_plus1(2,1),molecules_wt_plus2(2,1)])+1.3257e+05)/6.0368e+03)+1) 0 mean([molecules_wt_plus1(1,4),molecules_wt_plus2(1,4)])];
moleculesrne_minus = [0 mean(((([molecules_rne_minus1(2,1),molecules_rne_minus2(2,1)])+1.3257e+05)/6.0368e+03)+1) mean([molecules_rne_minus1(1,4),molecules_rne_minus2(1,4)])];
moleculeswt_minus = [0 mean(((([molecules_wt_minus1(2,1),molecules_wt_minus2(2,1)])+1.3257e+05)/6.0368e+03)+1) mean([molecules_wt_minus1(1,4),molecules_wt_minus2(1,4)])];

std_rne_plus = [std((([molecules_rne_plus1(3,1),molecules_rne_plus2(3,1)])+5.2879e+04)/7.3302e+03) 0 0  std(((([molecules_rne_plus1(2,1),molecules_rne_plus2(2,1)])+1.3257e+05)/6.0368e+03)+1) 0 std([molecules_rne_plus1(1,4),molecules_rne_plus2(1,4)])];
std_wt_plus = [std((([molecules_wt_plus1(3,1),molecules_wt_plus2(3,1)])+5.2879e+04)/7.3302e+03) 0 0 std(((([molecules_wt_plus1(2,1),molecules_wt_plus2(2,1)])+1.3257e+05)/6.0368e+03)+1) 0 std([molecules_wt_plus1(1,4),molecules_wt_plus2(1,4)])];
std_rne_minus = [0 std(((([molecules_rne_minus1(2,1),molecules_rne_minus2(2,1)])+1.3257e+05)/6.0368e+03)+1) std([molecules_rne_minus1(1,4),molecules_rne_minus2(1,4)])];
std_wt_minus = [0 std(((([molecules_wt_minus1(2,1),molecules_wt_minus2(2,1)])+1.3257e+05)/6.0368e+03)+1) std([molecules_wt_minus1(1,4),molecules_wt_minus2(1,4)])];

rneplusp = zeros(100,120);
rneplusm = zeros(100,120);
rnepluss = zeros(100,120);

wtplusp = zeros(100,120);
wtplusm = zeros(100,120);
wtpluss = zeros(100,120);

rneminusp = zeros(100,60);
rneminusm = zeros(100,60);
rneminuss = zeros(100,60);

wtminusp = zeros(100,60);
wtminusm = zeros(100,60);
wtminuss = zeros(100,60);

rneplus_pspread = zeros(100,1);
rneplus_mspread = zeros(100,1);
rneplus_sspread = zeros(100,1);

wtplus_pspread = zeros(100,1);
wtplus_mspread = zeros(100,1);
wtplus_sspread = zeros(100,1);

rneminus_pspread = zeros(100,1);
rneminus_mspread = zeros(100,1);
rneminus_sspread = zeros(100,1);

wtminus_pspread = zeros(100,1);
wtminus_mspread = zeros(100,1);
wtminus_sspread = zeros(100,1);


% rneminusp(:,1:3) = initiation_fit([k_init_rne,kx_rne_minus],tvp,k_elon,b_m,b_p,n_plasmids,moleculesrne_minus);
% rneminusm(:,1:3) = initiation_fit([k_init_rne,kx_rne_minus],tv,k_elon,b_m,b_p,n_plasmids,moleculesrne_minus);
% 
wtminusp(:,1:3) = initiation_fit([k_init_wt,kx_wt_minus],tvp,k_elon,b_m,b_p,n_plasmids,moleculeswt_minus);
wtminusm(:,1:3) = initiation_fit([k_init_wt,kx_wt_minus],tv,k_elon,b_m,b_p,n_plasmids,moleculeswt_minus);

% figure(1)
% plot(tv,(wtplusm(:,4)+wtplusm(:,5)),'LineWidth',1,'Color','k')
% hold on
% plot(tv,wtminusm(:,2),'LineWidth',1,'Color','k')
% 
% figure(2)
% plot(tvp+360,wtplusp(:,6),'LineWidth',1,'Color','k')
% hold on
% plot(tvp+360,wtminusp(:,3),'LineWidth',1,'Color','k')
% 
% figure(3)
% plot(tv,wtplusm(:,1)+wtplusm(:,3)+wtplusm(:,5))
rep = zeros(11,11,11);

for i = 1:11
    for j = 1:11
        for k = 1:11
            
            wtplusp(:,1:6) = bestInit_fit(tvp,a_s_wt,b_s_wt,k_init_wt,k_elon,k_elon_prime,b_m,k_on(j),k_off,kx_wt_plus,kx_s*kx_wt_plus,b_ms,b_e(i),b_p,k_nuc(j),b_nuc_wt,n_plasmids,k_off_nuc,Pt(k),moleculeswt_plus);
            wtplusm(:,1:6) = bestInit_fit(tv,a_s_wt,b_s_wt,k_init_wt,k_elon,k_elon_prime,b_m,k_on(j),k_off,kx_wt_plus,kx_s*kx_wt_plus,b_ms,b_e(i),b_p,k_nuc(j),b_nuc_wt,n_plasmids,k_off_nuc,Pt(k),moleculeswt_plus);

            protein_rep = 1-(wtplusp(100,6)/wtminusp(100,3));
            mRNA_rep = 1-((wtplusm(100,4)+wtplusm(100,5))/wtminusm(100,2));
            
            rep(i,j,k) = protein_rep;
        end
    end
end
            
% 
% copy_array = [k_on,k_off,kx_s,b_e,b_ms,k_elon_prime,k_init_wt,k_init_rne,kx_wt_plus,kx_rne_plus,b_m,kx_wt_minus,kx_rne_minus,k_elon,k_nuc,b_nuc_wt,b_nuc_rne,k_off_nuc;
%     k_on_sd,k_off_sd,kxs_sd,be_sd,b_ms_sd,k_elon_prime_sd,k_init_wt_sd,k_init_rne_sd,kxw_sd,kxr_sd,b_m_sd,kxw_sd,kxr_sd,0,k_nuc_sd,b_nuc_wt_sd,b_nuc_rne_sd,k_off_nuc_sd];
% 
% copy_table = array2table(copy_array);copy_table.Properties.VariableNames = {'k_on' 'k_off' 'kx_s' 'b_e' 'b_ms' 'k_elon_prime' 'k_init_wt' 'k_init_rne_' 'kx_wt_plus' 'kx_rne_plus' 'b_m' 'kx_wt_minus' 'kx_rne_minus' 'k_elon' 'k_nuc' 'b_nuc_wt' 'b_nuc_rne' 'k_off_nuc'};
% 
% copy_table_file = strcat([folderTitle,'parameters.csv']);
% 
% writetable(copy_table,copy_table_file);

figure(1);
pcolor(k_on,Pt,reshape(rep(6,:,:),11,11)');
colormap jet
shading interp
colorbar;
xlabel('k_{on}','FontSize',14);
ylabel('P','FontSize',14);
title('{\beta_{e}} = 1.0x10^{-3}','FontSize',16);
caxis([0.0,0.8])
set(gca, 'FontName', 'Arial')
set(gca, 'XTick', [3E-5,3E-4,3.15E-3])
set(gca, 'XTickLabel', {'3E-5','3E-4','3E-3'})
set(gca, 'XScale', 'log')
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[435,953,396,346])
file1 = strcat([folderTitle,'k_on_P']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'k_on_P.fig']);
savefig(gcf,file1_fig)


