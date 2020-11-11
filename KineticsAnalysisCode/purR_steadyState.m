folderTitle = '/Users/reyer/Data/SingleCellEpi/MR188/Post8420/';

b_e = .569; %std = 2.23E-3
be_sd = 0.124;

a_s_wt = 0.3507; %a_s = 2.29E+03;
a_s_wt_sd = 0.02771858582;
a_s_rne = 0.36245; %a_s = 2.29E+03; 
a_s_rne_sd = 0.05064298767;
b_s_wt = 1.46E-03; 
b_s_wt_sd = 0.0002935200249;
b_s_rne = 0.85E-03; 
b_s_rne_sd = 0.0001349089028;

b_m = 0.01648033333;
b_m_sd = 0.0002649647775;

k_on = 3.33E-05;
k_on_sd = 0.000012374;
k_off = 5.69E-01;
k_off_sd = 0.197;

k_nuc = 0;
k_nuc_sd = 0;
k_off_nuc = 0;
k_off_nuc_sd = 0;

b_ms = 0.025945; % b_ms = 4.88E-03;
b_ms_sd = 0.007641;
b_p = .00018;

b_nuc_wt = b_e + b_ms;
b_nuc_wt_sd = 0.001594820993;
b_nuc_rne = b_ms;
b_nuc_rne_sd = 0;

kx_wt_minus =  1.5424;
kx_wt_plus = 1.646555556; % kx_wt = 0.0040191;
kxw_sd = 0.6590235201;

kx_rne_minus = 0.929;
kx_rne_plus = 1.01; %kx_rne = 0.0023639;
kxr_sd = 0.100464274;

kx_s = 0.8715; %std = .15
kxs_sd = 0.197;

k_init_wt = 0.0830215; %std = 310.41
k_init_wt_sd = 0.01727667989;
k_init_rne = 0.0899; %std = 1164.1
k_init_rne_sd = 0.01542453962;

k_elon = 0.0642;
k_elon_prime = 2.18E-02;
k_elon_prime_sd = .0001;

n_plasmids = 20;


close all
tv = linspace(0, 18000);
tvp = linspace(0, 108000);
tvm = linspace(0,1620);
%tvm = linspace(time(1), max(time)-time(3));

timex = [0,60,180,360,720,1080,1440,1800];

%Protein
%mRNA
%sRNA

molecules_rne_plus1= [26322.6146, 27390.38083, 21141.64156, 73754.83088, 79464.32493, 130074.8845, 170335.4906, 299159.4048;
                      24808.99509, 39363.23847, 168920.4282, 673637.5844, 603118.9613, 664133.2157, 642762.3578, 310908.97;
                      2088246.914, 2527482.1, 2412518.711, 2263463.165, 2414464.341, 2236072.449, 2299019.028, 923660.3895];

molecules_rne_plus2= [27710.95837, 25128.26139, 20590.76373, 38990.34986, 96988.16135, 147640.5684, 249792.2853, 201059.3225;
                      26004.23602, 35248.53293, 146375.2838, 615256.9138, 682525.7425, 642808.9317, 470719.4891, 649571.2048;
                      1692658.081, 2681364.609, 2204954.965, 2781083.839, 2326883.032, 2358047.508, 1364166.918, 1897987.52];
               
molecules_wt_plus1 = [20974.00014, 35006.82586, 24690.27744, 37263.60772, 81716.00854, 146370.0259, 96350.39233, 213774.9463;
                      21212.35091, 150312.015, 385766.2164, 342300.7812, 328246.2699, 356708.0122, 291218.9171, 486851.3321;
                      1241870.265, 1360925.379, 1440207.094, 1244487.329, 790049.0588, 953916.3595, 851013.6354, 708576.5192];
                  
molecules_wt_plus2 = [21010.87254, 23413.2199, 31874.47393, 43130.28607, 103514.2814, 168703.1014, 197693.4926, 172718.9147;
                      22667.51327, 254573.3561, 411830.7148, 376269.9698, 447103.5372, 333831.7481, 348654.2354, 644584.5322;
                      1690646.757, 1697123.941, 1401503.339, 1408826.706, 1184115.118, 690799.2367, 835454.4176, 749178.9447];
                  
molecules_rne_minus1 = [23991.3162,10432.40197,16223.15083,15940.80184,28717.68909,73112.88137,117449.6955,122362.238;
                        26701.82174,44417.53592,306766.8431,513891.2415,680667.4278,414183.6315,531892.4074,425906.8605;
                        0,0,0,0,0,0,0,0];
                    
molecules_rne_minus2 = [23991.3162,10432.40197,16223.15083,15940.80184,28717.68909,73112.88137,117449.6955,122362.238;
                        26701.82174,44417.53592,306766.8431,513891.2415,680667.4278,414183.6315,531892.4074,425906.8605;
                        0,0,0,0,0,0,0,0];
                    
                    
molecules_wt_minus1 = [18002.82139,52128.2963,29679.66145,83005.02853,92207.08011,159079.2027,177877.4467,284030.386;
                       24490.27815,75637.16402,561018.8064,481982.4034,656617.3933,362695.5841,415049.9444,442464.6317;
                       0,0,0,0,0,0,0,0];
                   
molecules_wt_minus2 = [16769.50244,29319.67688,26557.92336,66450.46066,104952.1526,178393.3847,232741.7802,283348.0441;
                       39563.15041,16118.45855,275596.0389,361665.8009,419283.7848,304101.002,342408.0657,362914.2187;
                       0,0,0,0,0,0,0,0]; 
 
                   
moleculesrne_plus = [mean((([molecules_rne_plus1(3,1),molecules_rne_plus2(3,1)])+5.2879e+04)/7.3302e+03) 0 0 mean(((([molecules_rne_plus1(2,1),molecules_rne_plus2(2,1)])+1.3257e+05)/6.0368e+03)+1) 0 mean([molecules_rne_plus1(1,4),molecules_rne_plus2(1,4)])];
moleculeswt_plus = [mean((([molecules_wt_plus1(3,1),molecules_wt_plus2(3,1)])+5.2879e+04)/7.3302e+03) 0 0 mean(((([molecules_wt_plus1(2,1),molecules_wt_plus2(2,1)])+1.3257e+05)/6.0368e+03)+1)/20 0 mean([molecules_wt_plus1(1,4),molecules_wt_plus2(1,4)])];
moleculesrne_minus = [0 mean(((([molecules_rne_minus1(2,1),molecules_rne_minus2(2,1)])+1.3257e+05)/6.0368e+03)+1) mean([molecules_rne_minus1(1,4),molecules_rne_minus2(1,4)])];
moleculeswt_minus = [0 mean(((([molecules_wt_minus1(2,1),molecules_wt_minus2(2,1)])+1.3257e+05)/6.0368e+03)+1)/20 mean([molecules_wt_minus1(1,4),molecules_wt_minus2(1,4)])];

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

rneplusp(:,1:6) = bestInit_fit(tvp,a_s_rne,b_s_rne,k_init_rne,k_elon,k_elon_prime,b_m,k_on,k_off,kx_rne_plus,kx_s*kx_rne_plus,b_ms,0,b_p,k_nuc,b_nuc_rne,n_plasmids,k_off_nuc,moleculesrne_plus);
rneplusm(:,1:6) = bestInit_fit(tv,a_s_rne,b_s_rne,k_init_rne,k_elon,k_elon_prime,b_m,k_on,k_off,kx_rne_plus,kx_s*kx_rne_plus,b_ms,0,b_p,k_nuc,b_nuc_rne,n_plasmids,k_off_nuc,moleculesrne_plus);

wtplusp(:,1:6) = bestInit_fit(tvp,a_s_wt,b_s_wt,k_init_wt,k_elon,k_elon_prime,b_m,k_on,k_off,kx_wt_plus,kx_s*kx_wt_plus,b_ms,b_e,b_p,k_nuc,b_nuc_wt,n_plasmids,k_off_nuc,moleculeswt_plus);
wtplusm(:,1:6) = bestInit_fit(tv,a_s_wt,b_s_wt,k_init_wt,k_elon,k_elon_prime,b_m,k_on,k_off,kx_wt_plus,kx_s*kx_wt_plus,b_ms,b_e,b_p,k_nuc,b_nuc_wt,n_plasmids,k_off_nuc,moleculeswt_plus);

rneminusp(:,1:3) = initiation_fit([k_init_rne,kx_rne_minus],tvp,k_elon,b_m,b_p,n_plasmids,moleculesrne_minus);
rneminusm(:,1:3) = initiation_fit([k_init_rne,kx_rne_minus],tv,k_elon,b_m,b_p,n_plasmids,moleculesrne_minus);

wtminusp(:,1:3) = initiation_fit([k_init_wt,kx_wt_minus],tvp,k_elon,b_m,b_p,n_plasmids,moleculeswt_minus);
wtminusm(:,1:3) = initiation_fit([k_init_wt,kx_wt_minus],tv,k_elon,b_m,b_p,n_plasmids,moleculeswt_minus);

figure(1)
plot(tv,(wtplusm(:,4)+wtplusm(:,5)),'LineWidth',1,'Color','k')
hold on
plot(tv,wtminusm(:,2),'LineWidth',1,'Color','k')

figure(2)
plot(tvp+360,wtplusp(:,6),'LineWidth',1,'Color','k')
hold on
plot(tvp+360,wtminusp(:,3),'LineWidth',1,'Color','k')

figure(3)
plot(tv,wtplusm(:,1)+wtplusm(:,3)+wtplusm(:,5))

protein_rep = 1-(wtplusp(100,6)/wtminusp(100,3))
mRNA_rep = 1-((wtplusm(100,4)+wtplusm(100,5))/wtminusm(100,2))