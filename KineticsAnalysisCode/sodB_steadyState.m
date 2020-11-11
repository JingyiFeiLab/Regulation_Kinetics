folderTitle = '/Users/reyer/Data/SingleCellEpi/MR243/Post8420/';

b_e = 2.23E-3; %std = 2.23E-3
be_sd = 0.0658;

a_s_wt = 0.263325; %a_s = 2.29E+03;
a_s_wt_sd = 0.067;
a_s_rne = 0.244975; %a_s = 2.29E+03; 
a_s_rne_sd = 0.0734;
b_s_wt = 0.00279305; 
b_s_wt_sd = 0.000233;
b_s_rne = 0.0019644; 
b_s_rne_sd = 0.000345;

b_m = 0.005404666667;
b_m_sd = 0.000314038745;

k_on = 1.38E-3;
k_on_sd = 0.00001707;
k_off = 9.57;
k_off_sd = 0.884;

k_nuc = 1.38E-3;
k_nuc_sd = 0;
k_off_nuc = 9.57;
k_off_nuc_sd = 0;

b_ms = 0.00560; % b_ms = 4.88E-03;
b_ms_sd = 0.0000645;
b_p = .00018;

b_nuc_wt = b_e + b_ms;
b_nuc_wt_sd = 0;
b_nuc_rne = b_ms;
b_nuc_rne_sd = 0;

kx_wt_minus =  28.707;
kx_wt_plus = 28.4625; % kx_wt = 0.0040191;
kxw_sd = 2.19;

kx_rne_minus = 10.355625;
kx_rne_plus = 9.105; %kx_rne = 0.0023639;
kxr_sd = 1.757;

kx_s = 0.364; %std = .15
kxs_sd = 0.148;

k_init_wt = 0.1317; %std = 310.41
k_init_wt_sd = 0.02577;
k_init_rne = 0.1067; %std = 1164.1
k_init_rne_sd = 0.0122;

k_elon = 0.0592;
k_elon_prime = 5.38E-03;
k_elon_prime_sd = .0001;

n_plasmids = 20;


close all
tv = linspace(0, 18000);
tvp = linspace(0, 50000);
tvm = linspace(0,1620);
%tvm = linspace(time(1), max(time)-time(3));

timex = [0,60,180,360,720,1080,1440,1800];

%Protein
%mRNA
%sRNA

molecules_rne_plus1= [537540.133,489861.7415,495303.3963,636534.158,1067081.43,1469663.453,2574256.314,2543087.016;
                      15372.71674,28146.29857,440707.2556,1226991.127,1666440.507,1320825.78,1559866.142,344620.4189;
                      593375.0032,455116.0876,568228.2777,511043.6822,508376.8531,461330.1373,395811.775,230814.0326];

molecules_rne_plus2= [536904.3621,461661.7087,495223.7801,678235.2065,1227608.695,1859059.735,2832320.07,3746419.851;
                      19898.01018,123286.8915,428240.4877,1243484.439,1724017.649,1690148.639,1684720.668,1434147.074;
                      461429.9867,478100.4001,485829.6041,455485.2018,518848.1603,573561.3758,465289.356,276181.5914];

molecules_wt_plus1 = [1130948.676,892617.3287,1074207.349,1375535.856,2229638.274,2729098.729,3567882.043,5582180.8;
                      10050.51839,305731.8887,803354.734,860573.0497,284796.5291,881811.9806,918839.9995,920764.8788;
                      430758.4504,399694.8762,385861.0661,425008.5325,117141.0959,317523.1091,392639.7222,311909.1325];                  

molecules_wt_plus2 = [1077656.5,1109629.183,1219469.886,1555835.48,2574231.765,3197925.634,4967996.936,5582180.8;
                      4809.786906,405355.402,570542.6576,739668.6922,788078.7262,859802.2183,782590.881,739851.4549;
                      231573.46,373192.0521,272348.6771,291664.8198,234679.4617,264677.7101,288519.334,210534.048]; 

molecules_wt_plus3 = [1130948.676,892617.3287,1074207.349,1375535.856,2229638.274,2729098.729,3567882.043,5582180.8;
                      10050.51839,305731.8887,803354.734,860573.0497,284796.5291,881811.9806,918839.9995,920764.8788;
                      430758.4504,399694.8762,385861.0661,425008.5325,117141.0959,317523.1091,392639.7222,311909.1325];                  

molecules_wt_plus4 = [1077656.5,1109629.183,1219469.886,1555835.48,2574231.765,3197925.634,4967996.936,5582180.8;
                      4809.786906,405355.402,570542.6576,739668.6922,788078.7262,859802.2183,782590.881,739851.4549;
                      231573.46,373192.0521,272348.6771,291664.8198,234679.4617,264677.7101,288519.334,210534.048];                  

                  
molecules_rne_minus1 = [433530.7916,431222.126,475938.4362,889377.6006,2074811.897,2709363.979,3824368.658,4063640.3139;
                        39417.71584,357559.7684,1405776.874,1801445.103,2190769.752,1816149.963,1737195.49,1669841.175;
                        0,0,0,0,0,0,0,0];
                    
molecules_rne_minus2 = [420793.839,398148.9589,447063.8845,712223.7411,1500282.596,2534672.953,3767130.147,4063640.3139;
                        29207.58194,371634.2204,1230121.393,1145491.732,1973455.387,1958699.105,2026503.723,2097355.703;
                        0,0,0,0,0,0,0,0]; 
                    
molecules_rne_minus3 = [129365.3871,62095.01169,58053.90591,270562.8565,1844490.394,3198326.792,4518674.053,4063640.3139;
                        25111.01036,22779.52205,402811.9358,1889523.498,3011722.555,3276609.047,2461194.629,2097355.703;
                        0,0,0,0,0,0,0,0];  
                    
molecules_rne_minus4 = [122189.2402,132940.12,189326.1,261839.7756,1093472.717,2072173.959,2607536.232,4063640.3139;
                        31736.66866,20709.34803,377412.5896,1678655.123,2537454.223,2470225.672,2053520.843,2097355.703;
                        0,0,0,0,0,0,0,0]; 
                    
molecules_wt_minus1 = [1330698.707,1415038.705,1658305.18,2926270.67,7015635.852,11537931.43,15246471.35,15619887.57;
                       24747.30501,529791.8018,1658133.555,2803290.984,2469521.543,1515026.528,3382199.233,1385334.766;
                       0,0,0,0,0,0,0,0];
                   
molecules_wt_minus2 = [361029.6549,193487.9339,608695.9576,2260052.762,6523354.828,9384318.365,11603718.83,14994865.29;
                       28758.37827,837864.1126,2009476.752,1939459.261,1849894.038,3720234.173,2781804.817,169693.2643;
                       0,0,0,0,0,0,0,0];  
                                          
                   
                   
moleculesrne_plus = [(9/4)*mean((([molecules_rne_plus1(3,1),molecules_rne_plus2(3,1)])+5.2879e+04)/7.3302e+03) 0 0 mean(((([molecules_rne_plus1(2,1),molecules_rne_plus2(2,1)])+1.3257e+05)/6.0368e+03)+1) 0 mean([molecules_rne_plus1(1,4),molecules_rne_plus2(1,4)])];
moleculeswt_plus = [(9/4)*mean((([molecules_wt_plus1(3,1),molecules_wt_plus2(3,1),molecules_wt_plus3(3,1),molecules_wt_plus4(3,1)])+5.2879e+04)/7.3302e+03) 0 0 mean(((([molecules_wt_plus1(2,1),molecules_wt_plus2(2,1),molecules_wt_plus3(2,1),molecules_wt_plus4(2,1)])+1.3257e+05)/6.0368e+03)+1)/20 0 mean([molecules_wt_plus1(1,4),molecules_wt_plus2(1,4),molecules_wt_plus3(1,4),molecules_wt_plus4(1,4)])];
moleculesrne_minus = [0 mean(((([molecules_rne_minus1(2,1),molecules_rne_minus2(2,1)])+1.3257e+05)/6.0368e+03)+1) mean([molecules_rne_minus1(1,4),molecules_rne_minus2(1,4)])];
moleculeswt_minus = [0 mean(((([molecules_wt_minus1(2,1),molecules_wt_minus2(2,1)])+1.3257e+05)/6.0368e+03)+1)/20 mean([molecules_wt_minus1(1,4),molecules_wt_minus2(1,4)])];

std_rne_plus = [(9/4)*std((([molecules_rne_plus1(3,1),molecules_rne_plus2(3,1)])+5.2879e+04)/7.3302e+03) 0 0  std(((([molecules_rne_plus1(2,1),molecules_rne_plus2(2,1)])+1.3257e+05)/6.0368e+03)+1) 0 std([molecules_rne_plus1(1,4),molecules_rne_plus2(1,4)])];
std_wt_plus = [(9/4)*std((([molecules_wt_plus1(3,1),molecules_wt_plus2(3,1),molecules_wt_plus3(3,1),molecules_wt_plus4(3,1)])+5.2879e+04)/7.3302e+03) 0 0 std(((([molecules_wt_plus1(2,1),molecules_wt_plus3(2,1),molecules_wt_plus4(2,1),molecules_wt_plus2(2,1)])+1.3257e+05)/6.0368e+03)+1) 0 std([molecules_wt_plus1(1,4),molecules_wt_plus2(1,4),molecules_wt_plus3(1,4),molecules_wt_plus4(1,4)])];
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