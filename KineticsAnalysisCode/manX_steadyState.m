folderTitle = '/Users/reyer/Data/SingleCellEpi/MR187/Fits8420/';

b_e = 2.23E-3; %std = 2.23E-3
be_sd = 0.00140;

a_s_wt = 0.3646; %a_s = 2.29E+03;
a_s_wt_sd = 0.02771858582;
a_s_rne = 0.37649; %a_s = 2.29E+03; 
a_s_rne_sd = 0.05064298767;
b_s_wt = 1.46E-03; 
b_s_wt_sd = 0.0002935200249;
b_s_rne = 0.85E-03; 
b_s_rne_sd = 0.0001349089028;

b_m = 3.25E-03;
b_m_sd = 0.0000833446659;

k_on = 3.40E-4;
k_on_sd = 0.0001186;
k_off = .211;
k_off_sd = 0.0629;

k_nuc = 3.40E-4;
k_nuc_sd = 0.0001186;
k_off_nuc = .211;
k_off_nuc_sd = 0.0629;

b_ms = 0.0037255; % b_ms = 4.88E-03;
b_ms_sd = 0.000304;
b_p = .00018;

b_nuc_wt = b_e + b_ms;
b_nuc_wt_sd = 0.001594820993;
b_nuc_rne = b_ms;
b_nuc_rne_sd = 0;

kx_wt_minus =  10.5897;
kx_wt_plus = 10.06; % kx_wt = 0.0040191;
kxw_sd = 0.9456;

kx_rne_minus = 5.4598;
kx_rne_plus = 5.4492; %kx_rne = 0.0023639;
kxr_sd = 0.287;

kx_s = 0.2288; %std = .15
kxs_sd = 0.121;

k_init_wt = 0.0582535; %std = 310.41
k_init_wt_sd = 0.006884137549;
k_init_rne = 0.0432; %std = 1164.1
k_init_rne_sd = 0.006404163056;

k_elon = 0.0641;
k_elon_prime = 2.18E-02;
k_elon_prime_sd = .0001;

n_plasmids = 20;


close all
tv = linspace(0, 18000);
tvp = linspace(0, 100000);
tvm = linspace(0,1620);
%tvm = linspace(time(1), max(time)-time(3));

timex = [0,60,180,360,720,1080,1440,1800];

%Protein
%mRNA
%sRNA

molecules_rne_plus1= [30740.06404,7611.392039,28438.54689,23048.9751,167085.4364,296431.9934,591670.9207,712535.0848;
                      22969.08289,13337.21322,162812.0361,435224.8164,905086.9004,1067163.977,1082526.701,1055862.227;
                      2085525.097,2072528.385,2038819.034,2209784.252,1869579.711,2355808.953,1676761.391,1123533.284];

molecules_rne_plus2= [25053.85663,27674.565,30010.12027,67783.35453,159716.7395,416788.5436,627153.735,323987.3492;
                      21402.5558,470096.1962,114914.2647,502208.8933,797456.0465,1006640.656,859554.7947,666903.2247;
                      2056380.255,2030448.259,2268022.237,2531172.333,2419731.03,2262055.679,1935417.938,831907.646];

molecules_wt_plus1 = [347513.2853,367107.5708,367229.4779,526899.8252,765383.2836,893353.6218,1194864.312,1611808.456;
                      126759.8746,174201.1888,861566.6491,829035.7939,1255372.267,1307704.871,1447041.427,1047774.86;
                      1846825.091,2026286.474,1517518.905,1095534.472,1302512.952,1455675.663,1423814.456,1391623.404];
                  
molecules_wt_plus2 = [435500.661,412698.7449,429219.6278,663807.8444,1035433.729,2144265.792,2486745.291,2684783.042;
                      55235.83625,123960.7343,1164299.635,1652964.543,1516875.402,1577627.268,1183959.582,1500131.281;
                      1371638.973,1479740.334,2300059.658,1877851.897,1420127.483,1090447.296,1096737.351,1324216.114];

molecules_wt_plus3 = [103121.4438,119740.8521,178169.0566,173472.7545,344063.0791,520752.0862,820034.7537,1184189.091;
                      18532.05393,42313.70869,650846.2755,1082916.612,1025421.802,1344767.177,1548645.684,1403606.941;
                      1142440.97,949293.0491,1536117.023,1213753.668,664036.2236,1182033.928,1183425.272,1304065.765];                  
                  
molecules_rne_minus1 = [154146.7479,145661.2402,245064.8148,258998.7581,437089.8109,402480.9372,1578219.724,1915866.24;
                        77546.75143,66964.48,214489.4754,876146.1855,1206895.125,555078.3933,1668292.438,1749240.381;
                        0,0,0,0,0,0,0,0];
                    
molecules_rne_minus2 = [117257.5384,128602.1801,134953.9028,195812.0289,463044.4978,794316.4453,1148494.574,1511621.315;
                        31169.97458,30772.43829,111422.0712,621586.7071,1236975.931,1393271.623,1025940.846,1018813.884;
                        0,0,0,0,0,0,0,0]; 

molecules_wt_minus1 = [294702.7318,291576.5226,346858.3621,448533.8108,1196744.873,1976856.739,3234456.669,3688022.195;
                       92213.41176,102255.6987,1002772.753,1367403.667,1694510.768,2176574.682,1479485.13,1785058.464;
                       0,0,0,0,0,0,0,0];
                   
molecules_wt_minus2 = [258457.3074,373669.483,313491.9673,512902.479,1151679.312,2036671.374,3146130.093,3466090.969;
                       69374.65311,161183.5205,912699.742,1928275.948,2058956.949,1544598.234,2087732.002,1633831.709;
                       0,0,0,0,0,0,0,0]; 

                   

moleculesrne_plus = [mean((([molecules_rne_plus1(3,1),molecules_rne_plus2(3,1)])+5.2879e+04)/7.3302e+03) 0 0 mean(((([molecules_rne_plus1(2,1),molecules_rne_plus2(2,1)])+1.3257e+05)/6.0368e+03)+1) 0 mean([molecules_rne_plus1(1,4),molecules_rne_plus2(1,4)])];
moleculeswt_plus = [mean((([molecules_wt_plus1(3,1),molecules_wt_plus2(3,1),molecules_wt_plus3(3,1)])+5.2879e+04)/7.3302e+03) 0 0 mean(((([molecules_wt_plus1(2,1),molecules_wt_plus2(2,1),molecules_wt_plus3(2,1)])+1.3257e+05)/6.0368e+03)+1)/20 0 mean([molecules_wt_plus1(1,4),molecules_wt_plus2(1,4),molecules_wt_plus3(1,4)])];
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