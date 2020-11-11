folderTitle = '/Users/reyer/Data/SingleCellEpi/MR156/UncoupledPost8520/';

b_e = 1.87E-02; %std = 2.23E-3
be_sd = 0.002823860721;

a_m_plus_wt = 1.60515;% std =589.18; 
amp_wt_sd = 0.004797568829;
a_m_plus_rne = 1.419; %std = 796.64
amp_rne_sd = 0.005597618541;

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

k_on = 7.35E-04;
k_on_sd = 0.0006410343481;
k_off = 5.21E-01;
k_off_sd = 0.2378706161;

kx_wt_minus = 13.7033;
kx_wt_plus = 12.24; % kx_wt = 0.0040191;
kxw_sd = 1.698;
% kx_wt_plus=mode(kx_fit); %kx_wt_plus=0.0033;
kx_rne_minus = 5.4504;
kx_rne_plus = 5.184; %kx_rne = 0.0023639;
kxr_sd = .437;

kx_s = 0.25875; %std = .15
kxs_sd = 0.2114211831;

b_ms_rne = 0.336386; % b_ms = 4.88E-03;
b_ms_wt = 0.336386; % b_ms = 4.88E-03;
b_ms_sd = 0.02918078865;
b_p = .00018;

a_m_minus_wt = 1.905; %std = 310.41
amw_sd = 0.305;
a_m_minus_rne = 1.431; %std = 1164.1
amm_sd = 0.031;

% b_e = 2.82E-03; %std = 2.23E-3
% be_sd = 0.0003917371568;
% a_m_plus_wt = 5.2318;% std =589.18; 
% amp_wt_sd = 0.7506645589;
% a_m_plus_rne = 5.8198; %std = 796.64
% amp_rne_sd = 0.2046367025;
% a_s_wt = 0.29845; %a_s = 2.29E+03;
% a_s_wt_sd = 0.09369164851;
% a_s_rne = 0.4121; %a_s = 2.29E+03; 
% a_s_rne_sd = 0.04369919908;
% b_s_wt = 1.44E-03; 
% b_s_wt_sd = 0.00005939696962;
% b_s_rne = 0.000070903679; 
% b_s_rne_sd = 0.0001447447581;
% b_m = 0.003279;
% b_m_sd = 0.00007212489168;
% k_on = 1.14E-04;
% k_on_sd = 0.0001287245469;
% k_off = 5.49E-02;
% k_off_sd = 0.01347745525;
% kx_wt_minus = 2.5196;
% kx_wt_plus = 1.251; % kx_wt = 0.0040191;
% kxw_sd = 0.1046518036;
% % kx_wt_plus=mode(kx_fit); %kx_wt_plus=0.0033;
% kx_rne_minus = 2.1407;
% kx_rne_plus = 1.66; %kx_rne = 0.0023639;
% kxr_sd = 0.03889087297;
% kx_s = 0.3885; %std = .15
% kxs_sd = 0.1463711037;
% b_ms_rne = 0.004015; % b_ms = 4.88E-03;
% b_ms_wt = 0.004015; % b_ms = 4.88E-03;
% b_ms_sd = 0.0001414213562;
% b_p = .00018;
% a_m_minus_wt = 5.5318; %std = 310.41
% amw_sd = 0.7506645589;
% a_m_minus_rne = 5.8198; %std = 1164.1
% amm_sd = 0.2046367025;

close all
tv = linspace(0, 1800);
tvp = linspace(0, 1440);
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
                   

                   
moleculesrne_plus = [mean((([molecules_rne_plus1(3,1),molecules_rne_plus2(3,1),molecules_rne_plus3(3,1),molecules_rne_plus4(3,1)])+5.2879e+04)/7.3302e+03) mean(((([molecules_rne_plus1(2,1),molecules_rne_plus2(2,1),molecules_rne_plus3(2,1),molecules_rne_plus4(2,1)])+1.3257e+05)/6.0368e+03)+1) 0 mean([molecules_rne_plus1(1,4),molecules_rne_plus2(1,4),molecules_rne_plus3(1,4),molecules_rne_plus4(1,4)])];
moleculeswt_plus = [mean((([molecules_wt_plus1(3,1),molecules_wt_plus2(3,1)])+5.2879e+04)/7.3302e+03) mean(((([molecules_wt_plus1(2,1),molecules_wt_plus2(2,1)])+1.3257e+05)/6.0368e+03)+1) 0 mean([molecules_wt_plus1(1,4),molecules_wt_plus2(1,4)])];
moleculesrne_minus = [0 mean(((([molecules_rne_minus1(2,1),molecules_rne_minus2(2,1)])+1.3257e+05)/6.0368e+03)+1) 0 mean([molecules_rne_minus1(1,4),molecules_rne_minus2(1,4)])];
moleculeswt_minus = [0 mean(((([molecules_wt_minus1(2,1),molecules_wt_minus2(2,1)])+1.3257e+05)/6.0368e+03)+1) 0 mean([molecules_wt_minus1(1,4),molecules_wt_minus2(1,4)])];

std_rne_plus = [std((([molecules_rne_plus1(3,1),molecules_rne_plus2(3,1)])+5.2879e+04)/7.3302e+03) std(((([molecules_rne_plus1(2,1),molecules_rne_plus2(2,1)])+1.3257e+05)/6.0368e+03)+1) 0 std([molecules_rne_plus1(1,4),molecules_rne_plus2(1,4)])];
std_wt_plus = [std((([molecules_wt_plus1(3,1),molecules_wt_plus2(3,1)])+5.2879e+04)/7.3302e+03) std(((([molecules_wt_plus1(2,1),molecules_wt_plus2(2,1)])+1.3257e+05)/6.0368e+03)+1) 0 std([molecules_wt_plus1(1,4),molecules_wt_plus2(1,4)])];
std_rne_minus = [0 std(((([molecules_rne_minus1(2,1),molecules_rne_minus2(2,1)])+1.3257e+05)/6.0368e+03)+1) 0 std([molecules_rne_minus1(1,4),molecules_rne_minus2(1,4)])];
std_wt_minus = [0 std(((([molecules_wt_minus1(2,1),molecules_wt_minus2(2,1)])+1.3257e+05)/6.0368e+03)+1) 0 std([molecules_wt_minus1(1,4),molecules_wt_minus2(1,4)])];

rneplusp = zeros(100,80);
rneplusm = zeros(100,80);
rnepluss = zeros(100,80);

wtplusp = zeros(100,80);
wtplusm = zeros(100,80);
wtpluss = zeros(100,80);

rneminusp = zeros(100,80);
rneminusm = zeros(100,80);
rneminuss = zeros(100,80);

wtminusp = zeros(100,80);
wtminusm = zeros(100,80);
wtminuss = zeros(100,80);

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

rneplusp(:,1:4) = roughPlussRNA_fit(tvp,0,a_m_plus_rne,a_s_rne,b_s_rne,b_m,k_on,k_off,kx_rne_plus,kx_s*kx_rne_plus,b_ms_rne,b_p,moleculesrne_plus);
rneplusm(:,1:4) = roughPlussRNA_fit(tv,0,a_m_plus_rne,a_s_rne,b_s_rne,b_m,k_on,k_off,kx_rne_plus,kx_s*kx_rne_plus,b_ms_rne,b_p,moleculesrne_plus);

wtplusp(:,1:4) = roughPlussRNA_fit(tvp,b_e,a_m_plus_wt,a_s_wt,b_s_wt,b_m,k_on,k_off,kx_wt_plus,kx_s*kx_wt_plus,b_ms_wt,b_p,moleculeswt_plus);
wtplusm(:,1:4) = roughPlussRNA_fit(tv,b_e,a_m_plus_wt,a_s_wt,b_s_wt,b_m,k_on,k_off,kx_wt_plus,kx_s*kx_wt_plus,b_ms_wt,b_p,moleculeswt_plus);

rneminusp(:,1:4) = nosRNA_fit([a_m_minus_rne,kx_rne_minus],tvp,0,0,0,0,0,b_m,0,0,b_p,moleculesrne_minus);
rneminusm(:,1:4) = nosRNA_fit([a_m_minus_rne,kx_rne_minus],tv,0,0,0,0,0,b_m,0,0,b_p,moleculesrne_minus);

wtminusp(:,1:4) = nosRNA_fit([a_m_minus_wt,kx_wt_minus],tvp,0,0,0,0,0,b_m,0,0,b_p,moleculeswt_minus);
wtminusm(:,1:4) = nosRNA_fit([a_m_minus_wt,kx_wt_minus],tv,0,0,0,0,0,b_m,0,0,b_p,moleculeswt_minus);

i_rand = 5:8;
i_rand_check = [];
for i = 1:19
    
    a1 = roughPlussRNA_fit(tvp,0,normrnd(a_m_plus_rne,amp_rne_sd),normrnd(a_s_rne,a_s_rne_sd),normrnd(b_s_rne,b_s_rne_sd),normrnd(b_m,b_m_sd),normrnd(k_on,k_on_sd),normrnd(k_off,k_off_sd),normrnd(kx_rne_plus,kxr_sd),normrnd(kx_s,kxs_sd)*kx_rne_plus,normrnd(b_ms_rne,b_ms_sd),b_p,normrnd(moleculesrne_plus,std_rne_plus));
    if length(a1(:,1)) == 100
        rneplusp(:,i_rand) = a1;
    else
        i_rand_check = [i_rand_check i_rand];
    end
    
    a2 = roughPlussRNA_fit(tv,0,normrnd(a_m_plus_rne,amp_rne_sd),normrnd(a_s_rne,a_s_rne_sd),normrnd(b_s_rne,b_s_rne_sd),normrnd(b_m,b_m_sd),normrnd(k_on,k_on_sd),normrnd(k_off,k_off_sd),normrnd(kx_rne_plus,kxr_sd),normrnd(kx_s,kxs_sd)*kx_rne_plus,normrnd(b_ms_rne,b_ms_sd),b_p,normrnd(moleculesrne_plus,std_rne_plus));
    if length(a2(:,1)) == 100
        rneplusm(:,i_rand) = a2;
    else
        i_rand_check = [i_rand_check i_rand];
    end
    
    a3 = roughPlussRNA_fit(tvp,normrnd(b_e,be_sd),normrnd(a_m_plus_wt,amp_wt_sd),normrnd(a_s_wt,a_s_wt_sd),normrnd(b_s_wt,b_s_wt_sd),normrnd(b_m,b_m_sd),normrnd(k_on,k_on_sd),normrnd(k_off,k_off_sd),normrnd(kx_wt_plus,kxw_sd),normrnd(kx_s,kxs_sd)*kx_wt_plus,normrnd(b_ms_wt,b_ms_sd),b_p,normrnd(moleculeswt_plus,std_wt_plus));
    if length(a3(:,1)) == 100
        wtplusp(:,i_rand) = a3;
    else
        i_rand_check = [i_rand_check i_rand];
    end
    
    a4 = roughPlussRNA_fit(tv,normrnd(b_e,be_sd),normrnd(a_m_plus_wt,amp_wt_sd),normrnd(a_s_wt,a_s_wt_sd),normrnd(b_s_wt,b_s_wt_sd),normrnd(b_m,b_m_sd),normrnd(k_on,k_on_sd),normrnd(k_off,k_off_sd),normrnd(kx_wt_plus,kxw_sd),normrnd(kx_s,kxs_sd)*kx_wt_plus,normrnd(b_ms_wt,b_ms_sd),b_p,normrnd(moleculeswt_plus,std_wt_plus));
    if length(a4(:,1)) == 100
        wtplusm(:,i_rand) = a4;
    else
        i_rand_check = [i_rand_check i_rand];
    end
   
    a5 = nosRNA_fit([normrnd(a_m_minus_rne,amm_sd),normrnd(kx_rne_minus,kxr_sd)],tvp,0,0,0,0,0,normrnd(b_m,b_m_sd),0,0,b_p,normrnd(moleculesrne_minus,std_rne_minus));
    if length(a5(:,1)) == 100
        rneminusp(:,i_rand) = a5;
    else
        i_rand_check = [i_rand_check i_rand];
    end
    
    a6 = nosRNA_fit([normrnd(a_m_minus_rne,amm_sd),normrnd(kx_rne_minus,kxr_sd)],tv,0,0,0,0,0,normrnd(b_m,b_m_sd),0,0,b_p,normrnd(moleculesrne_minus,std_rne_minus));
    if length(a6(:,1)) == 100
        rneminusm(:,i_rand) = a6;
    else
        i_rand_check = [i_rand_check i_rand];
    end
    
    a7 = nosRNA_fit([normrnd(a_m_minus_wt,amw_sd),normrnd(kx_wt_minus,kxw_sd)],tvp,0,0,0,0,0,normrnd(b_m,b_m_sd),0,0,b_p,normrnd(moleculeswt_minus,std_wt_minus));
    if length(a7(:,1)) == 100
        wtminusp(:,i_rand) = a7;
    else
        i_rand_check = [i_rand_check i_rand];
    end
    
    a8 = nosRNA_fit([normrnd(a_m_minus_wt,amw_sd),normrnd(kx_wt_minus,kxw_sd)],tv,0,0,0,0,0,normrnd(b_m,b_m_sd),0,0,b_p,normrnd(moleculeswt_minus,std_wt_minus));
    if length(a8(:,1)) == 100
        wtminusm(:,i_rand) = a8;
    else
        i_rand_check = [i_rand_check i_rand];
    end
    
    i_rand = i_rand + 4;
    
end

rne_minus_protein_error = zeros(100,20);
rne_plus_protein_error = zeros(100,20);
rne_minus_mRNA_error = zeros(100,20);
rne_plus_mRNA_error = zeros(100,20);
rne_minus_sRNA_error = zeros(100,20);
rne_plus_sRNA_error = zeros(100,20);
wt_minus_protein_error = zeros(100,20);
wt_plus_protein_error = zeros(100,20);
wt_minus_mRNA_error = zeros(100,20);
wt_plus_mRNA_error = zeros(100,20);
wt_minus_sRNA_error = zeros(100,20);
wt_plus_sRNA_error = zeros(100,20);

for i = 1:100
    ij = 1;
    for j = 2:4:80
        rne_minus_mRNA_error(i,ij) = rneminusm(i,j)+rneminusm(i,j+1);
        rne_plus_mRNA_error(i,ij) = rneplusm(i,j) + rneplusm(i,j+1);
        wt_minus_mRNA_error(i,ij) = wtminusm(i,j) + wtminusm(i,j+1);
        wt_plus_mRNA_error(i,ij) = wtplusm(i,j) + wtplusm(i,j+1);
        ij = ij + 1;
        
    end

    ik = 1;
    for k = 4:4:80
        rne_minus_protein_error(i,k) = rneminusp(i,k);
        rne_plus_protein_error(i,k) = rneplusp(i,k);
        wt_minus_protein_error(i,k) = wtminusp(i,k);
        wt_plus_protein_error(i,k) = wtplusp(i,k);
        ik = ik + 1;
    end
    
    il = 1;
    for j = 1:4:80
        rne_minus_sRNA_error(i,il) = rneminusm(i,j)+rneminusm(i,j+2);
        rne_plus_sRNA_error(i,il) = rneplusm(i,j) + rneplusm(i,j+2);
        wt_minus_sRNA_error(i,il) = wtminusm(i,j) + wtminusm(i,j+2);
        wt_plus_sRNA_error(i,il) = wtplusm(i,j) + wtplusm(i,j+2);
        il = il + 1;
        
    end
    
end

for i = 1:100
    
    percentiles = [];
    nonOutliers = [];
    percentiles = prctile(rne_plus_protein_error(i,:),[5,95]);
    outlierIndexes = rne_plus_protein_error(i,:) < percentiles(1) | rne_plus_protein_error(i,:) > percentiles(2);
    nonOutliers = rne_plus_protein_error(i,~outlierIndexes);
    
    rneplus_pspread(i) = std(nonOutliers);
    
    percentiles = [];
    nonOutliers = [];
    percentiles = prctile(rne_plus_mRNA_error(i,:),[5,95]);
    outlierIndexes = rne_plus_mRNA_error(i,:) < percentiles(1) | rne_plus_mRNA_error(i,:) > percentiles(2);
    nonOutliers = rne_plus_mRNA_error(i,~outlierIndexes);
    
    rneplus_mspread(i) = std(nonOutliers)/1;
    
    percentiles = [];
    nonOutliers = [];
    percentiles = prctile(rne_plus_sRNA_error(i,:),[5,95]);
    outlierIndexes = rne_plus_sRNA_error(i,:) < percentiles(1) | rne_plus_sRNA_error(i,:) > percentiles(2);
    nonOutliers = rne_plus_sRNA_error(i,~outlierIndexes);
    
    rneplus_sspread(i) = std(nonOutliers)/2;
    
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
    
    percentiles = [];
    nonOutliers = [];
    percentiles = prctile(rne_minus_protein_error(i,:),[5,95]);
    outlierIndexes = rne_minus_protein_error(i,:) < percentiles(1) | rne_minus_protein_error(i,:) > percentiles(2);
    nonOutliers = rne_minus_protein_error(i,~outlierIndexes);
    
    rneminus_pspread(i) = std(nonOutliers);
    
    percentiles = [];
    nonOutliers = [];
    percentiles = prctile(rne_minus_mRNA_error(i,:),[5,95]);
    outlierIndexes = rne_minus_mRNA_error(i,:) < percentiles(1) | rne_minus_mRNA_error(i,:) > percentiles(2);
    nonOutliers = rne_minus_mRNA_error(i,~outlierIndexes);
    
    rneminus_mspread(i) = std(nonOutliers);
    
    
    rneminus_sspread(i) = 0;
    
    percentiles = [];
    nonOutliers = [];
    percentiles = prctile(wt_minus_protein_error(i,:),[5,95]);
    outlierIndexes = rne_minus_protein_error(i,:) < percentiles(1) | wt_minus_protein_error(i,:) > percentiles(2);
    nonOutliers = wt_minus_protein_error(i,~outlierIndexes);
    
    wtminus_pspread(i) = std(nonOutliers);
    
    percentiles = [];
    nonOutliers = [];
    percentiles = prctile(wt_minus_mRNA_error(i,:),[5,95]);
    outlierIndexes = wt_minus_mRNA_error(i,:) < percentiles(1) | wt_minus_mRNA_error(i,:) > percentiles(2);
    nonOutliers = wt_minus_mRNA_error(i,~outlierIndexes);
    
    wtminus_mspread(i) = std(nonOutliers);
    
    
    wtminus_sspread(i) = 0;
end  

molecules_rne_plus = zeros(6,8);
molecules_rne_minus = zeros(6,8);
molecules_wt_plus = zeros(6,8);
molecules_wt_minus = zeros(6,8);
for i = 1:8
    molecules_rne_plus(1,i) = mean([molecules_rne_plus1(1,i),molecules_rne_plus2(1,i),molecules_rne_plus3(1,i),molecules_rne_plus4(1,i)]);
    molecules_rne_minus(1,i) = mean([molecules_rne_minus1(1,i),molecules_rne_minus2(1,i)]); 
    molecules_wt_plus(1,i) = mean([molecules_wt_plus1(1,i),molecules_wt_plus2(1,i)]);
    molecules_wt_minus(1,i) = mean([molecules_wt_minus1(1,i),molecules_wt_minus2(1,i)]);
    
    molecules_rne_plus(2,i) = std([molecules_rne_plus1(1,i),molecules_rne_plus2(1,i),molecules_rne_plus3(1,i),molecules_rne_plus4(1,i)])/1.414;
    molecules_rne_minus(2,i) = std([molecules_rne_minus1(1,i),molecules_rne_minus2(1,i)])/1.414; 
    molecules_wt_plus(2,i) = std([molecules_wt_plus1(1,i),molecules_wt_plus2(1,i)])/1.414;
    molecules_wt_minus(2,i) = std([molecules_wt_minus1(1,i),molecules_wt_minus2(1,i)])/1.723;
    
    molecules_rne_plus(3,i) = mean((([molecules_rne_plus1(2,i),molecules_rne_plus2(2,i),molecules_rne_plus3(2,i),molecules_rne_plus4(2,i)])+1.3257e+05)/6.0368e+03)+1;
    molecules_rne_minus(3,i) = mean((([molecules_rne_minus1(2,i),molecules_rne_minus2(2,i)])+1.3257e+05)/6.0368e+03)+1; 
    molecules_wt_plus(3,i) = mean((([molecules_wt_plus1(2,i),molecules_wt_plus2(2,i)])+1.3257e+05)/6.0368e+03)+1;
    molecules_wt_minus(3,i) = mean((([molecules_wt_minus1(2,i),molecules_wt_minus2(2,i)])+1.3257e+05)/6.0368e+03)+1;
    
    molecules_rne_plus(4,i) = std(((([molecules_rne_plus1(2,i),molecules_rne_plus2(2,i),molecules_rne_plus3(2,i),molecules_rne_plus4(2,i)])+1.3257e+05)/6.0368e+03)+1)/1.414;
    molecules_rne_minus(4,i) = std(((([molecules_rne_minus1(2,i),molecules_rne_minus2(2,i)])+1.3257e+05)/6.0368e+03)+1)/1.414; 
    molecules_wt_plus(4,i) = std(((([molecules_wt_plus1(2,i),molecules_wt_plus2(2,i)])+1.3257e+05)/6.0368e+03)+1)/1.414;
    molecules_wt_minus(4,i) = std(((([molecules_wt_minus1(2,i),molecules_wt_minus2(2,i)])+1.3257e+05)/6.0368e+03)+1)/1.723;
    
    molecules_rne_plus(5,i) = mean((([molecules_rne_plus1(3,i),molecules_rne_plus2(3,i),molecules_rne_plus3(3,i),molecules_rne_plus4(3,i)])+5.2879e+04)/7.3302e+03);
    molecules_rne_minus(5,i) = mean((([molecules_rne_minus1(3,i),molecules_rne_minus2(3,i)])+5.2879e+04)/7.3302e+033); 
    molecules_wt_plus(5,i) = mean((([molecules_wt_plus1(3,i),molecules_wt_plus2(3,i)])+5.2879e+04)/7.3302e+03);
    molecules_wt_minus(5,i) = mean((([molecules_wt_minus1(3,i),molecules_wt_minus2(3,i)])+5.2879e+04)/7.3302e+03);
    
    molecules_rne_plus(6,i) = std((([molecules_rne_plus1(3,i),molecules_rne_plus2(3,i),molecules_rne_plus3(3,i),molecules_rne_plus4(3,i)])+5.2879e+04)/7.3302e+03)/1.414;
    molecules_rne_minus(6,i) = std((([molecules_rne_minus1(3,i),molecules_rne_minus2(3,i)])+5.2879e+04)/7.3302e+03)/1.414; 
    molecules_wt_plus(6,i) = std((([molecules_wt_plus1(3,i),molecules_wt_plus2(3,i)])+5.2879e+04)/7.3302e+03)/1.414;
    molecules_wt_minus(6,i) = std((([molecules_wt_minus1(3,i),molecules_wt_minus2(3,i)])+5.2879e+04)/7.3302e+03)/1.72;
    
    
end


figure(1)

e1 = errorbar(timex,molecules_rne_plus(3,:),molecules_rne_plus(4,:),'o');
e1.MarkerFaceColor = [0,.3,0];
e1.MarkerSize = 5;
e1.Color = [0,.3,0];
e1.LineWidth = 2;
hold on
e2 = errorbar(timex,molecules_rne_minus(3,:),molecules_rne_minus(4,:),'d');
e2.MarkerFaceColor = [0,.65,.1];
e2.MarkerSize = 5;
e2.Color = [0,.65,.1];
e2.LineWidth = 2;
hold on
shadedErrorBar(tv,rneplusm(:,2)+rneplusm(:,3),rneplus_mspread,'lineProps', {'Color',[0,.40,0]})
hold on
plot(tv,rneplusm(:,2)+rneplusm(:,3),'LineWidth',1,'Color','k')
hold on
shadedErrorBar(tv,rneminusm(:,2)+rneminusm(:,3),rneminus_mspread,'lineProps', {'Color',[0,.95,0]})
hold on
plot(tv,rneminusm(:,2)+rneminusm(:,3),'LineWidth',1,'Color','k')
hold on
e1 = errorbar(timex,molecules_rne_plus(3,:),molecules_rne_plus(4,:),'o');
e1.MarkerFaceColor = [0,.3,0];
e1.MarkerSize = 5;
e1.Color = [0,.3,0];
e1.LineWidth = 2;
hold on
e2 = errorbar(timex,molecules_rne_minus(3,:),molecules_rne_minus(4,:),'d');
e2.MarkerFaceColor = [0,.65,.1];
e2.MarkerSize = 5;
e2.Color = [0,.65,.1];
e2.LineWidth = 2;
lgd = legend('WT SgrS','{\Delta SgrS}');
lgd.FontSize = 10;
lgd.Location = 'northwest';
title('rne701{\it ptsG-sfGFP}','FontSize',16,'Color',[0,.40,0],'FontName','Arial')
xlabel('Time after Induction (min)','FontSize',14,'FontName','Arial')
ylabel('Copy Number ','FontSize',14,'FontName','Arial')
set(gca,'YLim',[-100 3.5E6/3025.7])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10)
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[626,281,248,201])
file1 = strcat([folderTitle,'rneptsGmRNA_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'rneptsGmRNA_fit.fig']);
savefig(gcf,file1_fig)

figure(2)
e3 = errorbar(timex,molecules_wt_plus(3,:),molecules_wt_plus(4,:),'s');
e3.MarkerFaceColor = [0,.3,0];
e3.MarkerSize = 5;
e3.Color = [0,.3,0];
e3.LineWidth = 2;
hold on
e4 = errorbar(timex,molecules_wt_minus(3,:),molecules_wt_minus(4,:),'*');
e4.MarkerFaceColor = [0,.65,.1];
e4.MarkerSize = 5;
e4.Color = [0,.65,.1];
e4.LineWidth = 2;
hold on
shadedErrorBar(tv,(wtplusm(:,2)+wtplusm(:,3)),wtplus_mspread,'lineProps', {'Color',[0,.40,0]})
hold on
plot(tv,(wtplusm(:,2)+wtplusm(:,3)),'LineWidth',1,'Color','k')
hold on
shadedErrorBar(tv,wtminusm(:,2)+wtminusm(:,3),wtminus_mspread,'lineProps', {'Color',[0,.95,0]})
hold on
plot(tv,wtminusm(:,2)+wtminusm(:,3),'LineWidth',1,'Color','k')
hold on
e3 = errorbar(timex,molecules_wt_plus(3,:),molecules_wt_plus(4,:),'s');
e3.MarkerFaceColor = [0,.3,0];
e3.MarkerSize = 5;
e3.Color = [0,.3,0];
e3.LineWidth = 2;
hold on
e4 = errorbar(timex,molecules_wt_minus(3,:),molecules_wt_minus(4,:),'*');
e4.MarkerFaceColor = [0,.65,.1];
e4.MarkerSize = 5;
e4.Color = [0,.65,.1];
e4.LineWidth = 2;
lgd = legend('WT SgrS','{\Delta SgrS}');
lgd.FontSize = 10;
lgd.Location = 'northwest';
title('WT{\it ptsG-sfGFP}','FontSize',16,'Color',[0,.40,0],'FontName','Arial')
xlabel('Time after Induction (min)','FontSize',14,'FontName','Arial')
ylabel('Copy Number','FontSize',14,'FontName','Arial')
set(gca,'YLim',[-100 3.5E6/3025.7])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10)
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


figure(3)
e5 = errorbar(timex,molecules_rne_plus(1,:),molecules_rne_plus(2,:),'o');
e5.MarkerFaceColor = [.05,.05,.53];
e5.MarkerSize = 5;
e5.Color = [.05,.05,.53];
e5.LineWidth = 2;
hold on
e6 = errorbar(timex,molecules_rne_minus(1,:),molecules_rne_minus(2,:),'d');
e6.MarkerFaceColor = [0,.6,.6];
e6.MarkerSize = 5;
e6.Color = [0,.6,.6];
e6.LineWidth = 2;
hold on
shadedErrorBar(tvp+360,rneplusp(:,4),rneplus_pspread,'lineProps', {'Color',[.05,.05,.38]})
hold on
plot(tvp+360,rneplusp(:,4),'LineWidth',1,'Color','k')
hold on
shadedErrorBar(tvp+360,rneminusp(:,4),rneminus_pspread,'lineProps', {'Color',[0,.95,.95]})
hold on
plot(tvp+360,rneminusp(:,4),'LineWidth',1,'Color','k')
hold on
e5 = errorbar(timex,molecules_rne_plus(1,:),molecules_rne_plus(2,:),'o');
e5.MarkerFaceColor = [.05,.05,.53];
e5.MarkerSize = 5;
e5.Color = [.05,.05,.53];
e5.LineWidth = 2;
hold on
e6 = errorbar(timex,molecules_rne_minus(1,:),molecules_rne_minus(2,:),'d');
e6.MarkerFaceColor = [0,.6,.6];
e6.MarkerSize = 5;
e6.Color = [0,.6,.6];
e6.LineWidth = 2;
lgd = legend('WT SgrS','{\Delta SgrS}');
lgd.FontSize = 10;
lgd.Location = 'northwest';
title('{\it rne}701 sfGFP','FontSize',16,'Color','b','FontName','Arial')
xlabel('Time after Induction (min)','FontSize',14,'FontName','Arial')
ylabel('Fluorescence (A.U)','FontSize',14,'FontName','Arial')
set(gca,'YLim',[-10000 6E6])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10)
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[626,281,248,201])
file1 = strcat([folderTitle,'rneptsGProtein_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'rneptsGProtein_fit.fig']);
savefig(gcf,file1_fig)

figure(4)
e7 = errorbar(timex,molecules_wt_plus(1,:),molecules_wt_plus(2,:),'s');
e7.MarkerFaceColor = [.05,.05,.53];
e7.MarkerSize = 5;
e7.Color = [.05,.05,.53];
e7.LineWidth = 2;
hold on
e8 = errorbar(timex,molecules_wt_minus(1,:),molecules_wt_minus(2,:),'*');
e8.MarkerFaceColor = [0,.6,.6];
e8.MarkerSize = 5;
e8.Color = [0,.6,.6];
e8.LineWidth = 2;
hold on
shadedErrorBar(tvp+360,wtplusp(:,4),wtplus_pspread,'lineProps', {'Color',[.05,.05,.38]})
hold on
plot(tvp+360,wtplusp(:,4),'LineWidth',1,'Color','k')
hold on
shadedErrorBar(tvp+360,wtminusp(:,4),wtminus_pspread,'lineProps', {'Color',[0,.95,.95]})
hold on
plot(tvp+360,wtminusp(:,4),'LineWidth',1,'Color','k')
hold on
e7 = errorbar(timex,molecules_wt_plus(1,:),molecules_wt_plus(2,:),'s');
e7.MarkerFaceColor = [.05,.05,.53];
e7.MarkerSize = 5;
e7.Color = [.05,.05,.53];
e7.LineWidth = 2;
hold on
e8 = errorbar(timex,molecules_wt_minus(1,:),molecules_wt_minus(2,:),'*');
e8.MarkerFaceColor = [0,.6,.6];
e8.MarkerSize = 5;
e8.Color = [0,.6,.6];
e8.LineWidth = 2;
lgd = legend('WT SgrS','{\Delta SgrS}');
lgd.FontSize = 10;
lgd.Location = 'northwest';
title('WT sfGFP','FontSize',16,'FontName','Arial','Color','b')
xlabel('Time after Induction (min)','FontSize',14,'FontName','Arial')
ylabel('Fluorescence (A.U)','FontSize',14,'FontName','Arial')
set(gca,'YLim',[-10000 8E6])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10)
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

figure(5)

e9 = errorbar(timex,molecules_rne_plus(5,:),molecules_rne_plus(6,:),'o');
e9.MarkerFaceColor = [1,.75,.79];
e9.MarkerSize = 5;
e9.Color = [1,.75,.79];
e9.LineWidth = 2;
hold on
e10 = errorbar(timex,molecules_rne_minus(5,:),molecules_rne_minus(6,:),'d');
e10.MarkerFaceColor = [1,.27,0];
e10.MarkerSize = 5;
e10.Color = [1,.27,0];
e10.LineWidth = 2;
hold on
shadedErrorBar(tv,rneplusm(:,1)+rneplusm(:,3),rneplus_sspread,'lineProps', {'Color',[1,.42,.71]})
hold on
plot(tv,rneplusm(:,1)+rneplusm(:,3),'LineWidth',1,'Color','k')
hold on
shadedErrorBar(tv,zeros(100,1),rneminus_sspread,'lineProps', {'Color',[1,0,0]})
hold on
plot(tv,zeros(100,1),'LineWidth',1,'Color','k')
hold on
e9 = errorbar(timex,molecules_rne_plus(5,:),molecules_rne_plus(6,:),'o');
e9.MarkerFaceColor = [1,.75,.79];
e9.MarkerSize = 5;
e9.Color = [1,.75,.79];
e9.LineWidth = 2;
hold on
e10 = errorbar(timex,molecules_rne_minus(5,:),molecules_rne_minus(6,:),'d');
e10.MarkerFaceColor = [1,.27,0];
e10.MarkerSize = 5;
e10.Color = [1,.27,0];
e10.LineWidth = 2;
lgd = legend('WT SgrS','{\Delta SgrS}');
lgd.FontSize = 10;
lgd.FontWeight = 'bold';
lgd.Location = 'northwest';
title('{\it rne}701 SgrS','FontSize',16,'FontName','Arial','Color','r','FontName','Arial')
xlabel('Time after Induction (min)','FontSize',14,'FontName','Arial')
ylabel('Copy Number','FontSize',14,'FontName','Arial')
set(gca,'YLim',[-10 3.5E6/7096.25])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10)
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[626,281,248,201])
file1 = strcat([folderTitle,'rneptsGsRNA_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'rneptsGsRNA_fit.fig']);
savefig(gcf,file1_fig)

figure(6)
e11 = errorbar(timex,molecules_wt_plus(5,:),molecules_wt_plus(6,:),'s');
e11.MarkerFaceColor = [1,.75,.79];
e11.MarkerSize = 5;
e11.Color = [1,.75,.79];
e11.LineWidth = 2;
hold on
e12 = errorbar(timex,molecules_wt_minus(5,:),molecules_wt_minus(6,:),'*');
e12.MarkerFaceColor = [0,.65,.1];
e12.MarkerSize = 5;
e12.Color = [1,.27,0];
e12.LineWidth = 2;
hold on
shadedErrorBar(tv,wtplusm(:,1)+wtplusm(:,3),wtplus_sspread,'lineProps', {'Color',[1,.42,.71]})
hold on
plot(tv,wtplusm(:,1)+wtplusm(:,3),'LineWidth',1,'Color','k')
hold on
shadedErrorBar(tv,zeros(100,1),wtminus_sspread,'lineProps', {'Color',[1,0,0]})
hold on
plot(tv,zeros(100,1),'LineWidth',1,'Color','k')
hold on
e11 = errorbar(timex,molecules_wt_plus(5,:),molecules_wt_plus(6,:),'s');
e11.MarkerFaceColor = [1,.75,.79];
e11.MarkerSize = 5;
e11.Color = [1,.75,.79];
e3.LineWidth = 2;
hold on
e12 = errorbar(timex,molecules_wt_minus(5,:),molecules_wt_minus(6,:),'*');
e12.MarkerFaceColor = [1,.27,0];
e12.MarkerSize = 5;
e12.Color = [1,.27,0];
e12.LineWidth = 2;
lgd = legend('WT SgrS','{\Delta SgrS}');
lgd.FontSize = 10;
lgd.Location = 'northwest';
title('WT SgrS','FontSize',16,'FontName','Arial','Color','r')
xlabel('Time after Induction (min)','FontSize',14,'FontName','Arial')
ylabel('Copy Number','FontSize',14,'FontName','Arial')
set(gca,'YLim',[-10 3.5E6/7096.25])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10)
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
