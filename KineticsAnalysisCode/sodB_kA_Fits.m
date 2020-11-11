folderTitle = '/Users/reyer/Data/SingleCellEpi/MR243/FixedP51820/';

b_e = .276E-1; %std = 2.23E-3
be_sd = .130E-5;
a_m_plus_wt = 8.771;% std =589.18; 
amp_wt_sd = 0.6542609062;
a_m_plus_rne = 6.683333333; %std = 796.64
amp_rne_sd = 0.1810230188;
a_s_wt = 0.2391; %a_s = 2.29E+03;
a_s_wt_sd = 0.01852608971;
a_s_rne = 0.1683; %a_s = 2.29E+03; 
a_s_rne_sd = 0.01018233765;
b_s_wt = 2.98E-03; 
b_s_wt_sd = 0.000007071067812;
b_s_rne = 1.65E-03; 
b_s_rne_sd = 0.00006687550623;
b_m = 0.005404666667;
b_m_sd = 0.000314038745;
kA = .8575;
kA_sd = 1.060E-9;
k_off = 1.35E+01;
k_off_sd = 10.86426875;
kx_wt_minus = 7.740833333;
kx_wt_plus = 7.262333333; % kx_wt = 0.0040191;
kxw_sd = 1.013540389;
% kx_wt_plus=mode(kx_fit); %kx_wt_plus=0.0033;
kx_rne_minus = 3.293125;
kx_rne_plus = 3.078333333; %kx_rne = 0.0023639;
kxr_sd = 0.3233888269;
kx_s = 0.5938333333; %std = .15
kxs_sd = 	0.5953333333;
b_ms_rne = 0.012117666667; % b_ms = 4.88E-03;
b_ms_wt = 0.012117666667; % b_ms = 4.88E-03;
b_ms_sd = 0.0005819180068;
b_p = .00018;
a_m_minus_wt = 8.781066667; %std = 310.41
amw_sd = 1.202949726;
a_m_minus_rne = 6.827125; %std = 1164.1
amm_sd = 0.9099164797;
% 
% b_e = 2.35E-02; %std = 2.23E-3
% be_sd = 0.009265820705;
% a_m_plus_wt = 8.753333333;% std =589.18; 
% amp_wt_sd = 0.6542609062;
% a_m_plus_rne = 9.030333333; %std = 796.64
% amp_rne_sd = 0.1810230188;
% a_s_wt = 0.2145333333; %a_s = 2.29E+03;
% a_s_wt_sd = 0.009066605392;
% a_s_rne = 0.1920666667; %a_s = 2.29E+03; 
% a_s_rne_sd = 0.01420715782;
% b_s_wt = 2.09E-03; 
% b_s_wt_sd = 0.00001178982612;
% b_s_rne = 1.59E-03; 
% b_s_rne_sd = 0.00006687550623;
% b_m = 0.005676333333;
% b_m_sd = 0.0000830321223;
% kA = 3.23E-05;
% kA_sd = 0.00001910452913;
% k_off = 2.69E-01;
% k_off_sd = 0.2438269127;
% kx_wt_minus = 5.0991;
% kx_wt_plus = 3.6991; % kx_wt = 0.0040191;
% kxw_sd = 0.4656547362;
% % kx_wt_plus=mode(kx_fit); %kx_wt_plus=0.0033;
% kx_rne_minus = 2.257766667;
% kx_rne_plus = 1.091333333; %kx_rne = 0.0023639;
% kxr_sd = 0.04250098038;
% kx_s = 0.0336; %std = .15
% kxs_sd = 	0.1409290602;
% b_ms_rne = 0.006663333333; % b_ms = 4.88E-03;
% b_ms_wt = 0.006663333333; % b_ms = 4.88E-03;
% b_ms_sd = 0.0008603100216;
% b_p = .00018;
% a_m_minus_wt = 12.4053; %std = 310.41
% amw_sd = 1.261931228;
% a_m_minus_rne = 10.72636667; %std = 1164.1
% amm_sd = 1.258103972;

close all
tv = linspace(0, 1800);
tvp = linspace(0, 1440);
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

molecules_wt_plus1 = [709601.6922,537503.346,764685.903,1046192.294,1673136.66,2479907.037,3291319.359,4597902.17;
                      85785.11586,89696.35246,1069219.666,1410525.687,1929558.835,2032196.313,1953665.542,920764.8788;
                      195136.4736,153818.3165,152979.2042,136083.2954,105430.2528,164738.7969,145182.471,311909.1325];
                  
molecules_wt_plus2 = [1313397.907,1125670.049,1425321.897,1153036.313,2077538.263,2336094.806,4485842.258,5582180.8;
                      25159.09882,593584.0755,1400372.742,1178942.461,1596514.889,1398077.813,1723092.936,739851.4549;
                      747298.835,639277.5936,638182.4591,520444.2552,594376.9852,408138.3882,513049.5264,210534.048];

molecules_wt_plus3 = [606596.0203,594100.1408,752934.4109,807808.2311,1692234.02,2537962.205,2634054.939,5582180.8;
                      174776.1687,133085.5829,1290067.901,1112613.131,2082298.936,2162712.435,1787054.665,739851.4549;
                      170751.0349,175518.4587,174960.2044,114106.3714,164581.7131,185604.6023,115970.7151,210534.048];                  

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
                                     
molecules_wt_minus3 = [1607328.013,1760187.181,2132163.243,3417556.799,6328210.867,7808997.984,9082251.066,9935174.519;
                       37928.82483,959568.1077,2625031.333,1126449.883,2791824.164,2045500.094,1777356.339,1933814.577;
                       0,0,0,0,0,0,0,0];      
                   
                   
moleculesrne_plus = [(9/4)*mean((([molecules_rne_plus1(3,1),molecules_rne_plus2(3,1)])+1.7582e+04)/7.8067e+03) mean(((([molecules_rne_plus1(2,1),molecules_rne_plus2(2,1)])+4.3472e+04)/1.6831e+03)+1) 0  mean([molecules_rne_plus1(1,4),molecules_rne_plus2(1,4)])];
moleculeswt_plus = [(9/4)*mean((([molecules_wt_plus1(3,1),molecules_wt_plus2(3,1),molecules_wt_plus3(3,1)])+1.7582e+04)/7.8067e+03) mean(((([molecules_wt_plus1(2,1),molecules_wt_plus2(2,1),molecules_wt_plus3(2,1)])+4.3472e+04)/1.6831e+03)+1) 0  mean([molecules_wt_plus1(1,4),molecules_wt_plus2(1,4),molecules_wt_plus3(1,4)])];
moleculesrne_minus = [(9/4)*mean((([molecules_rne_minus1(3,1),molecules_rne_minus2(3,1),molecules_rne_minus3(3,1),molecules_rne_minus4(3,1)])+1.7582e+04)/7.8067e+03) mean(((([molecules_rne_minus1(2,1),molecules_rne_minus2(2,1),molecules_rne_minus3(2,1),molecules_rne_minus4(2,1)])+4.3472e+04)/1.6831e+03)+1) 0 mean([molecules_rne_minus1(1,4),molecules_rne_minus2(1,4),molecules_rne_minus3(1,4),molecules_rne_minus4(1,4)])];
moleculeswt_minus = [(9/4)*mean((([molecules_wt_minus1(3,1),molecules_wt_minus2(3,1),molecules_wt_minus3(3,1)])+1.7582e+04)/7.8067e+03) mean(((([molecules_wt_minus1(2,1),molecules_wt_minus2(2,1),molecules_wt_minus3(2,1)])+4.3472e+04)/1.6831e+03)+1) 0 mean([molecules_wt_minus1(1,4),molecules_wt_minus2(1,4),molecules_wt_minus3(1,4)])];

std_rne_plus = [(9/4)*std((([molecules_rne_plus1(3,1),molecules_rne_plus2(3,1)])+1.7582e+04)/7.8067e+03) std(((([molecules_rne_plus1(2,1),molecules_rne_plus2(2,1)])+4.3472e+04)/1.6831e+03)+1) 0  std([molecules_rne_plus1(1,4),molecules_rne_plus2(1,4)])];
std_wt_plus = [(9/4)*std((([molecules_wt_plus1(3,1),molecules_wt_plus2(3,1),molecules_wt_plus3(3,1)])+1.7582e+04)/7.8067e+03) std(((([molecules_wt_plus1(2,1),molecules_wt_plus2(2,1),molecules_wt_plus3(2,1)])+4.3472e+04)/1.6831e+03)+1) 0  std([molecules_wt_plus1(1,4),molecules_wt_plus2(1,4),molecules_wt_plus3(1,4)])];
std_rne_minus = [(9/4)*std((([molecules_rne_minus1(3,1),molecules_rne_minus2(3,1),molecules_rne_minus3(3,1),molecules_rne_minus4(3,1)])+1.7582e+04)/7.8067e+03) std(((([molecules_rne_minus1(2,1),molecules_rne_minus2(2,1),molecules_rne_minus3(2,1),molecules_rne_minus4(2,1)])+4.3472e+04)/1.6831e+03)+1) 0 std([molecules_rne_minus1(1,4),molecules_rne_minus2(1,4),molecules_rne_minus3(1,4),molecules_rne_minus4(1,4)])];
std_wt_minus = [(9/4)*std((([molecules_wt_minus1(3,1),molecules_wt_minus2(3,1),molecules_wt_minus3(3,1)])+1.7582e+04)/7.8067e+03) std(((([molecules_wt_minus1(2,1),molecules_wt_minus2(2,1),molecules_wt_minus3(2,1)])+4.3472e+04)/1.6831e+03)+1) 0 std([molecules_wt_minus1(1,4),molecules_wt_minus2(1,4),molecules_wt_minus3(1,4)])];

rneplusp = zeros(100,500);
rneplusm = zeros(100,500);
rnepluss = zeros(100,500);

wtplusp = zeros(100,500);
wtplusm = zeros(100,500);
wtpluss = zeros(100,500);

rneminusp = zeros(100,400);
rneminusm = zeros(100,400);
rneminuss = zeros(100,400);

wtminusp = zeros(100,400);
wtminusm = zeros(100,400);
wtminuss = zeros(100,400);

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


rneplusp(:,1:4) = KA_fit(tvp,0,a_m_plus_rne,a_s_rne,b_s_rne,b_m,kA,kx_rne_plus,kx_s*kx_rne_plus,b_ms_rne,b_p,moleculesrne_plus);
rneplusm(:,1:4) = KA_fit(tv,0,a_m_plus_rne,a_s_rne,b_s_rne,b_m,kA,kx_rne_plus,kx_s*kx_rne_plus,b_ms_rne,b_p,moleculesrne_plus);

wtplusp(:,1:4) = KA_fit(tvp,b_e,a_m_plus_wt,a_s_wt,b_s_wt,b_m,kA,kx_wt_plus,kx_s*kx_wt_plus,b_ms_wt,b_p,moleculeswt_plus);
wtplusm(:,1:4) = KA_fit(tv,b_e,a_m_plus_wt,a_s_wt,b_s_wt,b_m,kA,kx_wt_plus,kx_s*kx_wt_plus,b_ms_wt,b_p,moleculeswt_plus);

rneminusp(:,1:4) = nosRNA_fit([a_m_minus_rne,kx_rne_minus],tvp,0,0,0,0,0,b_m,0,0,b_p,moleculesrne_minus);
rneminusm(:,1:4) = nosRNA_fit([a_m_minus_rne,kx_rne_minus],tv,0,0,0,0,0,b_m,0,0,b_p,moleculesrne_minus);

wtminusp(:,1:4) = nosRNA_fit([a_m_minus_wt,kx_wt_minus],tvp,0,0,0,0,0,b_m,0,0,b_p,moleculeswt_minus);
wtminusm(:,1:4) = nosRNA_fit([a_m_minus_wt,kx_wt_minus],tv,0,0,0,0,0,b_m,0,0,b_p,moleculeswt_minus);

i_rand = 5:8;
i_rand_check = [];
for i = 1:19
    
    a1 = KA_fit(tvp,0,normrnd(a_m_plus_rne,amp_rne_sd),normrnd(a_s_rne,a_s_rne_sd),normrnd(b_s_rne,b_s_rne_sd),normrnd(b_m,b_m_sd),normrnd(kA,kA_sd),normrnd(kx_rne_plus,kxr_sd),normrnd(kx_s,kxs_sd)*kx_rne_plus,normrnd(b_ms_rne,b_ms_sd),b_p,normrnd(moleculesrne_plus,std_rne_plus));
    if length(a1(:,1)) == 100
        rneplusp(:,i_rand) = a1;
    else
        i_rand_check = [i_rand_check i_rand];
    end
    
    a2 = KA_fit(tv,0,normrnd(a_m_plus_rne,amp_rne_sd),normrnd(a_s_rne,a_s_rne_sd),normrnd(b_s_rne,b_s_rne_sd),normrnd(b_m,b_m_sd),normrnd(kA,kA_sd),normrnd(kx_rne_plus,kxr_sd),normrnd(kx_s,kxs_sd)*kx_rne_plus,normrnd(b_ms_rne,b_ms_sd),b_p,normrnd(moleculesrne_plus,std_rne_plus));
    if length(a2(:,1)) == 100
        rneplusm(:,i_rand) = a2;
    else
        i_rand_check = [i_rand_check i_rand];
    end
    
    a3 = KA_fit(tvp,normrnd(b_e,be_sd),normrnd(a_m_plus_wt,amp_wt_sd),normrnd(a_s_wt,a_s_wt_sd),normrnd(b_s_wt,b_s_wt_sd),normrnd(b_m,b_m_sd),normrnd(kA,kA_sd),normrnd(kx_wt_plus,kxw_sd),normrnd(kx_s,kxs_sd)*kx_wt_plus,normrnd(b_ms_wt,b_ms_sd),b_p,normrnd(moleculeswt_plus,std_wt_plus));
    if length(a3(:,1)) == 100
        wtplusp(:,i_rand) = a3;
    else
        i_rand_check = [i_rand_check i_rand];
    end
    
    a4 = KA_fit(tv,normrnd(b_e,be_sd),normrnd(a_m_plus_wt,amp_wt_sd),normrnd(a_s_wt,a_s_wt_sd),normrnd(b_s_wt,b_s_wt_sd),normrnd(b_m,b_m_sd),normrnd(kA,kA_sd),normrnd(kx_wt_plus,kxw_sd),normrnd(kx_s,kxs_sd)*kx_wt_plus,normrnd(b_ms_wt,b_ms_sd),b_p,normrnd(moleculeswt_plus,std_wt_plus));
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
wt_plus_sRNA_error = zeros(100,120);

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
    
    rneplus_mspread(i) = std(nonOutliers);
    
    percentiles = [];
    nonOutliers = [];
    percentiles = prctile(rne_plus_sRNA_error(i,:),[5,95]);
    outlierIndexes = rne_plus_sRNA_error(i,:) < percentiles(1) | rne_plus_sRNA_error(i,:) > percentiles(2);
    nonOutliers = rne_plus_sRNA_error(i,~outlierIndexes);
    
    rneplus_sspread(i) = std(nonOutliers);
    
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
    
    percentiles = [];
    nonOutliers = [];
    percentiles = prctile(rne_minus_sRNA_error(i,:),[5,95]);
    outlierIndexes = rne_minus_sRNA_error(i,:) < percentiles(1) | rne_minus_sRNA_error(i,:) > percentiles(2);
    nonOutliers = rne_minus_sRNA_error(i,~outlierIndexes);
    
    rneminus_sspread(i) = std(nonOutliers);
    
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
    
    percentiles = [];
    nonOutliers = [];
    percentiles = prctile(wt_minus_sRNA_error(i,:),[5,95]);
    outlierIndexes = wt_minus_sRNA_error(i,:) < percentiles(1) | wt_minus_sRNA_error(i,:) > percentiles(2);
    nonOutliers = wt_minus_sRNA_error(i,~outlierIndexes);
    
    wtminus_sspread(i) = std(nonOutliers);
end  


molecules_rne_plus = zeros(6,8);
molecules_rne_minus = zeros(6,8);
molecules_wt_plus = zeros(6,8);
molecules_wt_minus = zeros(6,8);
for i = 1:8
    molecules_rne_plus(1,i) = mean([molecules_rne_plus1(1,i),molecules_rne_plus2(1,i)]);
    molecules_rne_minus(1,i) = mean([molecules_rne_minus1(1,i),molecules_rne_minus2(1,i),molecules_rne_minus3(1,i),molecules_rne_minus4(1,i)]); 
    molecules_wt_plus(1,i) = mean([molecules_wt_plus1(1,i),molecules_wt_plus2(1,i),molecules_wt_plus3(1,i)]);
    molecules_wt_minus(1,i) = mean([molecules_wt_minus1(1,i),molecules_wt_minus2(1,i),molecules_wt_minus3(1,i)]);
    
    molecules_rne_plus(2,i) = std([molecules_rne_plus1(1,i),molecules_rne_plus2(1,i)])/1.723;
    molecules_rne_minus(2,i) = std([molecules_rne_minus1(1,i),molecules_rne_minus2(1,i),molecules_rne_minus3(1,i),molecules_rne_minus4(1,i)])/2; 
    molecules_wt_plus(2,i) = std([molecules_wt_plus1(1,i),molecules_wt_plus2(1,i),molecules_wt_plus3(1,i)])/1.723;
    molecules_wt_minus(2,i) = std([molecules_wt_minus1(1,i),molecules_wt_minus2(1,i),molecules_wt_minus3(1,i)])/2;
    
    molecules_rne_plus(3,i) = (mean((([molecules_rne_plus1(2,i),molecules_rne_plus2(2,i)])+4.3472e+04)/1.6831e+03)+1);
    molecules_rne_minus(3,i) =( mean((([molecules_rne_minus1(2,i),molecules_rne_minus2(2,i),molecules_rne_minus3(2,i),molecules_rne_minus4(2,i)])+4.3472e+04)/1.6831e+03)+1); 
    molecules_wt_plus(3,i) = (mean((([molecules_wt_plus1(2,i),molecules_wt_plus2(2,i),molecules_wt_plus3(2,i)])+4.3472e+04)/1.6831e+03)+1);
    molecules_wt_minus(3,i) =( mean((([molecules_wt_minus1(2,i),molecules_wt_minus2(2,i),molecules_wt_minus3(2,i)])+4.3472e+04)/1.6831e+03)+1);
    
    molecules_rne_plus(4,i) = (std(((([molecules_rne_plus1(2,i),molecules_rne_plus2(2,i)])+4.3472e+04)/1.6831e+03)+1))/1.723;
    molecules_rne_minus(4,i) = (std(((([molecules_rne_minus1(2,i),molecules_rne_minus2(2,i),molecules_rne_minus3(2,i),molecules_rne_minus4(2,i)])+4.3472e+04)/1.6831e+03)+1))/2; 
    molecules_wt_plus(4,i) = (std(((([molecules_wt_plus1(2,i),molecules_wt_plus2(2,i),molecules_wt_plus3(2,i)])+4.3472e+04)/1.6831e+03)+1))/1.723;
    molecules_wt_minus(4,i) = (std(((([molecules_wt_minus1(2,i),molecules_wt_minus2(2,i),molecules_wt_minus3(2,i)])+4.3472e+04)/1.6831e+03)+1))/2;
    
    molecules_rne_plus(5,i) = (9/4)*mean((([molecules_rne_plus1(3,i),molecules_rne_plus2(3,i)])+1.7582e+04)/7.8067e+03);
    molecules_rne_minus(5,i) = (9/4)*mean((([molecules_rne_minus1(3,i),molecules_rne_minus2(3,i),molecules_rne_minus3(3,i),molecules_rne_minus4(3,i)])+1.7582e+04)/7.8067e+03); 
    molecules_wt_plus(5,i) = (9/4)*mean((([molecules_wt_plus1(3,i),molecules_wt_plus2(3,i),molecules_wt_plus3(3,i)])+1.7582e+04)/7.8067e+03);
    molecules_wt_minus(5,i) = (9/4)*mean((([molecules_wt_minus1(3,i),molecules_wt_minus2(3,i),molecules_wt_minus3(3,i)])+1.7582e+04)/7.8067e+03);
    
    molecules_rne_plus(6,i) = std((([molecules_rne_plus1(3,i),molecules_rne_plus2(3,i)])+1.7582e+04)/7.8067e+03);
    molecules_rne_minus(6,i) = std((([molecules_rne_minus1(3,i),molecules_rne_minus2(3,i),molecules_rne_minus3(3,i),molecules_rne_minus4(3,i)])+1.7582e+04)/7.8067e+03)/2; 
    molecules_wt_plus(6,i) = std((([molecules_wt_plus1(3,i),molecules_wt_plus2(3,i),molecules_wt_plus3(3,i)])+1.7582e+04)/7.8067e+03)/1.723;
    molecules_wt_minus(6,i) = std((([molecules_wt_minus1(3,i),molecules_wt_minus2(3,i),molecules_wt_minus3(3,i)])+1.7582e+04)/7.8067e+03)/2;
    
end


time_points = round([0,1,3,6,12,18,24]*3.3333)+1;
time_points = [time_points,100];

time_points_p = round([6,12,18]*4.1667)+1;
time_points_p = [time_points_p,100];




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
lgd = legend('WT RyhB','{\Delta RyhB}');
lgd.FontSize = 10;
lgd.FontWeight = 'bold';
lgd.Location = 'northwest';
title('rne701{\it sodB130+30-sfGFP}, fixed {\alpha_{m}}','FontSize',14,'FontWeight','bold','Color','g','FontName','Arial')
xlabel('Time after Induction (min)','FontSize',12,'FontName','Arial')
ylabel('Copy Number ','FontSize',18,'FontName','Arial')
set(gca,'YLim',[-100 3.5E6/1125.7])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10,'FontWeight','bold')
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[626,281,248,201])
file1 = strcat([folderTitle,'rnesodB130+30mRNA_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'rnesodB130+30mRNA_fit.fig']);
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
shadedErrorBar(tv,wtplusm(:,2)+wtplusm(:,3),wtplus_mspread,'lineProps', {'Color',[0,.40,0]})
hold on
plot(tv,wtplusm(:,2)+wtplusm(:,3),'LineWidth',1,'Color','k')
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
lgd = legend('WT RyhB','{\Delta RyhB}');
lgd.FontSize = 10;
lgd.FontWeight = 'bold';
lgd.Location = 'northwest';
title('WT{\it sodB130+30-sfGFP}, fixed {\alpha_{m}}','FontSize',18,'FontWeight','bold','Color','g','FontName','Arial')
xlabel('Time after Induction (min)','FontSize',18,'FontName','Arial')
ylabel('Copy Number','FontSize',18,'FontName','Arial')
set(gca,'YLim',[-100 3.5E6/1125.7])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10,'FontWeight','bold')
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[626,281,248,201])
file1 = strcat([folderTitle,'WTsodB130+30mRNA_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'WTsodB130+30mRNA_fit.fig']);
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
lgd = legend('WT RyhB','{\Delta RyhB}');
lgd.FontSize = 10;
lgd.FontWeight = 'bold';
lgd.Location = 'northwest';
title('rne701 sfGFP, fixed {\alpha_{m}}','FontSize',18,'FontWeight','bold','Color','b','FontName','Arial')
xlabel('Time after Induction (min)','FontSize',18,'FontWeight','bold')
ylabel('Fluorescence (A.U)','FontSize',18,'FontName','Arial')
set(gca,'YLim',[-10000 14E6])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10,'FontWeight','bold')
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[626,281,248,201])
file1 = strcat([folderTitle,'rnesodB130+30Protein_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'rnesodB130+30Protein_fit.fig']);
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
lgd = legend('WT RyhB','{\Delta RyhB}');
lgd.FontSize = 10;
lgd.FontWeight = 'bold';
lgd.Location = 'northwest';
title('WT sfGFP, fixed {\alpha_{m}}','FontSize',14,'FontWeight','bold','Color','b')
xlabel('Time after Induction (min)','FontSize',12,'FontWeight','bold')
ylabel('Fluorescence (A.U)','FontSize',12,'FontWeight','bold')
set(gca,'YLim',[-10000 14E6])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10,'FontWeight','bold')
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[626,281,248,201])
file1 = strcat([folderTitle,'WTsodB130+30Protein_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'WTsodB130+30Protein_fit.fig']);
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
shadedErrorBar(tv,rneminusm(:,1)+rneminusm(:,3),rneminus_sspread,'lineProps', {'Color',[1,0,0]})
hold on
plot(tv,rneminusm(:,1)+rneminusm(:,3),'LineWidth',1,'Color','k')
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
lgd = legend('WT RyhB','{\Delta RyhB}');
lgd.FontSize = 10;
lgd.FontWeight = 'bold';
lgd.Location = 'northwest';
title('rne701 RyhB, fixed {\alpha_{m}}','FontSize',18,'FontWeight','bold','Color','r','FontName','Arial')
xlabel('Time after Induction (min)','FontSize',18,'FontName','Arial')
ylabel('Copy Number','FontSize',18,'FontName','Arial')
set(gca,'YLim',[-10 300])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10,'FontWeight','bold')
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[626,281,248,201])
file1 = strcat([folderTitle,'rnesodB130+30sRNA_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'rnesodB130+30sRNA_fit.fig']);
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
shadedErrorBar(tv,wtminusm(:,1)+wtminusm(:,3),wtminus_sspread,'lineProps', {'Color',[1,0,0]})
hold on
plot(tv,wtminusm(:,1)+wtminusm(:,3),'LineWidth',1,'Color','k')
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
lgd = legend('WT RyhB','{\Delta RyhB}');
lgd.FontSize = 10;
lgd.FontWeight = 'bold';
lgd.Location = 'northwest';
title('WT RyhB, fixed {\alpha_{m}}','FontSize',18,'FontWeight','bold','Color','r')
xlabel('Time after Induction (min)','FontSize',18,'FontName','Arial')
ylabel('Copy Number','FontSize',18,'FontName','Arial')
set(gca,'YLim',[-10 300])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10,'FontWeight','bold')
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[626,281,248,201])
file1 = strcat([folderTitle,'WTsodB130+30sRNA_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'WTsodB130+30sRNA_fit.fig']);
savefig(gcf,file1_fig)

copy_array = [kA,kx_s,b_e,a_m_plus_wt,a_m_plus_rne,kx_wt_plus,kx_rne_plus,b_m,b_ms_wt,a_m_minus_wt,a_m_minus_rne,kx_wt_minus,kx_rne_minus;
    kA_sd,kxs_sd,be_sd,amp_wt_sd,amp_rne_sd,kxw_sd,kxr_sd,b_m_sd,b_ms_sd,amw_sd,amm_sd,kxw_sd,kxr_sd];

copy_table = array2table(copy_array);copy_table.Properties.VariableNames = {'kA' 'kx_s' 'b_e' 'a_m_wt_plus' 'a_m_rne_plus' 'kx_wt_plus' 'kx_rne_plus' 'b_m' 'b_ms' 'a_m_wt_minus' 'a_m_rne_minus' 'kx_wt_minus' 'kx_rne_minus'};

copy_table_file = strcat([folderTitle,'parameters.csv']);

writetable(copy_table,copy_table_file);
