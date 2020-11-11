folderTitle = '/Users/reyer/Data/SingleCellEpi/MR243/Post8420/';

b_e = 8.58E-01; %std = 2.23E-3
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

k_on = 8.64E-05;
k_on_sd = 0.00001707;
k_off = 7.56;
k_off_sd = 0.884;

k_nuc = 0;
k_nuc_sd = 0;
k_off_nuc = 0;
k_off_nuc_sd = 0;

b_ms = 0.006055; % b_ms = 4.88E-03;
b_ms_sd = 0.0000645;
b_p = .00018;

b_nuc_wt = b_e + b_ms;
b_nuc_wt_sd = 0;
b_nuc_rne = b_ms;
b_nuc_rne_sd = 0;

kx_wt_minus =  28.707;
kx_wt_plus = 26.535; % kx_wt = 0.0040191;
kxw_sd = 2.19;

kx_rne_minus = 10.355625;
kx_rne_plus = 9.105; %kx_rne = 0.0023639;
kxr_sd = 1.757;

kx_s = 0.469; %std = .15
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

molecules_wt_plus1 = [1130948.676,892617.3287,1074207.349,1375535.856,2229638.274,2729098.729,3567882.043,5582180.8;
                      10050.51839,305731.8887,803354.734,860573.0497,284796.5291,881811.9806,918839.9995,920764.8788;
                      430758.4504,399694.8762,385861.0661,425008.5325,117141.0959,317523.1091,392639.7222,311909.1325];                  

molecules_wt_plus2 = [1077656.5,1109629.183,1219469.886,1555835.48,2574231.765,3197925.634,4967996.936,5582180.8;
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
moleculeswt_plus = [(9/4)*mean((([molecules_wt_plus1(3,1),molecules_wt_plus2(3,1)])+5.2879e+04)/7.3302e+03) 0 0 mean(((([molecules_wt_plus1(2,1),molecules_wt_plus2(2,1),molecules_wt_plus3(2,1)])+1.3257e+05)/6.0368e+03)+1)/20 0 mean([molecules_wt_plus1(1,4),molecules_wt_plus2(1,4)])];
moleculesrne_minus = [0 mean(((([molecules_rne_minus1(2,1),molecules_rne_minus2(2,1)])+1.3257e+05)/6.0368e+03)+1) mean([molecules_rne_minus1(1,4),molecules_rne_minus2(1,4)])];
moleculeswt_minus = [0 mean(((([molecules_wt_minus1(2,1),molecules_wt_minus2(2,1)])+1.3257e+05)/6.0368e+03)+1)/20 mean([molecules_wt_minus1(1,4),molecules_wt_minus2(1,4)])];

std_rne_plus = [(9/4)*std((([molecules_rne_plus1(3,1),molecules_rne_plus2(3,1)])+5.2879e+04)/7.3302e+03) 0 0  std(((([molecules_rne_plus1(2,1),molecules_rne_plus2(2,1)])+1.3257e+05)/6.0368e+03)+1) 0 std([molecules_rne_plus1(1,4),molecules_rne_plus2(1,4)])];
std_wt_plus = [(9/4)*std((([molecules_wt_plus1(3,1),molecules_wt_plus2(3,1)])+5.2879e+04)/7.3302e+03) 0 0 std(((([molecules_wt_plus1(2,1),molecules_wt_plus2(2,1)])+1.3257e+05)/6.0368e+03)+1) 0 std([molecules_wt_plus1(1,4),molecules_wt_plus2(1,4)])];
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

i_rand = 7:12;
i_rand2 = 4:6;
i_rand_check = [];
for i = 1:19
    
    a1 = bestInit_fit(tvp,normrnd(a_s_rne,a_s_rne_sd),normrnd(b_s_rne,b_s_rne_sd),normrnd(k_init_rne,k_init_rne_sd),k_elon,normrnd(k_elon_prime,k_elon_prime_sd),normrnd(b_m,b_m_sd),normrnd(k_on,k_on_sd),normrnd(k_off,k_off_sd),normrnd(kx_rne_plus,kxr_sd),normrnd(kx_s,kxs_sd)*kx_rne_plus,normrnd(b_ms,b_ms_sd),0,b_p,normrnd(k_nuc,k_nuc_sd),normrnd(b_nuc_rne,b_nuc_rne_sd),n_plasmids,normrnd(k_off_nuc,k_off_nuc_sd),normrnd(moleculesrne_plus,std_rne_plus));
    if length(a1(:,1)) == 100
        rneplusp(:,i_rand) = a1;
    else
        i_rand_check = [i_rand_check i_rand];
    end
    
    a2 = bestInit_fit(tv,normrnd(a_s_rne,a_s_rne_sd),normrnd(b_s_rne,b_s_rne_sd),normrnd(k_init_rne,k_init_rne_sd),k_elon,normrnd(k_elon_prime,k_elon_prime_sd),normrnd(b_m,b_m_sd),normrnd(k_on,k_on_sd),normrnd(k_off,k_off_sd),normrnd(kx_rne_plus,kxr_sd),normrnd(kx_s,kxs_sd)*kx_rne_plus,normrnd(b_ms,b_ms_sd),0,b_p,normrnd(k_nuc,k_nuc_sd),normrnd(b_nuc_rne,b_nuc_rne_sd),n_plasmids,normrnd(k_off_nuc,k_off_nuc_sd),normrnd(moleculesrne_plus,std_rne_plus));
    if length(a2(:,1)) == 100
        rneplusm(:,i_rand) = a2;
    else
        i_rand_check = [i_rand_check i_rand];
    end
    
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
   
    a5 = initiation_fit([normrnd(k_init_rne,k_init_rne_sd),normrnd(kx_rne_minus,kxr_sd)],tvp,k_elon,normrnd(b_m,b_m_sd),b_p,n_plasmids,normrnd(moleculesrne_minus,std_rne_minus));
    if length(a5(:,1)) == 100
        rneminusp(:,i_rand2) = a5;
    else
        i_rand_check = [i_rand_check i_rand2];
    end
    
    a6 = initiation_fit([normrnd(k_init_rne,k_init_rne_sd),normrnd(kx_rne_minus,kxr_sd)],tv,k_elon,normrnd(b_m,b_m_sd),b_p,n_plasmids,normrnd(moleculesrne_minus,std_rne_minus));
    if length(a6(:,1)) == 100
        rneminusm(:,i_rand2) = a6;
    else
        i_rand_check = [i_rand_check i_rand2];
    end
    
    a7 = initiation_fit([normrnd(k_init_wt,k_init_wt_sd),normrnd(kx_wt_minus,kxr_sd)],tvp,k_elon,normrnd(b_m,b_m_sd),b_p,n_plasmids,normrnd(moleculeswt_minus,std_wt_minus));
    if length(a7(:,1)) == 100
        wtminusp(:,i_rand2) = a7;
    else
        i_rand_check = [i_rand_check i_rand2];
    end
    
    a8 = initiation_fit([normrnd(k_init_wt,k_init_wt_sd),normrnd(kx_wt_minus,kxr_sd)],tv,k_elon,normrnd(b_m,b_m_sd),b_p,n_plasmids,normrnd(moleculeswt_minus,std_wt_minus));
    if length(a8(:,1)) == 100
        wtminusm(:,i_rand2) = a8;
    else
        i_rand_check = [i_rand_check i_rand];
    end
    
    i_rand = i_rand + 6;
    i_rand2 = i_rand2 + 3;
    
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
    for j = 4:6:120
        
        rne_plus_mRNA_error(i,ij) = rneplusm(i,j) + rneplusm(i,j+1);
        
        wt_plus_mRNA_error(i,ij) = wtplusm(i,j) + wtplusm(i,j+1);
        ij = ij + 1;
        
    end

    ik = 1;
    for k = 6:6:120
        
        rne_plus_protein_error(i,k) = rneplusp(i,k);
        
        wt_plus_protein_error(i,k) = wtplusp(i,k);
        ik = ik + 1;
    end
    
    il = 1;
    for j = 1:6:120
        
        rne_plus_sRNA_error(i,il) = rneplusm(i,j) + rneplusm(i,j+2) + rneplusm(i,j+4);
        
        wt_plus_sRNA_error(i,il) = wtplusm(i,j) + wtplusm(i,j+2) + wtplusm(i,j+4);
        il = il + 1;
        
    end
    
end

for i = 1:100
    ij = 1;
    for j = 2:3:60
        rne_minus_mRNA_error(i,ij) = rneminusm(i,j);
        
        wt_minus_mRNA_error(i,ij) = wtminusm(i,j);
        
        ij = ij + 1;
        
    end

    ik = 1;
    for k = 3:3:60
        rne_minus_protein_error(i,k) = rneminusp(i,k);
        
        wt_minus_protein_error(i,k) = wtminusp(i,k);
        
        ik = ik + 1;
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
    
    rneplus_sspread(i) = std(nonOutliers)/1;
    
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
    molecules_rne_plus(1,i) = mean([molecules_rne_plus1(1,i),molecules_rne_plus2(1,i)]);
    molecules_rne_minus(1,i) = mean([molecules_rne_minus1(1,i),molecules_rne_minus2(1,i),molecules_rne_minus3(1,i),molecules_rne_minus4(1,i)]); 
    molecules_wt_plus(1,i) = mean([molecules_wt_plus1(1,i),molecules_wt_plus2(1,i)]);
    molecules_wt_minus(1,i) = mean([molecules_wt_minus1(1,i),molecules_wt_minus2(1,i)]);
    
    molecules_rne_plus(2,i) = std([molecules_rne_plus1(1,i),molecules_rne_plus2(1,i)])/1.414;
    molecules_rne_minus(2,i) = std([molecules_rne_minus1(1,i),molecules_rne_minus2(1,i),molecules_rne_minus3(1,i),molecules_rne_minus4(1,i)])/2; 
    molecules_wt_plus(2,i) = std([molecules_wt_plus1(1,i),molecules_wt_plus2(1,i)])/1.414;
    molecules_wt_minus(2,i) = std([molecules_wt_minus1(1,i),molecules_wt_minus2(1,i)])/1.414;
    
    molecules_rne_plus(3,i) = (mean((([molecules_rne_plus1(2,i),molecules_rne_plus2(2,i)])+1.3257e+05)/6.0368e+03)+1);
    molecules_rne_minus(3,i) =( mean((([molecules_rne_minus1(2,i),molecules_rne_minus2(2,i),molecules_rne_minus3(2,i),molecules_rne_minus4(2,i)])+1.3257e+05)/6.0368e+03)+1); 
    molecules_wt_plus(3,i) = (mean((([molecules_wt_plus1(2,i),molecules_wt_plus2(2,i)])+1.3257e+05)/6.0368e+03)+1);
    molecules_wt_minus(3,i) =( mean((([molecules_wt_minus1(2,i),molecules_wt_minus2(2,i)])+1.3257e+05)/6.0368e+03)+1);
    
    molecules_rne_plus(4,i) = (std(((([molecules_rne_plus1(2,i),molecules_rne_plus2(2,i)])+1.3257e+05)/6.0368e+03)+1))/414;
    molecules_rne_minus(4,i) = (std(((([molecules_rne_minus1(2,i),molecules_rne_minus2(2,i),molecules_rne_minus3(2,i),molecules_rne_minus4(2,i)])+1.3257e+05)/6.0368e+03)+1))/2; 
    molecules_wt_plus(4,i) = (std(((([molecules_wt_plus1(2,i),molecules_wt_plus2(2,i)])+1.3257e+05)/6.0368e+03)+1))/1.414;
    molecules_wt_minus(4,i) = (std(((([molecules_wt_minus1(2,i),molecules_wt_minus2(2,i)])+1.3257e+05)/6.0368e+03)+1))/1.414;
    
    molecules_rne_plus(5,i) = (9/4)*mean((([molecules_rne_plus1(3,i),molecules_rne_plus2(3,i)])+5.2879e+04)/7.3302e+03);
    molecules_rne_minus(5,i) = (9/4)*mean((([molecules_rne_minus1(3,i),molecules_rne_minus2(3,i),molecules_rne_minus3(3,i),molecules_rne_minus4(3,i)])+5.2879e+04)/7.3302e+03); 
    molecules_wt_plus(5,i) = (9/4)*mean((([molecules_wt_plus1(3,i),molecules_wt_plus2(3,i)])+5.2879e+04)/7.3302e+03);
    molecules_wt_minus(5,i) = (9/4)*mean((([molecules_wt_minus1(3,i),molecules_wt_minus2(3,i)])+5.2879e+04)/7.3302e+03);
    
    molecules_rne_plus(6,i) = std((([molecules_rne_plus1(3,i),molecules_rne_plus2(3,i)])+5.2879e+04)/7.3302e+03)/1.414;
    molecules_rne_minus(6,i) = std((([molecules_rne_minus1(3,i),molecules_rne_minus2(3,i),molecules_rne_minus3(3,i),molecules_rne_minus4(3,i)])+5.2879e+04)/7.3302e+03)/2; 
    molecules_wt_plus(6,i) = std((([molecules_wt_plus1(3,i),molecules_wt_plus2(3,i)])+5.2879e+04)/7.3302e+03)/1.414;
    molecules_wt_minus(6,i) = std((([molecules_wt_minus1(3,i),molecules_wt_minus2(3,i)])+5.2879e+04)/7.3302e+03)/1.414;
    
end


% figure(1)
% 
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
shadedErrorBar(tv,rneplusm(:,4)+rneplusm(:,5),rneplus_mspread,'lineProps', {'Color',[0,.40,0]})
hold on
plot(tv,rneplusm(:,4)+rneplusm(:,5),'LineWidth',1,'Color','k')
hold on
shadedErrorBar(tv,rneminusm(:,2),rneminus_mspread,'lineProps', {'Color',[0,.95,0]})
hold on
plot(tv,rneminusm(:,2),'LineWidth',1,'Color','k')
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
lgd.Location = 'northwest';
title('{\it rne701_{130+30}-sfGFP}','FontSize',16,'Color',[0,.40,0],'FontName','Arial')
xlabel('Time after Induction (min)','FontSize',14,'FontName','Arial')
ylabel('Copy Number ','FontSize',14,'FontName','Arial')
set(gca,'YLim',[-10 800])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10)
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
shadedErrorBar(tv,(wtplusm(:,4)+wtplusm(:,5)),wtplus_mspread,'lineProps', {'Color',[0,.40,0]})
hold on
plot(tv,(wtplusm(:,4)+wtplusm(:,5)),'LineWidth',1,'Color','k')
hold on
shadedErrorBar(tv,wtminusm(:,2),wtminus_mspread,'lineProps', {'Color',[0,.95,0]})
hold on
plot(tv,wtminusm(:,2),'LineWidth',1,'Color','k')
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
lgd.Location = 'northwest';
title('WT{\it sodB_{130+30}-sfGFP}','FontSize',16,'Color',[0,.40,0],'FontName','Arial')
xlabel('Time after Induction (min)','FontSize',14,'FontName','Arial')
ylabel('Copy Number','FontSize',14,'FontName','Arial')
set(gca,'YLim',[-10 800])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10)
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
set(gcf,'position',[835,883,868,667])
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
shadedErrorBar(tvp+360,rneplusp(:,6),rneplus_pspread,'lineProps', {'Color',[.05,.05,.38]})
hold on
plot(tvp+360,rneplusp(:,6),'LineWidth',1,'Color','k')
hold on
shadedErrorBar(tvp+360,rneminusp(:,3),rneminus_pspread,'lineProps', {'Color',[0,.95,.95]})
hold on
plot(tvp+360,rneminusp(:,3),'LineWidth',1,'Color','k')
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
lgd.Location = 'northwest';
title('{\it rne}701 sfGFP','FontSize',16,'Color','b','FontName','Arial')
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
shadedErrorBar(tvp+360,wtplusp(:,6),wtplus_pspread,'lineProps', {'Color',[.05,.05,.38]})
hold on
plot(tvp+360,wtplusp(:,6),'LineWidth',1,'Color','k')
hold on
shadedErrorBar(tvp+360,wtminusp(:,3),wtminus_pspread,'lineProps', {'Color',[0,.95,.95]})
hold on
plot(tvp+360,wtminusp(:,3),'LineWidth',1,'Color','k')
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
lgd.Location = 'northwest';
title('WT sfGFP','FontSize',16,'Color','b','FontName','Arial')
xlabel('Time after Induction (min)','FontSize',14,'FontName','Arial')
ylabel('Fluorescence (A.U)','FontSize',14,'FontName','Arial')
set(gca,'YLim',[-10000 16E6])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10)
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
set(gcf,'position',[835,883,868,667])
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
shadedErrorBar(tv,rneplusm(:,1)+rneplusm(:,3)+rneplusm(:,5),rneplus_sspread,'lineProps', {'Color',[1,.42,.71]})
hold on
plot(tv,rneplusm(:,1)+rneplusm(:,3)+rneplusm(:,5),'LineWidth',1,'Color','k')
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
lgd = legend('WT RyhB','{\Delta RyhB}');
lgd.FontSize = 10;
lgd.Location = 'northwest';
title('{\it rne}701 RyhB','FontSize',16,'Color','r','FontName','Arial')
xlabel('Time after Induction (min)','FontSize',14,'FontName','Arial')
ylabel('Copy Number','FontSize',14,'FontName','Arial')
set(gca,'YLim',[-10 300])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10)
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
shadedErrorBar(tv,wtplusm(:,1)+wtplusm(:,3)+wtplusm(:,5),wtplus_sspread,'lineProps', {'Color',[1,.42,.71]})
hold on
plot(tv,wtplusm(:,1)+wtplusm(:,3)+wtplusm(:,5),'LineWidth',1,'Color','k')
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
lgd = legend('WT RyhB','{\Delta RyhB}');
lgd.FontSize = 10;
lgd.Location = 'northwest';
title('WT RyhB','FontSize',16,'FontWeight','bold','Color','r','FontName','Arial')
xlabel('Time after Induction (min)','FontSize',14,'FontName','Arial')
ylabel('Copy Number','FontSize',14,'FontName','Arial')
set(gca,'YLim',[-10 300])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10)
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
% 
% figure(7)
% plot(linspace(0,1440,100),wtplusm(:,1),'LineWidth',3); 
% title('Free sRNA','FontSize',18);
% file1 = strcat([folderTitle,'sRNA']);
% set(gcf,'PaperPositionMode','auto');
% print(file1,'-dpng','-r0')
% 
% figure(8)
% plot(linspace(0,1440,100),wtplusm(:,2),'LineWidth',3); 
% title('Elongating mRNA','FontSize',18);
% file1 = strcat([folderTitle,'elong_mRNA']);
% set(gcf,'PaperPositionMode','auto');
% print(file1,'-dpng','-r0')
% 
% figure(9)
% plot(linspace(0,1440,100),wtplusm(:,3),'LineWidth',3); 
% title('Co-Transcriptionally Bound mRNA','FontSize',18);
% file1 = strcat([folderTitle,'cobound_mRNA']);
% set(gcf,'PaperPositionMode','auto');
% print(file1,'-dpng','-r0')
% 
% figure(10)
% plot(linspace(0,1440,100),wtplusm(:,4),'LineWidth',3); 
% title('Full mRNA','FontSize',18);
% file1 = strcat([folderTitle,'full_mRNA']);
% set(gcf,'PaperPositionMode','auto');
% print(file1,'-dpng','-r0')
% 
% figure(11)
% plot(linspace(0,1440,100),wtplusm(:,5),'LineWidth',3); 
% title('Post sRNA-mRNA Complex','FontSize',18);
% file1 = strcat([folderTitle,'s_mRNA']);
% set(gcf,'PaperPositionMode','auto');
% print(file1,'-dpng','-r0')
% 
% figure(12)
% plot(linspace(0,1440,100),wtplusm(:,6),'LineWidth',3); 
% title('Protein','FontSize',18);
% file1 = strcat([folderTitle,'protein']);
% set(gcf,'PaperPositionMode','auto');
% print(file1,'-dpng','-r0')
% 
% figure(13)
% plot(linspace(0,1440,100),rneplusm(:,1),'LineWidth',3); 
% title('rne701 Free sRNA','FontSize',18);
% file1 = strcat([folderTitle,'rne701_sRNA']);
% set(gcf,'PaperPositionMode','auto');
% print(file1,'-dpng','-r0')
% 
% figure(14)
% plot(linspace(0,1440,100),rneplusm(:,2),'LineWidth',3); 
% title('rne701 Elongating mRNA','FontSize',18);
% file1 = strcat([folderTitle,'rne701_elong_mRNA']);
% set(gcf,'PaperPositionMode','auto');
% print(file1,'-dpng','-r0')
% 
% figure(15)
% plot(linspace(0,1440,100),rneplusm(:,3),'LineWidth',3); 
% title('rne701 Co-Transcriptionally Bound mRNA','FontSize',18);
% file1 = strcat([folderTitle,'rne701_cobound_mRNA']);
% set(gcf,'PaperPositionMode','auto');
% print(file1,'-dpng','-r0')
% 
% figure(16)
% plot(linspace(0,1440,100),rneplusm(:,4),'LineWidth',3); 
% title('rne701 Full mRNA','FontSize',18);
% file1 = strcat([folderTitle,'rne701_full_mRNA']);
% set(gcf,'PaperPositionMode','auto');
% print(file1,'-dpng','-r0')
% 
% figure(17)
% plot(linspace(0,1440,100),rneplusm(:,5),'LineWidth',3); 
% title('rne701 Post sRNA-mRNA Complex','FontSize',18);
% file1 = strcat([folderTitle,'rne701_s_mRNA']);
% set(gcf,'PaperPositionMode','auto');
% print(file1,'-dpng','-r0')
% 
% figure(18)
% plot(linspace(0,1440,100),rneplusm(:,6),'LineWidth',3); 
% title('rne701 Protein','FontSize',18);
% file1 = strcat([folderTitle,'rne701_protein']);
% set(gcf,'PaperPositionMode','auto');
% print(file1,'-dpng','-r0')
% 
copy_array = [k_on,k_off,kx_s,b_e,b_ms,k_elon_prime,k_init_wt,k_init_rne,kx_wt_plus,kx_rne_plus,b_m,kx_wt_minus,kx_rne_minus,k_elon,k_nuc,b_nuc_wt,b_nuc_rne,k_off_nuc;
    k_on_sd,k_off_sd,kxs_sd,be_sd,b_ms_sd,k_elon_prime_sd,k_init_wt_sd,k_init_rne_sd,kxw_sd,kxr_sd,b_m_sd,kxw_sd,kxr_sd,0,k_nuc_sd,b_nuc_wt_sd,b_nuc_rne_sd,k_off_nuc_sd];

copy_table = array2table(copy_array);copy_table.Properties.VariableNames = {'k_on' 'k_off' 'kx_s' 'b_e' 'b_ms' 'k_elon_prime' 'k_init_wt' 'k_init_rne_' 'kx_wt_plus' 'kx_rne_plus' 'b_m' 'kx_wt_minus' 'kx_rne_minus' 'k_elon' 'k_nuc' 'b_nuc_wt' 'b_nuc_rne' 'k_off_nuc'};

copy_table_file = strcat([folderTitle,'parameters.csv']);

writetable(copy_table,copy_table_file);
