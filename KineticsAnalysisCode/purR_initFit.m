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
tv = linspace(0, 1800);
tvp = linspace(0, 1440);
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
                    
molecules_rne_minus2 = [23258.43053,33031.81622,20087.82158,15510.2511,71133.80532,228951.7932,395241.2526,598925.2416;
                        22684.94363,25823.41754,68785.51533,282973.9619,714643.3316,634445.9883,1332342.717,1450938.474;
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
    molecules_rne_minus(1,i) = mean([molecules_rne_minus1(1,i),molecules_rne_minus2(1,i)]); 
    molecules_wt_plus(1,i) = mean([molecules_wt_plus1(1,i),molecules_wt_plus2(1,i)]);
    molecules_wt_minus(1,i) = mean([molecules_wt_minus1(1,i),molecules_wt_minus2(1,i)]);
    
    molecules_rne_plus(2,i) = std([molecules_rne_plus1(1,i),molecules_rne_plus2(1,i)])/1.414;
    molecules_rne_minus(2,i) = std([molecules_rne_minus1(1,i),molecules_rne_minus2(1,i)])/1.414; 
    molecules_wt_plus(2,i) = std([molecules_wt_plus1(1,i),molecules_wt_plus2(1,i)])/1.414;
    molecules_wt_minus(2,i) = std([molecules_wt_minus1(1,i),molecules_wt_minus2(1,i)])/1.414;
    
    molecules_rne_plus(3,i) = mean(((([molecules_rne_plus1(2,i),molecules_rne_plus2(2,i)])+1.3257e+05)/6.0368e+03)+1);
    molecules_rne_minus(3,i) = mean(((([molecules_rne_minus1(2,i),molecules_rne_minus2(2,i)])+1.3257e+05)/6.0368e+03)+1); 
    molecules_wt_plus(3,i) = mean(((([molecules_wt_plus1(2,i),molecules_wt_plus2(2,i)])+1.3257e+05)/6.0368e+03)+1);
    molecules_wt_minus(3,i) = mean(((([molecules_wt_minus1(2,i),molecules_wt_minus2(2,i)])+1.3257e+05)/6.0368e+03)+1);
    
    molecules_rne_plus(4,i) = std(((([molecules_rne_plus1(2,i),molecules_rne_plus2(2,i)]+1.3257e+05)/6.0368e+03)+1))/1.414;
    molecules_rne_minus(4,i) = std(((([molecules_rne_minus1(2,i),molecules_rne_minus2(2,i)]+1.3257e+05)/6.0368e+03)+1))/1.414; 
    molecules_wt_plus(4,i) = std(((([molecules_wt_plus1(2,i),molecules_wt_plus2(2,i)]+1.3257e+05)/6.0368e+03)+1))/1.414;
    molecules_wt_minus(4,i) = std(((([molecules_wt_minus1(2,i),molecules_wt_minus2(2,i)]+1.3257e+05)/6.0368e+03)+1))/1.414;
    
    molecules_rne_plus(5,i) = mean((([molecules_rne_plus1(3,i),molecules_rne_plus2(3,i)])+5.2879e+04)/7.3302e+03);
    molecules_rne_minus(5,i) = mean((([molecules_rne_minus1(3,i),molecules_rne_minus2(3,i)])+5.2879e+04)/7.3302e+03); 
    molecules_wt_plus(5,i) = mean((([molecules_wt_plus1(3,i),molecules_wt_plus2(3,i)])+5.2879e+04)/7.3302e+03);
    molecules_wt_minus(5,i) = mean((([molecules_wt_minus1(3,i),molecules_wt_minus2(3,i)])+5.2879e+04)/7.3302e+03);
    
    molecules_rne_plus(6,i) = std((([molecules_rne_plus1(3,i),molecules_rne_plus2(3,i)]+5.2879e+04)/7.3302e+03))/1.414;
    molecules_rne_minus(6,i) = std((([molecules_rne_minus1(3,i),molecules_rne_minus2(3,i)]+5.2879e+04)/7.3302e+03))/1.414; 
    molecules_wt_plus(6,i) = std((([molecules_wt_plus1(3,i),molecules_wt_plus2(3,i)]+5.2879e+04)/7.3302e+03))/1.414;
    molecules_wt_minus(6,i) = std((([molecules_wt_minus1(3,i),molecules_wt_minus2(3,i)]+5.2879e+04)/7.3302e+03))/1.414;
    
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
lgd = legend('WT SgrS','{\Delta SgrS}');
lgd.FontSize = 10;
lgd.Location = 'northwest';
title('{\it rne701 purR-sfGFP}','FontSize',14,'Color',[0,.40,0],'FontName','Arial')
xlabel('Time after Induction (min)','FontSize',14,'FontName','Arial')
ylabel('Copy Number ','FontSize',14,'FontName','Arial')
set(gca,'YLim',[-10 200])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10)
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[626,281,248,201])
file1 = strcat([folderTitle,'rnepurRmRNA_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'rnepurRmRNA_fit.fig']);
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
lgd = legend('WT SgrS','{\Delta SgrS}');
lgd.FontSize = 10;
lgd.Location = 'northwest';
title('WT{\it purR-sfGFP}','FontSize',16,'Color',[0,.40,0],'FontName','Arial')
xlabel('Time after Induction (min)','FontSize',14,'FontName','Arial')
ylabel('Copy Number','FontSize',14,'FontName','Arial')
set(gca,'YLim',[-10 200])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10)
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[626,281,248,201])
file1 = strcat([folderTitle,'WTpurRmRNA_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'WTpurRmRNA_fit.fig']);
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
lgd = legend('WT SgrS','{\Delta SgrS}');
lgd.FontSize = 10;
lgd.Location = 'northwest';
title('{\it rne}701 sfGFP','FontSize',16,'Color','b','FontName','Arial')
xlabel('Time after Induction (min)','FontSize',14,'FontName','Arial')
ylabel('Fluorescence (A.U)','FontSize',14,'FontName','Arial')
set(gca,'YLim',[-10000 5E5])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10)
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[626,281,248,201])
file1 = strcat([folderTitle,'rnepurRProtein_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'rnepurRProtein_fit.fig']);
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
lgd = legend('WT SgrS','{\Delta SgrS}');
lgd.FontSize = 10;
lgd.Location = 'northwest';
title('WT sfGFP','FontSize',16,'FontName','Arial','Color','b')
xlabel('Time after Induction (min)','FontSize',14,'FontName','Arial')
ylabel('Fluorescence (A.U)','FontSize',14,'FontName','Arial')
set(gca,'YLim',[-10000 5E5])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10)
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[626,281,248,201])
file1 = strcat([folderTitle,'WTpurRProtein_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'WTpurRProtein_fit.fig']);
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
lgd = legend('WT SgrS','{\Delta SgrS}');
lgd.FontSize = 10;
lgd.Location = 'northwest';
title('rne701 SgrS','FontSize',16,'Color','r','FontName','Arial')
xlabel('Time after Induction (min)','FontSize',14,'FontName','Arial')
ylabel('Copy Number','FontSize',14,'FontName','Arial')
set(gca,'YLim',[-10 550])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10)
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[626,281,248,201])
file1 = strcat([folderTitle,'rnepurRsRNA_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'rnepurRsRNA_fit.fig']);
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
file1 = strcat([folderTitle,'WTpurRsRNA_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'WTpurRsRNA_fit.fig']);
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