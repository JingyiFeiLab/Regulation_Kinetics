close all
folderTitle = '/Users/reyer/Data/SingleCellEpi/MR156/premRNA8520/';

b_e = 1.33E-03; %std = 2.23E-3
be_sd = 0.0009905149974;

a_s_wt = 0.3646; %a_s = 2.29E+03;
a_s_wt_sd = 0.02771858582;
a_s_rne = 0.378825; %a_s = 2.29E+03; 
a_s_rne_sd = 0.01802764081;
b_s_wt = 1.46E-03; 
b_s_wt_sd = 0.0002935200249;
b_s_rne = 0.85E-03; 
b_s_rne_sd = 0.0001349089028;

b_m = 0.003229666667;
b_m_sd = 0.00007212489168;

k_on = 1.34E-03;
k_on_sd = 0.0006514260861;
k_off = 3.03E-01;
k_off_sd = 0.04203270155;

k_nuc = 1.34E-03;
k_nuc_sd = 0.0006514260861;
k_off_nuc = 3.03E-01;
k_off_nuc_sd = 0.04203270155;

b_ms = 0.004987666667; % b_ms = 4.88E-03;
b_ms_sd = 0.001467482561;
b_p = .00018;

b_nuc_wt = b_e + b_ms;
b_nuc_wt_sd = 0.001594820993;
b_nuc_rne = b_ms;
b_nuc_rne_sd = 0;

kx_wt_minus =  3.650233333;
kx_wt_plus = [15.5582,14.15,13.7,13.6]; % kx_wt = 0.0040191;
kxw_sd = 0.7665063622;

kx_rne_minus = 2.5366;
kx_rne_plus = 2.16075; %kx_rne = 0.0023639;
kxr_sd = 0.9100464274;

kx_s = 0.1395; %std = .15
kxs_sd = 0.03322472724;

k_init_wt = [0.10996,0.08896,0.08796,0.08696]; %std = 310.41
k_init_wt_sd = 0.009295339334;
k_init_rne = 0.2165; %std = 1164.1
k_init_rne_sd = 0.01909188309;

k_elon = 0.0643;
k_elon_prime = 0.0218;
k_elon_prime_sd = .0001;

n_plasmids = 20;

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
    
    
    
    molecules_wt_plus(3,i) = mean(((([wmp1(i),wmp2(i),wmp3(i)])+1.3257e+05)/6.0368e+03)+1);
    
     
    molecules_wt_plus(4,i) = std(((([wmp1(i),wmp2(i),wmp3(i)])+1.3257e+05)/6.0368e+03)+1)/2.25;
    
    
    molecules_wt_plus(5,i) = mean((([wsp1(i),wsp2(i),wsp3(i)])+5.2879e+04)/7.3302e+03);
    
    
     
    molecules_wt_plus(6,i) = std((([wsp1(i),wsp2(i),wsp3(i)])+5.2879e+04)/7.3302e+03)/1;
    
end

moleculeswt_plus1 = [0 0 molecules_wt_plus(5,1) molecules_wt_plus(3,1) 0 molecules_wt_plus(1,2)];
std_wt_plus1 = [ 0 0 molecules_wt_plus(6,1) molecules_wt_plus(4,1) 0 molecules_wt_plus(2,2)];


wtplusp1 = zeros(100,120);
wtplusm1 = zeros(100,120);
wtpluss1 = zeros(100,120);

wtplusp2 = zeros(100,120);
wtplusm2 = zeros(100,120);
wtpluss2 = zeros(100,120);

wtplusp3 = zeros(100,120);
wtplusm3 = zeros(100,120);
wtpluss3 = zeros(100,120);

wtplusp4 = zeros(100,120);
wtplusm4 = zeros(100,120);
wtpluss4 = zeros(100,120);


wtplus_pspread = zeros(100,1);
wtplus_mspread = zeros(100,1);
wtplus_sspread = zeros(100,1);


wtplusp1(:,1:6) = bestInit_fit(tvp1,a_s_wt,b_s_wt,k_init_wt(1),k_elon,k_elon_prime,b_m,k_on,k_off,kx_wt_plus(1),kx_s*kx_wt_plus(1),b_ms,b_e,b_p,k_nuc,b_nuc_wt,n_plasmids,k_off_nuc,moleculeswt_plus1);
wtplusm1(:,1:6) = bestInit_fit(tv1,a_s_wt,b_s_wt,k_init_wt(1),k_elon,k_elon_prime,b_m,k_on,k_off,kx_wt_plus(1),kx_s*kx_wt_plus(1),b_ms,b_e,b_p,k_nuc,b_nuc_wt,n_plasmids,k_off_nuc,moleculeswt_plus1);


moleculeswt_plus2 = [wtplusm1(100,1) wtplusm1(100,2) wtplusm1(100,3) wtplusm1(100,4) wtplusm1(100,5) wtplusp1(100,6)];

wtplusp2(:,1:6) = bestInit_fit(tvp2,a_s_wt,b_s_wt,k_init_wt(2),k_elon,k_elon_prime,b_m,k_on,k_off,kx_wt_plus(2),kx_s*kx_wt_plus(2),b_ms,b_e,b_p,k_nuc,b_nuc_wt,n_plasmids,k_off_nuc,moleculeswt_plus2);
wtplusm2(:,1:6) = bestInit_fit(tv2,a_s_wt,b_s_wt,k_init_wt(2),k_elon,k_elon_prime,b_m,k_on,k_off,kx_wt_plus(2),kx_s*kx_wt_plus(2),b_ms,b_e,b_p,k_nuc,b_nuc_wt,n_plasmids,k_off_nuc,moleculeswt_plus2);

moleculeswt_plus3 = [wtplusm2(100,1) wtplusm2(100,2) wtplusm2(100,3) wtplusm2(100,4) wtplusm2(100,5) wtplusp2(100,6)];

wtplusp3(:,1:6) = bestInit_fit(tvp3,a_s_wt,b_s_wt,k_init_wt(3),k_elon,k_elon_prime,b_m,k_on,k_off,kx_wt_plus(3),kx_s*kx_wt_plus(3),b_ms,b_e,b_p,k_nuc,b_nuc_wt,n_plasmids,k_off_nuc,moleculeswt_plus3);
wtplusm3(:,1:6) = bestInit_fit(tv3,a_s_wt,b_s_wt,k_init_wt(3),k_elon,k_elon_prime,b_m,k_on,k_off,kx_wt_plus(3),kx_s*kx_wt_plus(3),b_ms,b_e,b_p,k_nuc,b_nuc_wt,n_plasmids,k_off_nuc,moleculeswt_plus3);

moleculeswt_plus4 = [wtplusm3(100,1) wtplusm3(100,2) wtplusm3(100,3) wtplusm3(100,4) wtplusm3(100,5) wtplusp3(100,6)];

wtplusp4(:,1:6) = bestInit_fit(tvp4,a_s_wt,b_s_wt,k_init_wt(4),k_elon,k_elon_prime,b_m,k_on,k_off,kx_wt_plus(4),kx_s*kx_wt_plus(4),b_ms,b_e,b_p,k_nuc,b_nuc_wt,n_plasmids,k_off_nuc,moleculeswt_plus4);
wtplusm4(:,1:6) = bestInit_fit(tv4,a_s_wt,b_s_wt,k_init_wt(4),k_elon,k_elon_prime,b_m,k_on,k_off,kx_wt_plus(4),kx_s*kx_wt_plus(4),b_ms,b_e,b_p,k_nuc,b_nuc_wt,n_plasmids,k_off_nuc,moleculeswt_plus4);



i_rand = 7:12;
i_rand_check = [];
for i = 1:19
    
    a31 = bestInit_fit(tvp1,normrnd(a_s_wt,a_s_wt_sd),normrnd(b_s_wt,b_s_wt_sd),normrnd(k_init_wt(1),k_init_wt_sd),k_elon,normrnd(k_elon_prime,k_elon_prime_sd),normrnd(b_m,b_m_sd),normrnd(k_on,k_on_sd),normrnd(k_off,k_off_sd),normrnd(kx_wt_plus(1),kxr_sd),normrnd(kx_s,kxs_sd)*kx_wt_plus(1),normrnd(b_ms,b_ms_sd),normrnd(b_e,be_sd),b_p,normrnd(k_nuc,k_nuc_sd),normrnd(b_nuc_wt,b_nuc_wt_sd),n_plasmids,normrnd(k_off_nuc,k_off_nuc_sd),normrnd(moleculeswt_plus1,std_wt_plus1));
    if length(a31(:,1)) == 100
        wtplusp1(:,i_rand) = a31;
    else
        i_rand_check = [i_rand_check i_rand2];
    end
    
    moleculeswt_plus2 = [a31(100,1) a31(100,2) a31(100,3) a31(100,4) a31(100,5) a31(100,6)];
    
    a32 = bestInit_fit(tvp2,normrnd(a_s_wt,a_s_wt_sd),normrnd(b_s_wt,b_s_wt_sd),normrnd(k_init_wt(2),k_init_wt_sd),k_elon,normrnd(k_elon_prime,k_elon_prime_sd),normrnd(b_m,b_m_sd),normrnd(k_on,k_on_sd),normrnd(k_off,k_off_sd),normrnd(kx_wt_plus(2),kxr_sd),normrnd(kx_s,kxs_sd)*kx_wt_plus(2),normrnd(b_ms,b_ms_sd),normrnd(b_e,be_sd),b_p,normrnd(k_nuc,k_nuc_sd),normrnd(b_nuc_wt,b_nuc_wt_sd),n_plasmids,normrnd(k_off_nuc,k_off_nuc_sd),normrnd(moleculeswt_plus2,std_wt_plus1));
    if length(a32(:,1)) == 100
        wtplusp2(:,i_rand) = a32;
    else
        i_rand_check = [i_rand_check i_rand];
        continue
    end
    
    moleculeswt_plus3 = [a32(100,1) a32(100,2) a32(100,3) a32(100,4) a32(100,5) a32(100,6)];
    
    a33 = bestInit_fit(tvp3,normrnd(a_s_wt,a_s_wt_sd),normrnd(b_s_wt,b_s_wt_sd),normrnd(k_init_wt(3),k_init_wt_sd),k_elon,normrnd(k_elon_prime,k_elon_prime_sd),normrnd(b_m,b_m_sd),normrnd(k_on,k_on_sd),normrnd(k_off,k_off_sd),normrnd(kx_wt_plus(3),kxr_sd),normrnd(kx_s,kxs_sd)*kx_wt_plus(3),normrnd(b_ms,b_ms_sd),normrnd(b_e,be_sd),b_p,normrnd(k_nuc,k_nuc_sd),normrnd(b_nuc_wt,b_nuc_wt_sd),n_plasmids,normrnd(k_off_nuc,k_off_nuc_sd),normrnd(moleculeswt_plus3,std_wt_plus1));
    if length(a33(:,1)) == 100
        wtplusp3(:,i_rand) = a33;
    else
        i_rand_check = [i_rand_check i_rand];
        continue
    end
    
    moleculeswt_plus4 = [a33(100,1) a33(100,2) a33(100,3) a33(100,4) a33(100,5) a33(100,6)];
    
    a34 = bestInit_fit(tvp4,normrnd(a_s_wt,a_s_wt_sd),normrnd(b_s_wt,b_s_wt_sd),normrnd(k_init_wt(4),k_init_wt_sd),k_elon,normrnd(k_elon_prime,k_elon_prime_sd),normrnd(b_m,b_m_sd),normrnd(k_on,k_on_sd),normrnd(k_off,k_off_sd),normrnd(kx_wt_plus(4),kxr_sd),normrnd(kx_s,kxs_sd)*kx_wt_plus(4),normrnd(b_ms,b_ms_sd),normrnd(b_e,be_sd),b_p,normrnd(k_nuc,k_nuc_sd),normrnd(b_nuc_wt,b_nuc_wt_sd),n_plasmids,normrnd(k_off_nuc,k_off_nuc_sd),normrnd(moleculeswt_plus4,std_wt_plus1));
    if length(a34(:,1)) == 100
        wtplusp4(:,i_rand) = a34;
    else
        i_rand_check = [i_rand_check i_rand];
        continue
    end
    
    a41 = bestInit_fit(tv1,normrnd(a_s_wt,a_s_wt_sd),normrnd(b_s_wt,b_s_wt_sd),normrnd(k_init_wt(1),k_init_wt_sd),k_elon,normrnd(k_elon_prime,k_elon_prime_sd),normrnd(b_m,b_m_sd),normrnd(k_on,k_on_sd),normrnd(k_off,k_off_sd),normrnd(kx_wt_plus(1),kxr_sd),normrnd(kx_s,kxs_sd)*kx_wt_plus(1),normrnd(b_ms,b_ms_sd),normrnd(b_e,be_sd),b_p,normrnd(k_nuc,k_nuc_sd),normrnd(b_nuc_wt,b_nuc_wt_sd),n_plasmids,normrnd(k_off_nuc,k_off_nuc_sd),normrnd(moleculeswt_plus1,std_wt_plus1));
    if length(a41(:,1)) == 100
        wtplusm1(:,i_rand) = a41;
    else
        i_rand_check = [i_rand_check i_rand];
        continue
    end
    
    moleculeswt_plus2 = [a41(100,1) a41(100,2) a41(100,3) a41(100,4) a41(100,5) a41(100,6)];
    
    a42 = bestInit_fit(tv2,normrnd(a_s_wt,a_s_wt_sd),normrnd(b_s_wt,b_s_wt_sd),normrnd(k_init_wt(2),k_init_wt_sd),k_elon,normrnd(k_elon_prime,k_elon_prime_sd),normrnd(b_m,b_m_sd),normrnd(k_on,k_on_sd),normrnd(k_off,k_off_sd),normrnd(kx_wt_plus(2),kxr_sd),normrnd(kx_s,kxs_sd)*kx_wt_plus(2),normrnd(b_ms,b_ms_sd),normrnd(b_e,be_sd),b_p,normrnd(k_nuc,k_nuc_sd),normrnd(b_nuc_wt,b_nuc_wt_sd),n_plasmids,normrnd(k_off_nuc,k_off_nuc_sd),normrnd(moleculeswt_plus2,std_wt_plus1));
    if length(a42(:,1)) == 100
        wtplusm2(:,i_rand) = a42;
    else
        i_rand_check = [i_rand_check i_rand];
        continue
    end
    
    moleculeswt_plus3 = [a42(100,1) a42(100,2) a42(100,3) a42(100,4) a42(100,5) a42(100,6)];
    
    a43 = bestInit_fit(tv3,normrnd(a_s_wt,a_s_wt_sd),normrnd(b_s_wt,b_s_wt_sd),normrnd(k_init_wt(3),k_init_wt_sd),k_elon,normrnd(k_elon_prime,k_elon_prime_sd),normrnd(b_m,b_m_sd),normrnd(k_on,k_on_sd),normrnd(k_off,k_off_sd),normrnd(kx_wt_plus(3),kxr_sd),normrnd(kx_s,kxs_sd)*kx_wt_plus(3),normrnd(b_ms,b_ms_sd),normrnd(b_e,be_sd),b_p,normrnd(k_nuc,k_nuc_sd),normrnd(b_nuc_wt,b_nuc_wt_sd),n_plasmids,normrnd(k_off_nuc,k_off_nuc_sd),normrnd(moleculeswt_plus3,std_wt_plus1));
    if length(a43(:,1)) == 100
        wtplusm3(:,i_rand) = a43;
    else
        i_rand_check = [i_rand_check i_rand];
        continue
    end
    
    moleculeswt_plus4 = [a43(100,1) a43(100,2) a43(100,3) a43(100,4) a43(100,5) a43(100,6)];
    
    a44 = bestInit_fit(tv4,normrnd(a_s_wt,a_s_wt_sd),normrnd(b_s_wt,b_s_wt_sd),normrnd(k_init_wt(4),k_init_wt_sd),k_elon,normrnd(k_elon_prime,k_elon_prime_sd),normrnd(b_m,b_m_sd),normrnd(k_on,k_on_sd),normrnd(k_off,k_off_sd),normrnd(kx_wt_plus(4),kxr_sd),normrnd(kx_s,kxs_sd)*kx_wt_plus(4),normrnd(b_ms,b_ms_sd),normrnd(b_e,be_sd),b_p,normrnd(k_nuc,k_nuc_sd),normrnd(b_nuc_wt,b_nuc_wt_sd),n_plasmids,normrnd(k_off_nuc,k_off_nuc_sd),normrnd(moleculeswt_plus4,std_wt_plus1));
    if length(a44(:,1)) == 100
        wtplusm4(:,i_rand) = a44;
        
    else
        i_rand_check = [i_rand_check i_rand];
        continue
    end
   
    
    
    i_rand = i_rand + 6;
    
end
    

wt_plus_protein_error = zeros(100,20);
wt_plus_mRNA_error = zeros(100,20);
wt_plus_sRNA_error = zeros(100,20);

for i = 1:400
    ij = 1;
    for j = 4:6:120
        
        if i<=100
            wt_plus_mRNA_error(i,ij) = wtplusm1(i,j) + wtplusm1(i,j+1);
            ij = ij + 1;
        elseif (i>100 & i<=200)
            wt_plus_mRNA_error(i,ij) = wtplusm2(i-100,j) + wtplusm2(i-100,j+1);
            ij = ij + 1;
        elseif (i>200 & i<=300)
            wt_plus_mRNA_error(i,ij) = wtplusm3(i-200,j) + wtplusm3(i-200,j+1);
            ij = ij + 1;
        elseif (i>300 & i<=400)
            wt_plus_mRNA_error(i,ij) = wtplusm4(i-300,j) + wtplusm4(i-300,j+1);
            ij = ij + 1;
        end
        
    end

    ik = 1;
    for k = 1:6:120
        
        
        
        if i<=100
            wt_plus_protein_error(i,k) = wtplusp1(i,k);
            ik = ik + 1;
        elseif (i>100 & i<=200)
            wt_plus_protein_error(i,k) = wtplusp2(i-100,k);
            ik = ik + 1;
        elseif (i>200 & i<=300)
            wt_plus_protein_error(i,k) = wtplusp3(i-200,k);
            ik = ik + 1;
        elseif (i>300 & i<=400)
            wt_plus_protein_error(i,k) = wtplusp4(i-300,k);
            ik = ik + 1;
        end
    end
    
    il = 1;
    for j = 1:6:120
        
        if i<=100
            wt_plus_sRNA_error(i,il) = wtplusm1(i,j) + wtplusm1(i,j+2) + wtplusm1(i,j+4);
            il = il + 1;
        elseif (i>100 & i<=200)
            wt_plus_sRNA_error(i,il) = wtplusm2(i-100,j) + wtplusm2(i-100,j+2) + wtplusm2(i-100,j+4);
            il = il + 1;
        elseif (i>200 & i<=300)
            wt_plus_sRNA_error(i,il) = wtplusm3(i-200,j) + wtplusm3(i-200,j+2) + wtplusm3(i-200,j+4);
            il = il + 1;
        elseif (i>300 & i<=400)
            wt_plus_sRNA_error(i,il) = wtplusm4(i-300,j) + wtplusm4(i-300,j+2) + wtplusm4(i-300,j+4);
            il = il + 1;
        end
        
    end
    
end

for i = 1:400
    
    if i == 101 || i == 201 || i == 301 
        wtplus_pspread(i) = wtplus_pspread(i-1);
        wtplus_mspread(i) = wtplus_mspread(i-1);
        wtplus_sspread(i) = wtplus_sspread(i-1);
        continue
    end
        
    
    percentiles = [];
    nonOutliers = [];
    percentiles = prctile(wt_plus_protein_error(i,:),[10,90]);
    outlierIndexes = wt_plus_protein_error(i,:) < percentiles(1) | wt_plus_protein_error(i,:) > percentiles(2);
    nonOutliers = wt_plus_protein_error(i,~outlierIndexes);
    
    wtplus_pspread(i) = std(nonOutliers);
    
    percentiles = [];
    nonOutliers = [];
    percentiles = prctile(wt_plus_mRNA_error(i,:),[10,90]);
    outlierIndexes = wt_plus_mRNA_error(i,:) < percentiles(1) | wt_plus_mRNA_error(i,:) > percentiles(2);
    nonOutliers = wt_plus_mRNA_error(i,~outlierIndexes);
    nonOutliers(nonOutliers>5E6) = [];
    
    wtplus_pmspread(i) = std(nonOutliers);
    
    percentiles = [];
    nonOutliers = [];
    percentiles = prctile(wt_plus_sRNA_error(i,:),[10,90]);
    outlierIndexes = wt_plus_sRNA_error(i,:) < percentiles(1) | wt_plus_sRNA_error(i,:) > percentiles(2);
    nonOutliers = wt_plus_sRNA_error(i,~outlierIndexes);
    nonOutliers(nonOutliers>5E6) = [];
    
    wtplus_sspread(i) = std(nonOutliers);
    
    
end  

wtplus_pspread = sgolayfilt(wtplus_pspread,11,51);
wtplus_mspread = sgolayfilt(wtplus_mspread,11,51);
wtplus_sspread = sgolayfilt(wtplus_sspread,11,51);

figure(2)
e3 = errorbar(timex,molecules_wt_plus(3,:),molecules_wt_plus(4,:),'s');
e3.MarkerFaceColor = [0,.3,0];
e3.MarkerSize = 5;
e3.Color = [0,.3,0];
e3.LineWidth = 2;
hold on
plot([tv1';(tv2+360)';(tv3+720)';(tv4+1200)'],[wtplusm1(:,4)+wtplusm1(:,5);wtplusm2(:,4)+wtplusm2(:,5);wtplusm3(:,4)+wtplusm3(:,5);wtplusm4(:,4)+wtplusm4(:,5)],'LineWidth',1,'Color','k')

% 
% legend('WT +SgrS','WT -SgrS')
shadedErrorBar([tv1';(tv2+360)';(tv3+720)';(tv4+1200)'],[wtplusm1(:,4)+wtplusm1(:,5);wtplusm2(:,4)+wtplusm2(:,5);wtplusm3(:,4)+wtplusm3(:,5);wtplusm4(:,4)+wtplusm4(:,5)],wtplus_mspread,'lineProps', {'Color',[0,.40,0]})
hold on
plot([tv1';(tv2+360)';(tv3+720)';(tv4+1200)'],[wtplusm1(:,4)+wtplusm1(:,5);wtplusm2(:,4)+wtplusm2(:,5);wtplusm3(:,4)+wtplusm3(:,5);wtplusm4(:,4)+wtplusm4(:,5)],'LineWidth',1,'Color','k')
hold on

e3 = errorbar(timex,molecules_wt_plus(3,:),molecules_wt_plus(4,:),'s');
e3.MarkerFaceColor = [0,.3,0];
e3.MarkerSize = 5;
e3.Color = [0,.3,0];
e3.LineWidth = 2;
% 
title('Pre-Induced mRNA, {\it ptsG-sfGFP} (post)','FontSize',14,'FontWeight','bold','Color','g')
xlabel('Time after Induction (min)','FontSize',12,'FontWeight','bold')
ylabel('Copy Number','FontSize',18,'FontWeight','bold')
set(gca,'YLim',[0 6.5E6/8125.7])
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
plot([(tvp1+360)';(tvp2+360+360)';(tvp3+720+360)';(tvp4+1200+360)'],[wtplusp1(:,6);wtplusp2(:,6);wtplusp3(:,6);wtplusp4(:,6)],'LineWidth',1,'Color','k')

% 
% legend('WT +SgrS','WT -SgrS')
shadedErrorBar([(tvp1+360)';(tvp2+360+360)';(tvp3+720+360)';(tvp4+1200+360)'],[wtplusp1(:,6);wtplusp2(:,6);wtplusp3(:,6);wtplusp4(:,6)],wtplus_pspread,'lineProps', {'Color',[.05,.05,.38]})
hold on
plot([(tvp1+360)';(tvp2+360+360)';(tvp3+720+360)';(tvp4+1200+360)'],[wtplusp1(:,6);wtplusp2(:,6);wtplusp3(:,6);wtplusp4(:,6)],'LineWidth',1,'Color','k')
hold on

e7 = errorbar(timex,molecules_wt_plus(1,:),molecules_wt_plus(2,:),'s');
e7.MarkerFaceColor = [.05,.05,.53];
e7.MarkerSize = 5;
e7.Color = [.05,.05,.53];
e7.LineWidth = 2;
% 
title('Pre-Induced mRNA, sfGFP (post)','FontSize',14,'FontWeight','bold','Color','b')
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
plot([tv1';(tv2+360)';(tv3+720)';(tv4+1200)'],[wtplusm1(:,1)+wtplusm1(:,3)+wtplusm1(:,5);wtplusm2(:,1)+wtplusm2(:,3)+wtplusm2(:,5);wtplusm3(:,1)+wtplusm3(:,3)+wtplusm3(:,5);wtplusm4(:,1)+wtplusm4(:,3)+wtplusm4(:,5)],'LineWidth',1,'Color','k')

% 
hold on
legend('WT +SgrS')
shadedErrorBar([tv1';(tv2+360)';(tv3+720)';(tv4+1200)'],[wtplusm1(:,1)+wtplusm1(:,3)+wtplusm1(:,5);wtplusm2(:,1)+wtplusm2(:,3)+wtplusm2(:,5);wtplusm3(:,1)+wtplusm3(:,3)+wtplusm3(:,5);wtplusm4(:,1)+wtplusm4(:,3)+wtplusm4(:,5)],wtplus_sspread,'lineProps', {'Color',[1,.42,.71]})
hold on
plot([tv1';(tv2+360)';(tv3+720)';(tv4+1200)'],[wtplusm1(:,1)+wtplusm1(:,3)+wtplusm1(:,5);wtplusm2(:,1)+wtplusm2(:,3)+wtplusm2(:,5);wtplusm3(:,1)+wtplusm3(:,3)+wtplusm3(:,5);wtplusm4(:,1)+wtplusm4(:,3)+wtplusm4(:,5)],'LineWidth',1,'Color','k')

hold on
e11 = errorbar(timex,molecules_wt_plus(5,:),molecules_wt_plus(6,:),'s');
e11.MarkerFaceColor = [1,.75,.79];
e11.MarkerSize = 5;
e11.Color = [1,.75,.79];
e3.LineWidth = 2;
hold on

title('Pre-Induced mRNA, SgrS (post)','FontSize',14,'FontWeight','bold','Color','r')
xlabel('Time after Induction (min)','FontSize',12,'FontWeight','bold')
ylabel('Copy Number','FontSize',18,'FontWeight','bold')
set(gca,'YLim',[0 1E6/8096.25])
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