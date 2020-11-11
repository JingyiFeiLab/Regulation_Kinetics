folderTitle = '/Users/reyer/Data/SingleCellEpi/MR187/Fits4720/';

b_e = 5.77E-05; %std = 2.23E-3
be_sd = 0.00004325619532;
a_m_plus_wt = 4.13;% std =589.18; 
amp_wt_sd = 0.1710477516;
a_m_plus_rne = 2.537; %std = 796.64
amp_rne_sd = 0.5705690142;
a_s_wt = 0.2849; %a_s = 2.29E+03;
a_s_wt_sd = 0.04640542353;
a_s_rne = 0.3392; %a_s = 2.29E+03; 
a_s_rne_sd = 0.06959959291;
b_s_wt = 1.44E-03; 
b_s_wt_sd = 0.00005939696962;
b_s_rne = 0.0009734; 
b_s_rne_sd = 0.0001447447581;
b_m = 3.25E-03;
b_m_sd = 0.0000833446659;
k_on = 2.91E-03;
k_on_sd = 0.001377784601;
k_off = 11.78E+00;
k_off_sd = 0.04727714071;
kx_wt_minus = 2.318;
kx_wt_plus = 2.16; % kx_wt = 0.0040191;
kxw_sd = 0.074;
% kx_wt_plus=mode(kx_fit); %kx_wt_plus=0.0033;
kx_rne_minus = 1.4615;
kx_rne_plus = 1.389; %kx_rne = 0.0023639;
kxr_sd = 0.1404314067;
kx_s = 0.04296666667; %std = .15
kxs_sd = 0.04100662548;
b_ms_rne = 0.00333; % b_ms = 4.88E-03;
b_ms_wt = 0.00333; % b_ms = 4.88E-03;
b_ms_sd = 0.001143780282;
b_p = .00018;
a_m_minus_wt = 5.164; %std = 310.41
amw_sd = 0.1173797257;
a_m_minus_rne = 2.786; %std = 1164.1
amm_sd = 0.1979898987;


close all
tv = linspace(0, 1800);
tvp = linspace(0, 1440);
tvm = linspace(0,1620);
%tvm = linspace(time(1), max(time)-time(3));

timex = [0,60,180,360,720,1080,1440,1800];

%Protein
%mRNA
%sRNA

molecules_rne_plus1= [19492.43071,19015.52382,6280.162752,26398.27162,185786.0509,320896.4976,598652.1846,712535.0848;
                      16958.022,40849.6894,184819.7605,658452.4128,963770.2941,1048011.941,1168056.078,1055862.227;
                      825931.2773,1017981.839,1124410.438,1134054.485,1101884.722,1104395.974,1197324.715,1123533.284];

molecules_rne_plus2= [29571.2389,58368.68478,35320.27356,63683.48226,195297.6019,170485.2472,352867.5466,323987.3492;
                      26210.10641,151102.3531,337443.3026,694726.4592,886201.9069,735498.2732,827714.8431,666903.2247;
                      1032494.695,1166508.416,928444.2251,1046656.993,890193.9413,685346.8487,976810.8585,831907.646];

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

                   

moleculesrne_plus = [((mean([molecules_rne_plus1(3,1),molecules_rne_plus2(3,1)])-48376)/8096.5) mean(((([molecules_rne_plus1(2,1),molecules_rne_plus2(2,1)])-20022)/1125.7)+1) 0 mean([molecules_rne_plus1(1,4),molecules_rne_plus2(1,4)])];
moleculeswt_plus = [((mean([molecules_wt_plus1(3,1),molecules_wt_plus2(3,1),molecules_wt_plus3(3,1)])-48376)/8096.5) mean(((([molecules_wt_plus1(2,1),molecules_wt_plus2(2,1),molecules_wt_plus3(2,1)])-20022)/1125.7)+1) 0 mean([molecules_wt_plus1(1,4),molecules_wt_plus2(1,4),molecules_wt_plus3(1,4)])];
moleculesrne_minus = [((mean([molecules_rne_minus1(3,1),molecules_rne_minus2(3,1)])-48376)/8096.5) mean(((([molecules_rne_minus1(2,1),molecules_rne_minus2(2,1)])-20022)/1125.7)+1) 0 mean([molecules_rne_minus1(1,4),molecules_rne_minus2(1,4)])];
moleculeswt_minus = [((mean([molecules_wt_minus1(3,1),molecules_wt_minus2(3,1)])-48376)/8096.5) mean(((([molecules_wt_minus1(2,1),molecules_wt_minus2(2,1)])-20022)/1125.7)+1) 0 mean([molecules_wt_minus1(1,4),molecules_wt_minus2(1,4)])];

std_rne_plus = [std((([molecules_rne_plus1(3,1),molecules_rne_plus2(3,1)])-48376)/8096.5) std(((([molecules_rne_plus1(2,1),molecules_rne_plus2(2,1)])-20022)/1125.7)+1) 0 std([molecules_rne_plus1(1,4),molecules_rne_plus2(1,4)])];
std_wt_plus = [std((([molecules_wt_plus1(3,1),molecules_wt_plus2(3,1),molecules_wt_plus3(3,1)])-48376)/8096.5) std(((([molecules_wt_plus1(2,1),molecules_wt_plus2(2,1),molecules_wt_plus3(2,1)])-20022)/1125.7)+1) 0 std([molecules_wt_plus1(1,4),molecules_wt_plus2(1,4),molecules_wt_plus3(1,4)])];
std_rne_minus = [std((([molecules_rne_minus1(3,1),molecules_rne_minus2(3,1)])-48376)/8096.5) std(((([molecules_rne_minus1(2,1),molecules_rne_minus2(2,1)])-20022)/1125.7)+1) 0 std([molecules_rne_minus1(1,4),molecules_rne_minus2(1,4)])];
std_wt_minus = [std((([molecules_wt_minus1(3,1),molecules_wt_minus2(3,1)])-48376)/8096.5) std(((([molecules_wt_minus1(2,1),molecules_wt_minus2(2,1)])-20022)/1125.7)+1) 0 std([molecules_wt_minus1(1,4),molecules_wt_minus2(1,4)])];

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
    molecules_rne_minus(1,i) = mean([molecules_rne_minus1(1,i),molecules_rne_minus2(1,i)]); 
    molecules_wt_plus(1,i) = mean([molecules_wt_plus1(1,i),molecules_wt_plus2(1,i),molecules_wt_plus3(1,i)]);
    molecules_wt_minus(1,i) = mean([molecules_wt_minus1(1,i),molecules_wt_minus2(1,i)]);
    
    molecules_rne_plus(2,i) = std([molecules_rne_plus1(1,i),molecules_rne_plus2(1,i)])/1.414;
    molecules_rne_minus(2,i) = std([molecules_rne_minus1(1,i),molecules_rne_minus2(1,i)])/1.414; 
    molecules_wt_plus(2,i) = std([molecules_wt_plus1(1,i),molecules_wt_plus2(1,i),molecules_wt_plus3(1,i)])/1.723;
    molecules_wt_minus(2,i) = std([molecules_wt_minus1(1,i),molecules_wt_minus2(1,i)])/1.414;
    
    molecules_rne_plus(3,i) = mean(((([molecules_rne_plus1(2,i),molecules_rne_plus2(2,i)])-20022)/1125.7)+1);
    molecules_rne_minus(3,i) = mean(((([molecules_rne_minus1(2,i),molecules_rne_minus2(2,i)])-20022)/1125.7)+1); 
    molecules_wt_plus(3,i) = mean(((([molecules_wt_plus1(2,i),molecules_wt_plus2(2,i),molecules_wt_plus3(2,i)])-20022)/1125.7)+1);
    molecules_wt_minus(3,i) = mean(((([molecules_wt_minus1(2,i),molecules_wt_minus2(2,i)])-20022)/1125.7)+1);
    
    molecules_rne_plus(4,i) = std(((([molecules_rne_plus1(2,i),molecules_rne_plus2(2,i)]-20022)/1125.7)+1))/1.414;
    molecules_rne_minus(4,i) = std(((([molecules_rne_minus1(2,i),molecules_rne_minus2(2,i)]-20022)/1125.7)+1))/1.414; 
    molecules_wt_plus(4,i) = std(((([molecules_wt_plus1(2,i),molecules_wt_plus2(2,i),molecules_wt_plus3(2,i)]-20022)/1125.7)+1))/1.723;
    molecules_wt_minus(4,i) = std(((([molecules_wt_minus1(2,i),molecules_wt_minus2(2,i)]-20022)/1125.7)+1))/1.414;
    
    molecules_rne_plus(5,i) = mean((([molecules_rne_plus1(3,i),molecules_rne_plus2(3,i)])-48376)/8096.5);
    molecules_rne_minus(5,i) = mean((([molecules_rne_minus1(3,i),molecules_rne_minus2(3,i)])-48376)/8096.5); 
    molecules_wt_plus(5,i) = mean((([molecules_wt_plus1(3,i),molecules_wt_plus2(3,i),molecules_wt_plus3(3,i)])-48376)/8096.5);
    molecules_wt_minus(5,i) = mean((([molecules_wt_minus1(3,i),molecules_wt_minus2(3,i)])-48376)/8096.5);
    
    molecules_rne_plus(6,i) = std((([molecules_rne_plus1(3,i),molecules_rne_plus2(3,i)]-48376)/8096.5))/1.414;
    molecules_rne_minus(6,i) = std((([molecules_rne_minus1(3,i),molecules_rne_minus2(3,i)]-48376)/8096.5))/1.414; 
    molecules_wt_plus(6,i) = std((([molecules_wt_plus1(3,i),molecules_wt_plus2(3,i),molecules_wt_plus3(3,i)]-48376)/8096.5))/1.723;
    molecules_wt_minus(6,i) = std((([molecules_wt_minus1(3,i),molecules_wt_minus2(3,i)]-48376)/8096.5))/1.414;
    
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
lgd = legend('WT SgrS','{\Delta SgrS}');
lgd.FontSize = 10;
lgd.FontWeight = 'bold';
lgd.Location = 'northwest';
title('rne701{\it manX-sfGFP}, fixed {\alpha_{m}}','FontSize',14,'FontWeight','bold','Color','g','FontName','Arial')
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
file1 = strcat([folderTitle,'rnemanXmRNA_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'rnemanXmRNA_fit.fig']);
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
lgd = legend('WT SgrS','{\Delta SgrS}');
lgd.FontSize = 10;
lgd.FontWeight = 'bold';
lgd.Location = 'northwest';
title('WT{\it manX-sfGFP}, fixed {\alpha_{m}}','FontSize',18,'FontWeight','bold','Color','g','FontName','Arial')
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
file1 = strcat([folderTitle,'WTmanXmRNA_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'WTmanXmRNA_fit.fig']);
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
lgd.FontWeight = 'bold';
lgd.Location = 'northwest';
title('rne701 sfGFP, fixed {\alpha_{m}}','FontSize',18,'FontWeight','bold','Color','b','FontName','Arial')
xlabel('Time after Induction (min)','FontSize',18,'FontWeight','bold')
ylabel('Fluorescence (A.U)','FontSize',18,'FontName','Arial')
set(gca,'YLim',[-10000 6E6])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10,'FontWeight','bold')
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[626,281,248,201])
file1 = strcat([folderTitle,'rnemanXProtein_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'rnemanXProtein_fit.fig']);
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
lgd.FontWeight = 'bold';
lgd.Location = 'northwest';
title('WT sfGFP, fixed {\alpha_{m}}','FontSize',14,'FontWeight','bold','Color','b')
xlabel('Time after Induction (min)','FontSize',12,'FontWeight','bold')
ylabel('Fluorescence (A.U)','FontSize',12,'FontWeight','bold')
set(gca,'YLim',[-10000 6E6])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10,'FontWeight','bold')
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[626,281,248,201])
file1 = strcat([folderTitle,'WTmanXProtein_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'WTmanXProtein_fit.fig']);
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
lgd = legend('WT SgrS','{\Delta SgrS}');
lgd.FontSize = 10;
lgd.FontWeight = 'bold';
lgd.Location = 'northwest';
title('rne701 SgrS, fixed {\alpha_{m}}','FontSize',18,'FontWeight','bold','Color','r','FontName','Arial')
xlabel('Time after Induction (min)','FontSize',18,'FontName','Arial')
ylabel('Copy Number','FontSize',18,'FontName','Arial')
set(gca,'YLim',[-10 3E6/8096.25])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10,'FontWeight','bold')
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[626,281,248,201])
file1 = strcat([folderTitle,'rnemanXsRNA_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'rnemanXsRNA_fit.fig']);
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
lgd = legend('WT SgrS','{\Delta SgrS}');
lgd.FontSize = 10;
lgd.FontWeight = 'bold';
lgd.Location = 'northwest';
title('WT SgrS, fixed {\alpha_{m}}','FontSize',18,'FontWeight','bold','Color','r')
xlabel('Time after Induction (min)','FontSize',18,'FontName','Arial')
ylabel('Copy Number','FontSize',18,'FontName','Arial')
set(gca,'YLim',[-10 3E6/8096.25])
set(gca,'XLim',[-1 1540])
set(gca,'XTick',[0,360,720,1080,1440])
set(gca,'XTickLabel',[0,6,12,18,24],'FontSize',10,'FontWeight','bold')
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[626,281,248,201])
file1 = strcat([folderTitle,'WTmanXsRNA_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'WTmanXsRNA_fit.fig']);
savefig(gcf,file1_fig)




