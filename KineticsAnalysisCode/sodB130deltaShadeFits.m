folderTitle = '/Users/reyer/Data/SingleCellEpi/MR241/FixedP51820/';

b_e = 2.17E-04; %std = 2.23E-3
be_sd = 0.005755776157;
a_m_plus_wt = 8.771;% std =589.18; 
amp_wt_sd = 0.3419571123;
a_m_plus_rne = 6.683333333; %std = 796.64
amp_rne_sd = 0.109853994;

a_s_wt = 0.2391; %a_s = 2.29E+03;
a_s_wt_sd = 0.01852608971;
a_s_rne = 0.1683; %a_s = 2.29E+03; 
a_s_rne_sd = 0.01018233765;
b_s_wt = 2.98E-03; 
b_s_wt_sd = 0.000007071067812;
b_s_rne = 1.65E-03; 
b_s_rne_sd = 0.00006687550623;

b_m = 0.0049828;
b_m_sd = 0.0001734848697;
k_on = 6.27E-05;
k_on_sd = 0.000004081457064;
k_off = 6.40E+00;
k_off_sd = 4.19419086;
kx_wt_minus = 3.578333333;
kx_wt_plus = 2.2945; % kx_wt = 0.0040191;
kxw_sd = 1.229775996;
% kx_wt_plus=mode(kx_fit); %kx_wt_plus=0.0033;
kx_rne_minus = 1.81885;
kx_rne_plus = 1.81885; %kx_rne = 0.0023639;
kxr_sd = 0.2420538576;
kx_s = 0.4976666667; %std = .15
kxs_sd = 0.2698100566;
b_ms_rne = 0.01188833333; % b_ms = 4.88E-03;
b_ms_wt = 0.01188833333; % b_ms = 4.88E-03;
b_ms_sd = 0.006834095161;
b_p = .00018;
a_m_minus_wt = 8.877333333; %std = 310.41
amw_sd = 1.056258505;
a_m_minus_rne = 7.6579; %std = 1164.1
amm_sd = 1.284592039;


close all
tv = linspace(0, 1800);
tvp = linspace(0, 1440);
tvm = linspace(0,1620);
%tvm = linspace(time(1), max(time)-time(3));

timex = [0,60,180,360,720,1080,1440,1800];

%Protein
%mRNA
%sRNA

molecules_rne_plus1= [312459.207,388269.1916,278414.2402,348870.4773,1148738.288,1995535.445,2962595.947,235154.4538;
                      21022.29923,68618.7328,220829.1215,906086.2146,1485455.014,1601980.683,1879151.126,500606.5183;
                      301214.1291,442893.0082,294316.6629,300909.6674,290941.0776,338182.4746,318599.624,294238.0149];

molecules_rne_plus2= [296417.9704,327376.9363,413570.8962,380526.0993,713314.1639,1220653.148,1838376.74,262801.7701;
                      23448.8567,22575.82339,451505.272,1017287.758,1230098.949,1235860.031,1572049.379,532382.8749;
                      422493.8283,536884.2943,528911.3977,366595.5553,480627.8756,315041.4438,386024.6813,276181.5914];
                  
molecules_rne_plus3= [90162.12657,61207.16208,65190.31855,158122.6599,223439.5954,452600.6235,651476.1381,262801.7701;
                      30507.24666,31862.90128,132319.6045,256098.6477,298893.0509,294603.1303,312200.7337,532382.8749;
                      243615.076,203211.127,227068.3776,217371.7269,188423.0321,144095.6365,168201.7129,276181.5914];
                  
molecules_wt_plus1 = [928751.0435,945539.7993,994458.6726,1146200.253,1469220.281,2045016.476,3012954.447,3918573.443;
                      66087.63366,609091.0265,1445399.595,1480062.129,1826466.807,1750047.466,2274210.835,1160142.451;
                      386293.9825,382661.9625,341246.2669,222210.0288,328059.835,301982.5563,381734.8808,147913.77];
                  
molecules_wt_plus2 = [895378.3744,836236.4094,1040130.802,1042638.447,1461809.728,1963211.599,2919667.011,4782463.444;
                      33857.67064,631980.0823,1530505.495,1825334.144,1824538.821,1854408.497,2260601.543,1004738.614;
                      437730.8487,429529.892,410606.4183,418690.8348,399038.0563,377294.9272,413253.0065,222173.9802];

molecules_wt_plus3 = [146411.234,314756.5486,90745.19395,130325.5833,421470.4854,781590.5732,1144718.148,4782463.444;
                      20937.9833,503496.4997,1308419.514,1743107.298,1855637.381,1682596.738,2018330.604,1004738.614;
                      401611.3152,391506.2667,365125.276,371091.6192,333748.0932,321565.9084,305801.9997,1004738.614];                  

molecules_rne_minus1 = [80545.57201,137147.1782,104911.1979,160131.3138,750325.6481,1415555.028,2005143.723,2628169.692;
                        21121.39408,45329.91755,296115.6317,931845.8213,1838717.57,1726826.007,1585517.156,2496506.387;
                        0,0,0,0,0,0,0,0];
                    
molecules_rne_minus2 = [68198.13104,69966.62482,98852.00677,271534.1983,1040700.168,1835223.601,2364348.863,1937922.479;
                        24076.11943,16762.53937,341881.1211,979351.1708,1755030.777,1709674.98,1970561.677,1547179.084;
                        0,0,0,0,0,0,0,0]; 

molecules_rne_minus4 = [854520.5634,22678.87053,22586.06235,64972.34671,641637.0878,1335501.946,2445090.21,1937922.479;
                        61982.99973,70176.7428,482825.6819,1268276.623,2220302.773,2302748.798,2383273.072,1547179.084;
                        0,0,0,0,0,0,0,0];
                    
molecules_rne_minus3 = [132211.7047,215780.8637,227198.1701,243387.2755,618467.483,1549829.57,2667134.664,2839824.559;
                        42454.28066,62260.73866,429131.978,1029525.732,1673209.165,2220903.048,2284088.89,1938930.582;
                        0,0,0,0,0,0,0,0];                    
                    
                   
molecules_wt_minus1 = [1191010.547,1452117.131,1280905.299,1747603.901,3789115.753,5525198.702,7704268.795,11324996.58;
                       46594.41678,1213907.757,634894.6381,1987886.549,2224138.887,2610187.406,2825902.349,3451885.389;
                       0,0,0,0,0,0,0,0]; 
                   
molecules_wt_minus2 = [202451.6866,147376.4006,161676.2002,271224.2695,1503055.708,2639928.363,4148614.04,5864496.348;
                       27084.46433,517504.3554,1776995.851,1818249.57,2513161.418,1876589.261,1935856.061,1484099.069;
                       0,0,0,0,0,0,0,0];                    
                                     
molecules_wt_minus3 = [242656.4358,159771.1847,230764.5679,577911.8746,3037260.077,5371558.552,7516932.304,8927913.632;
                       44724.45028,600252.8811,1601328.952,1869381.377,1520591.167,1707582.68,2188694.773,2506915.665;
                       0,0,0,0,0,0,0,0];                   

molecules_wt_minus4 = [202451.6866,147376.4006,161676.2002,271224.2695,1503055.708,2639928.363,4148614.04,5864496.348;
                       27084.46433,517504.3554,1776995.851,1818249.57,2513161.418,1876589.261,1935856.061,1484099.069;
                       0,0,0,0,0,0,0,0];                   
                                   
                   
moleculesrne_plus = [(9/4)*mean((([molecules_rne_plus1(3,1),molecules_rne_plus2(3,1)])-48376)/8096.5) mean(((([molecules_rne_plus1(2,1),molecules_rne_plus2(2,1)])-20022)/1125.7)+1) 0  mean([molecules_rne_plus1(1,4),molecules_rne_plus2(1,4)])];
moleculeswt_plus = [(9/4)*mean((([molecules_wt_plus1(3,1),molecules_wt_plus2(3,1)])-48376)/8096.5) mean(((([molecules_wt_plus1(2,1),molecules_wt_plus2(2,1)])-20022)/1125.7)+1) 0  mean([molecules_wt_plus1(1,4),molecules_wt_plus2(1,4)])];
moleculesrne_minus = [(9/4)*mean((([molecules_rne_minus1(3,1),molecules_rne_minus2(3,1),molecules_rne_minus3(3,1),molecules_rne_minus4(3,1)])-48376)/8096.5) mean(((([molecules_rne_minus1(2,1),molecules_rne_minus2(2,1),molecules_rne_minus3(2,1),molecules_rne_minus4(2,1)])-20022)/1125.7)+1) 0 mean([molecules_rne_minus1(1,4),molecules_rne_minus2(1,4),molecules_rne_minus3(1,4),molecules_rne_minus4(1,4)])];
moleculeswt_minus = [(9/4)*mean((([molecules_wt_minus1(3,1),molecules_wt_minus2(3,1),molecules_wt_minus3(3,1),molecules_wt_minus4(3,1)])-48376)/8096.5) mean(((([molecules_wt_minus1(2,1),molecules_wt_minus2(2,1),molecules_wt_minus3(2,1),molecules_wt_minus4(2,1)])-20022)/1125.7)+1) 0 mean([molecules_wt_minus1(1,4),molecules_wt_minus2(1,4),molecules_wt_minus3(1,4),molecules_wt_minus4(1,4)])];

std_rne_plus = [(9/4)*std((([molecules_rne_plus1(3,1),molecules_rne_plus2(3,1)])-48376)/8096.5) std(((([molecules_rne_plus1(2,1),molecules_rne_plus2(2,1)])-20022)/1125.7)+1) 0  std([molecules_rne_plus1(1,4),molecules_rne_plus2(1,4)])];
std_wt_plus = [(9/4)*std((([molecules_wt_plus1(3,1),molecules_wt_plus2(3,1)])-48376)/8096.5) std(((([molecules_wt_plus1(2,1),molecules_wt_plus2(2,1)])-20022)/1125.7)+1) 0  std([molecules_wt_plus1(1,4),molecules_wt_plus2(1,4)])];
std_rne_minus = [(9/4)*std((([molecules_rne_minus1(3,1),molecules_rne_minus2(3,1),molecules_rne_minus3(3,1),molecules_rne_minus4(3,1)])-48376)/8096.5) std(((([molecules_rne_minus1(2,1),molecules_rne_minus2(2,1),molecules_rne_minus3(2,1),molecules_rne_minus4(2,1)])-20022)/1125.7)+1) 0 std([molecules_rne_minus1(1,4),molecules_rne_minus2(1,4),molecules_rne_minus3(1,4),molecules_rne_minus4(1,4)])];
std_wt_minus = [(9/4)*std((([molecules_wt_minus1(3,1),molecules_wt_minus2(3,1),molecules_wt_minus3(3,1),molecules_wt_minus4(3,1)])-48376)/8096.5) std(((([molecules_wt_minus1(2,1),molecules_wt_minus2(2,1),molecules_wt_minus3(2,1),molecules_wt_minus4(2,1)])-20022)/1125.7)+1) 0 std([molecules_wt_minus1(1,4),molecules_wt_minus2(1,4),molecules_wt_minus3(1,4),molecules_wt_minus4(1,4)])];


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
    molecules_rne_minus(1,i) = mean([molecules_rne_minus1(1,i),molecules_rne_minus2(1,i),molecules_rne_minus3(1,i),molecules_rne_minus4(1,i)]); 
    molecules_wt_plus(1,i) = mean([molecules_wt_plus1(1,i),molecules_wt_plus2(1,i)]);
    molecules_wt_minus(1,i) = mean([molecules_wt_minus1(1,i),molecules_wt_minus2(1,i),molecules_wt_minus3(1,i),molecules_wt_minus4(1,i)]);
    
    molecules_rne_plus(2,i) = std([molecules_rne_plus1(1,i),molecules_rne_plus2(1,i)])/1.414;
    molecules_rne_minus(2,i) = std([molecules_rne_minus1(1,i),molecules_rne_minus2(1,i),molecules_rne_minus3(1,i),molecules_rne_minus4(1,i)])/2; 
    molecules_wt_plus(2,i) = std([molecules_wt_plus1(1,i),molecules_wt_plus2(1,i)])/1.414;
    molecules_wt_minus(2,i) = std([molecules_wt_minus1(1,i),molecules_wt_minus2(1,i),molecules_wt_minus3(1,i),molecules_wt_minus4(1,i)])/2;
    
    molecules_rne_plus(3,i) = mean(((([molecules_rne_plus1(2,i),molecules_rne_plus2(2,i)])-20022)/1125.7)+1);
    molecules_rne_minus(3,i) = mean(((([molecules_rne_minus1(2,i),molecules_rne_minus2(2,i),molecules_rne_minus3(2,i),molecules_rne_minus4(2,i)])-20022)/1125.7)+1); 
    molecules_wt_plus(3,i) = mean(((([molecules_wt_plus1(2,i),molecules_wt_plus2(2,i)])-20022)/1125.7)+1);
    molecules_wt_minus(3,i) = mean(((([molecules_wt_minus1(2,i),molecules_wt_minus2(2,i),molecules_wt_minus3(2,i),molecules_wt_minus4(2,i)])-20022)/1125.7)+1);
    
    molecules_rne_plus(4,i) = std(((([molecules_rne_plus1(2,i),molecules_rne_plus2(2,i)]-20022)/1125.7)+1))/1.414;
    molecules_rne_minus(4,i) = std(((([molecules_rne_minus1(2,i),molecules_rne_minus2(2,i),molecules_rne_minus3(2,i),molecules_rne_minus4(2,i)]-20022)/1125.7)+1))/2; 
    molecules_wt_plus(4,i) = std(((([molecules_wt_plus1(2,i),molecules_wt_plus2(2,i)]-20022)/1125.7)+1))/1.414;
    molecules_wt_minus(4,i) = std(((([molecules_wt_minus1(2,i),molecules_wt_minus2(2,i),molecules_wt_minus3(2,i),molecules_wt_minus4(2,i)]-20022)/1125.7)+1))/2;
    
    molecules_rne_plus(5,i) = (9/4)*mean((([molecules_rne_plus1(3,i),molecules_rne_plus2(3,i)])-48376)/8096.5);
    molecules_rne_minus(5,i) = (9/4)*mean((([molecules_rne_minus1(3,i),molecules_rne_minus2(3,i),molecules_rne_minus3(3,i),molecules_rne_minus4(3,i)])-48376)/8096.5); 
    molecules_wt_plus(5,i) = (9/4)*mean((([molecules_wt_plus1(3,i),molecules_wt_plus2(3,i)])-48376)/8096.5);
    molecules_wt_minus(5,i) = (9/4)*mean((([molecules_wt_minus1(3,i),molecules_wt_minus2(3,i),molecules_wt_minus3(3,i),molecules_wt_minus4(3,i)])-48376)/8096.5);
    
    molecules_rne_plus(6,i) = std((([molecules_rne_plus1(3,i),molecules_rne_plus2(3,i)])-48376)/8096.5)/1.414;
    molecules_rne_minus(6,i) = std((([molecules_rne_minus1(3,i),molecules_rne_minus2(3,i),molecules_rne_minus3(3,i),molecules_rne_minus4(3,i)])-48376)/8096.5)/2; 
    molecules_wt_plus(6,i) = std((([molecules_wt_plus1(3,i),molecules_wt_plus2(3,i)])-48376)/8096.5)/1.414;
    molecules_wt_minus(6,i) = std((([molecules_wt_minus1(3,i),molecules_wt_minus2(3,i),molecules_wt_minus3(3,i),molecules_wt_minus4(3,i)])-48376)/8096.5)/2;
    
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
title('rne701{\it sodB130-sfGFP}, fixed {\alpha_{m}}','FontSize',14,'FontWeight','bold','Color','g','FontName','Arial')
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
file1 = strcat([folderTitle,'rnesodB130mRNA_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'rnesodB130mRNA_fit.fig']);
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
title('WT{\it sodB130-sfGFP}, fixed {\alpha_{m}}','FontSize',18,'FontWeight','bold','Color','g','FontName','Arial')
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
file1 = strcat([folderTitle,'WTsodB130mRNA_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'WTsodB130mRNA_fit.fig']);
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
file1 = strcat([folderTitle,'rnesodB130Protein_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'rnesodB130Protein_fit.fig']);
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
file1 = strcat([folderTitle,'WTsodB130Protein_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'WTsodB130Protein_fit.fig']);
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
file1 = strcat([folderTitle,'rnesodB130sRNA_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'rnesodB130sRNA_fit.fig']);
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
file1 = strcat([folderTitle,'WTsodB130sRNA_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'WTsodB130sRNA_fit.fig']);
savefig(gcf,file1_fig)

