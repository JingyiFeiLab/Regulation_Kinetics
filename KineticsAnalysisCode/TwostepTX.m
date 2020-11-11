%% Minus sRNA a_m and kx fitting, 
close all

folderTitle = '/Users/reyer/Data/SingleCellEpi/MR156/TwoStep2/';
parentDir = '/Users/reyer/Data/SingleCellEpi/';
strains = [156];
time = {};
DatesV = {'June_6_2020_2','June_6_2020_2'};
graph_title = strcat([DatesV{1},' WT ptsG']);
timestamps = {[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8]};
n_plasmids = 20;
a_m = .09604;
x0 = [a_m, 6.5584];
lb = [.09604, .8084];
ub = [.09604, 16.2584];

GFP = 1;
PTSG = 1; % 1 if you want to make ptsG comparisons
SGRS = 0; %1 if you want to make SgrS comparisons

sigma_check = 0; % If you want to include standard Deviation

mean_gfp = {};
mean_blue = {};
mean_ptsg = {};
mean_green = {};
mean_sgrs = {};
mean_red = {};
mean_violet = {};

std_gfp = {};
std_blue = {};
std_ptsg = {};
std_green = {};
std_sgrs = {};
std_red = {};
std_violet = {};

se_blue = {};
se_green = {};
se_red = {};


for i = 1:length(DatesV)
    strain = strcat(['MR',num2str(strains)]);
        

    strainDir = strcat([parentDir,strain]);
    if exist(strainDir,'dir')~=7
        mkdir(parentDir,strain{i})
    end
    
   
    
    dataDir = strcat([parentDir,strain,'/',DatesV{i},'/']);
    compareDir = strcat([dataDir,'/Fits']);
    if exist(compareDir,'dir')~=7
        mkdir(dataDir,strcat(['Fits']))
    end
    dataFile = strcat([dataDir,'VThresh_BackgroundSummaryData.mat']);
    load(dataFile)
    time{i} = int_table{timestamps{i},1};
    
    if GFP == 1
        mean_gfp{i} = copy_table{timestamps{i},2};
        mean_blue{i} = int_table{timestamps{i},2};
        std_gfp{i} = copy_table{timestamps{i},3};
        std_blue{i} = int_table{timestamps{i},3};
        se_blue{i} = int_table{timestamps{i},4};
    end
    
    if SGRS == 1
        mean_sgrs{i} = copy_table{timestamps{i},6};
        mean_red{i} = int_table{timestamps{i},8};
        std_sgrs{i} = copy_table{timestamps{i},7};
        std_red{i} = int_table{timestamps{i},9};
        se_red{i} = int_table{timestamps{i},10};
    end
    
    if PTSG == 1
        mean_ptsg{i} = copy_table{timestamps{i},4};
        mean_green{i} = int_table{timestamps{i},5};
        std_ptsg{i} = copy_table{timestamps{i},5};
        std_green{i} = int_table{timestamps{i},6};
        se_green{i} = int_table{timestamps{i},7};
    end
    time
end

filename = strcat([compareDir,'/fit.mat']);

compare_green = zeros(length(DatesV),length(timestamps{1,1}));
compare_blue = zeros(length(DatesV),length(timestamps{1,1}));
compare_red = zeros(length(DatesV),length(timestamps{1,1}));

error_green = zeros(1,length(time{1,1}));
error_blue = zeros(1,length(time{1,1}));
error_red = zeros(1,length(time{1,1}));

replicate_green = zeros(1,length(time{1,1}));
replicate_blue = zeros(1,length(time{1,1}));
replicate_red = zeros(1,length(time{1,1}));

for i = 1:length(DatesV)
    %for j = 1:length(time{1,1})
    for j = 1:length(timestamps{1})
        if sum(j == [9] > 0) && i == 4
            compare_green(i,j) = mean_green{1,2}(j);
            compare_blue(i,j) = mean_blue{1,2}(j);
        else
            compare_green(i,j) = mean_green{1,i}(j);
            compare_blue(i,j) = mean_blue{1,i}(j);
        end
        
%         compare_green(i,j) = mean_green{1,i}(j);
%         compare_blue(i,j) = mean_blue{1,i}(j);
        
        if SGRS == 1
            compare_red(i,j) = mean_red{1,i}(j);
        else
            compare_red(i,j) = 0;
        end
    
    end
end

compare_green = ((compare_green+1.3257e+05)/6.0368e+03)+1;


for i = 1:length(timestamps{1})
    replicate_green(i) = mean(compare_green(:,i));
    replicate_blue(i) = mean(compare_blue(:,i));
    replicate_red(i) = mean(compare_red(:,i));
    
    error_green(i) = std(compare_green(:,i))/1.414;
    error_blue(i) = std(compare_blue(:,i))/1.414;
    error_red(i) = std(compare_red(:,i))/2;
    
end



time = time{1}'*60; % update time axis, unit in second 
t_m = 360; % GFP maturation time, assuming 10min, unit in second
time_shift = time - t_m;

protein_indices = [];
for i = 1:length(time_shift)
    if time_shift(i) >= 0
        protein_indices = [protein_indices i];
    end
end

% name = 'xxxx'; put the name of the data set here
%manX in Rne701
if SGRS == 1
    sRNA_data = replicate_red'; 
else
    sRNA_data = zeros(1,length(time))';
end
mRNA_data = replicate_green'; % mRNA data corresponding to each t in WT
protein_data = replicate_blue'; % protein data corresponding to each t in WT


% known parameters



if strains == 187 || strains == 184 || strains == 238 %manX
    b_m = .0032479;
    kx = 0.0267; %[P]/[m] at steady state
    k_elon = 0.0641;
    %k_elon = 0.0844;
elseif strains == 156 || strains == 162 || strains == 60 || strains == 226 || strains == 227 %ptsG
    b_m = 3.25E-03;
    kx = 0.002; %[P]/[m] at steady state
    %k_elon = 0.0643;
    k_elon = 0.0844;
    t_elon = 0.0574;
elseif strains == 188 || strains == 192 %purR
    b_m = .01579;
    kx = 0.006; %[P]/[m] at steady state
    k_elon = 0.0642;
elseif strains == 196 || strains == 201 || strains == 213 || strains == 208 %asd 
    b_m = 0.1636/60;
    kx = 0.0139; %[P]/[m] at steady state
elseif strains == 94 %asd
    b_m = 0.0172;
    kx = 0.0139; %[P]/[m] at steady state
elseif strains == 97 %asd
    b_m = 0.0304;
    kx = 1E-4; %[P]/[m] at steady state
elseif strains == 105 %asd
    b_m = .0267;
    kx = 1E-6; %[P]/[m] at steady state
elseif strains == 148 %asd
    b_m = 0.0139;
    kx = 0.0139; %[P]/[m] at steady state
elseif  strains == 243 || strains == 244 %asd
    b_m = 0.005404666667;
    kx = 0.0139; %[P]/[m] at steady state
    k_elon = 0.0592;
elseif  strains == 241 || strains == 242 %asd
    b_m = 0.0050;
    kx = 0.0139; %[P]/[m] at steady state
    k_elon = 0.0592;
end


if strains == 0 %strains == 184 || strains == 162 || strains == 192 || strains == 201 || strains == 208   % rne701 mutants
    b_s = 0.0017/2; % for RNAse E, b_s is reduced by 2 fold
    molecules0M = [0 mRNA_data(3) 0 mean([protein_data(protein_indices(1)),protein_data(protein_indices(2))])];
    molecules0P = [0 mRNA_data(1) 0 protein_data(protein_indices(1))];
    tv = linspace(time(3), max(time));
    tv2 = linspace(time_shift(protein_indices(1)), time_shift(protein_indices(length(protein_indices))));
    mRNA_indices = 3:8;
else
    b_s = 0.001;
    %molecules0 = [0 mRNA_data(1) 0 protein_data(protein_indices(1))];
    molecules0M = [0 mRNA_data(1) 0 protein_data(protein_indices(1))];
    molecules0P = [0 mRNA_data(1) 0 protein_data(protein_indices(1))];
    tv = linspace(time(1), max(time));
    tv2 = linspace(time_shift(protein_indices(1)), time_shift(protein_indices(length(protein_indices))));
    mRNA_indices = 1:8;
end


b_p = 0.00018; % protein degradation time, assuming it's the doubling time of the cell

b_e = 0;

% fitting parameters
 % set up the parameter space, degradation rate on nacked mRNA
 
%  ini = 1e-5;   % set up the parameter space, ribosome per mRNA
% rs = 0.44;
% ka = 0.0278;
% b_e = 0.0398;
if SGRS == 1
    ka = logspace (-10,-4,10);
    kd = logspace(-7,-1,10);
    rc = linspace (0,1,5);
    a_s = 1.6877E3;
elseif SGRS == 0
    ka = 0;
    kd = 0;
    rc = 0;
    a_s = 0;
end
 % mRNA transcription, A.U./sec

kx_s = 10000;
b_ms = 0.00457*exp(-kx_s/0.0047)+b_m;
ydataVar = [error_red',error_green',error_blue'];
ydata = [sRNA_data, mRNA_data, protein_data];


%x0 = kx;
options = optimoptions('fmincon','Algorithm','sqp','MaxFunEvals',10000000);
funcToFitM = @(x0,time) TwostepTX_fit(x0,time,k_elon,t_elon,b_m,b_p,n_plasmids,molecules0M);
funcToFitP = @(x0,time) TwostepTX_fit(x0,time,k_elon,t_elon,b_m,b_p,n_plasmids,molecules0P);
funcToFit = @(x0,time) TwostepTX_fit(x0,time,k_elon,t_elon,b_m,b_p,n_plasmids,molecules0M);
subindex = @(A, r, c) A(r, c); 

%funcWeightedLeastSquares = @(x) sum([sum(((subindex(funcToFitM(x, time(mRNA_indices)), ':', 2) - ydata(mRNA_indices,2)).^2)),sum(((subindex(funcToFitP(x, time_shift(protein_indices)), 1:length(protein_indices), 4) - ydata(protein_indices,3)).^2))]);
funcWeightedLeastSquares = @(x) sum([sum(((1/max(compare_green(:)))*(subindex(funcToFitM(x, time(mRNA_indices)), ':', 2) - ydata(mRNA_indices,2)).^2)),sum(((1/max(compare_blue(:)))*(subindex(funcToFitP(x, time_shift(protein_indices)), 1:length(protein_indices), 3) - ydata(protein_indices,3)).^2))]);
%funcWeightedLeastSquares = @(x) sum([sum(((1/1)*(subindex(funcToFitM(x, time(mRNA_indices)), ':', 2) - ydata(mRNA_indices,2)).^2)),sum(((1/1)*(subindex(funcToFitP(x, time_shift(protein_indices)), 1:length(protein_indices), 3) - ydata(protein_indices,3)).^2))]);




A = [];
b = [];
Aeq = [];
beq = [];
5
xMin = fmincon(funcWeightedLeastSquares,x0,A,b,Aeq,beq,lb,ub);
6

%xMin = [6611.3792.1842,0.0020078];
pbest = xMin;
%xMin = [0.1776 30.0000];
Cfit = TwostepTX_fit(xMin,tv,k_elon,t_elon,b_m,b_p,n_plasmids,molecules0M);
Cfit2 = TwostepTX_fit(xMin,tv2,k_elon,t_elon,b_m,b_p,n_plasmids,molecules0P);


% figure(1)
% plot(tv,Cfit(:,2),'--g','LineWidth',3)
% hold on
% for ig = 1:length(DatesV)
%     scatter(time,compare_green(ig,:),'og','LineWidth',3)
%     hold on
% end
% title(strcat([graph_title,'mRNA Fit',DatesV]),'FontSize',32)
% xlabel('Time after Induction (min)','FontSize',24)
% % set(gca,'YLim',[0 2E6])
% % set(gca,'XLim',[0 12])
% set(gca,'FontSize',18)
% set(gcf,'position',[835,883,868,667])
% lgd = legend(strcat(['a_m=',num2str(pbest(1))]),strcat(['kx=',num2str(pbest(2))]));
% lgd.Location = 'northeast';
% lgd.FontSize = 14;
% file1 = strcat([compareDir,'/',strain,'mRNA_fit']);
% set(gcf,'PaperPositionMode','auto')
% print(file1,'-painters','-depsc','-r0')
% set(gcf,'PaperPositionMode','auto')
% print(file1,'-dpng','-r0')
% file1_fig = strcat([compareDir,'/',strain,'mRNA_fit.fig']);
% savefig(gcf,file1_fig)

figure(2)
e3 = errorbar(time,replicate_green,error_green,'s');
e3.MarkerFaceColor = [0,.3,0];
e3.MarkerSize = 15;
e3.Color = [0,.3,0];
e3.LineWidth = 2;
hold on 

plot(tv,Cfit(:,2),'--g','LineWidth',3)

title(graph_title,'FontSize',32)
xlabel('Time after Induction (s)','FontSize',24)
set(gca,'YLim',[0 2E6/2025.7])
set(gca,'XLim',[0,1800])
set(gca,'FontSize',18)
set(gcf,'position',[835,883,868,667])
lgd = legend(strcat(['a_m=',num2str(pbest(1))]),strcat(['kx=',num2str(pbest(2))]));
lgd.Location = 'northeast';
lgd.FontSize = 14;
file1 = strcat([folderTitle,'WTptsGmRNA_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'WTptsGmRNA_fit.fig']);
savefig(gcf,file1_fig)

% figure(2)
% plot(tv2+t_m,Cfit2(:,4),'--b','LineWidth',3)
% hold on
% for ig = 1:length(DatesV)
%     scatter(time,compare_blue(ig,:),'ob','LineWidth',3)
%     hold on
% end
% title(strcat([graph_title,'Protein Fit',DatesV]),'FontSize',32)
% xlabel('Time after Induction (min)','FontSize',24)
% % set(gca,'YLim',[0 2E6])
% % set(gca,'XLim',[0 12])
% set(gca,'FontSize',18)
% set(gcf,'position',[835,883,868,667])
% lgd = legend(strcat(['a_m=',num2str(pbest(1))]),strcat(['kx=',num2str(pbest(2))]));
% lgd.Location = 'northeast';
% lgd.FontSize = 14;
% file1 = strcat([compareDir,'/',strain,'protein_fit']);
% set(gcf,'PaperPositionMode','auto')
% print(file1,'-painters','-depsc','-r0')
% set(gcf,'PaperPositionMode','auto')
% print(file1,'-dpng','-r0')
% file1_fig = strcat([compareDir,'/',strain,'protein_fit.fig']);
% savefig(gcf,file1_fig)

figure(4)
plot(tv2+t_m,Cfit2(:,4),'--b','LineWidth',3)
hold on
e7 = errorbar(time,replicate_blue,error_blue,'s');
e7.MarkerFaceColor = [.05,.05,.53];
e7.MarkerSize = 15;
e7.Color = [.05,.05,.53];
e7.LineWidth = 2;



title(graph_title,'FontSize',32)
xlabel('Time after Induction (s)','FontSize',24)
set(gca,'YLim',[0 1.8E7])
set(gca,'XLim',[0 1800])
set(gca,'FontSize',18)
set(gcf,'position',[835,883,868,667])
lgd = legend(strcat(['a_m=',num2str(pbest(1))]),strcat(['kx=',num2str(pbest(2))]));
lgd.Location = 'northeast';
lgd.FontSize = 14;
file1 = strcat([folderTitle,'WTptsGProtein_fit']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'WTptsGProtein_fit.fig']);
savefig(gcf,file1_fig)

save(filename)


