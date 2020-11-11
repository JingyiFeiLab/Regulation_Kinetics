folderTitle = '/Users/reyer/Data/SingleCellEpi/qPCR_Calibrations';
close all

fluor_ptsG = [18820.895
%150596.5095
1323813.47
748940.706
1887883.785
%2291633.69
1767069.87];

fluor_ptsG_e = [1279.813776
%45364.17745
60538.21168
299190.6055
744660.6779
%175291.1205
406989.648];

fluor_SgrS = [41228.51175
417710.0507
931980
1528089.239
1680250.807];

fluor_SgrS_e = [7252.810872
189914.6677
308320.0
53866.31837
170060.3221];

copy_ptsG = [12.6718203
%387.2489972
510.3442885
395.2501471
546.4683326
%11.651667
411.6890778];

copy_ptsG_e = [1.614771959
%22.14850118
69.48233593
81.04005066
104.5351164
%0.5588510894
37.78796732];

copy_SgrS = [0.1643092015
73.20369834
117.712673
184.0887551
225.2314869];

copy_SgrS_e = [0.09433489413
5.952734536
32.62096117
14.38005541
41.10310832];

x_ptsG = 0:600;
x_SgrS = 0:230;

fun = @(x,xdata) x(1)*xdata;
%fun = @(x,xdata) x(1)*(-1./(xdata)) - x(2);
x0 = 10000;


p1 = polyfitn(copy_ptsG,fluor_ptsG,1);
%p1.Coefficients(2) = 1;
y1_spread = zeros(1000,length(x_ptsG));
y1_fit = p1.Coefficients(1)*x_ptsG+p1.Coefficients(2);


p2 = polyfitn(copy_SgrS,fluor_SgrS,1);
y2_spread = zeros(1000,length(x_SgrS));
p2.Coefficients(2) = 0;
p1.Coefficients(1) = 7.5067e+03;
y2_fit = p2.Coefficients(1)*x_SgrS+p2.Coefficients(2);

for i = 2:1000
    
    y1_spread(i,:) = normrnd(p1.Coefficients(1),p1.ParameterStd(1)/2)*x_ptsG+p1.Coefficients(2);
    y2_spread(i,:) = normrnd(p2.Coefficients(1),p2.ParameterStd(1)/2)*x_SgrS+p2.Coefficients(2);
    
    
end

y1_error = zeros(1,length(x_ptsG));
y2_error = zeros(1,length(x_SgrS));

for i = 1:length(x_ptsG)
    y1_error(i) = std(y1_spread(:,i))*2;
end

for i = 1:length(x_SgrS)
    y2_error(i) = std(y2_spread(:,i))*2;
end


RegressionLine = y1_fit(int32(copy_ptsG')+1);
y_plus = fluor_ptsG;
RMSE = sqrt(mean((y_plus'-RegressionLine).^2));
% R2 between regression line and y
SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
SS_Y = sum((y_plus'-mean(y_plus')).^2);
SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y_plus'-mean(y_plus')));
R_squared_plus = SS_XY/sqrt(SS_X*SS_Y);

figure(1)
scatter(copy_ptsG,fluor_ptsG,20,'filled')
hold on
errorbarxy(copy_ptsG,fluor_ptsG,copy_ptsG_e,fluor_ptsG_e)
hold on
plot(x_ptsG,y1_fit)
% hold on
% shadedErrorBar(x_ptsG,y1_fit,y1_error,'lineProps', {'Color','r'})
% hold on
legend(strcat(['R^{2} = ',num2str(R_squared_plus)]),'Location','Northwest')
%set(gca,'XScale','log')
xlabel('ptsG Copy Number','FontSize',12)
ylabel('ptsG Fluorescence','FontSize',12)
title('ptsG CT vs Fluorescence','FontSize',18)
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
set(gca,'XLim',[-10 600])
set(gca,'YLim',[-100000 2.5E6])
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[626,281,248,201])
file1 = strcat([folderTitle,'copy_v_fluor_ptsG']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'copy_v_fluor_ptsG.fig']);
savefig(gcf,file1_fig)

RegressionLine = y2_fit(int32(copy_SgrS')+1);
y_plus = fluor_SgrS;
RMSE = sqrt(mean((y_plus'-RegressionLine).^2));
% R2 between regression line and y
SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
SS_Y = sum((y_plus'-mean(y_plus')).^2);
SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y_plus'-mean(y_plus')));
R_squared_plus = SS_XY/sqrt(SS_X*SS_Y);
R_squared_plus = SS_XY/sqrt(SS_X*SS_Y);


figure(2)
scatter(copy_SgrS,fluor_SgrS,20,'filled')
hold on
errorbarxy(copy_SgrS,fluor_SgrS,copy_SgrS_e,fluor_SgrS_e)
hold on
plot(x_SgrS,y2_fit)
hold on
% shadedErrorBar(x_SgrS,y2_fit,y2_error,'lineProps', {'Color','r'})
% hold on
legend(strcat(['R^{2} = ',num2str(R_squared_plus)]),'Location','Northwest')
%set(gca,'XScale','log')
xlabel('SgrS Copy Number','FontSize',12)
ylabel('SgrS Fluorescence','FontSize',12)
title('SgrS CT vs Fluorescence','FontSize',18)
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
set(gca,'XLim',[-10 270])
set(gca,'YLim',[-100000 2E6])
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[626,281,248,201])
file1 = strcat([folderTitle,'copy_v_fluor']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
print(file1,'-painters','-dpdf','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([folderTitle,'copy_v_fluor.fig']);
savefig(gcf,file1_fig)