close all
clear all
file = '/Users/reyer/Data/SingleCellEpi/growth/';
time1 = [0,20,30,65,100,120,140,160,180,205];
time2 = [0,15,30,45,60,75,90,105,120];

delta1 = [.085,.114,.124,.153,.186,.202,.235,.227,.233,.235];
wt1 = [.105,.118,.133,.164,.214,.235,.248,.254,.256,.252];
delta2 = [.194,.216,.237,.259,.296,.328,.353,.350,.380];
wt2 = [.172,.190,.208,.231,.257,.288,.318,.340,.353];

y1 = log2(delta1(1:6)/delta1(1));
y6 = log2(wt1(1:6)/wt1(1));
t5_line = 0:125;
p1 = polyfitn([time1(1:6)';time2(1:6)'],[y1';y6'],1);
y1_fit = p1.Coefficients(1)*t5_line;
y5_spread(1,:) = y1_fit;

% y2 = log2(delta1(7:10)/delta1(1));
% t5_line = 0:255;
% p2 = polyfitn(time1(7:10)',y2,1);
% y2_fit = p2.Coefficients(1)*t5_line+p2.Coefficients(2);
% y5_spread(1,:) = y2_fit;
% 
% y3 = log2(wt1(1:6)/wt1(1));
% t5_line = 0:255;
% p3 = polyfitn(time1(1:6)',y3,1);
% y3_fit = p3.Coefficients(1)*t5_line+p3.Coefficients(2);
% y5_spread(1,:) = y3_fit;
% 
% y4 = log2(wt1(6:10)/wt1(1));
% t5_line = 0:255;
% p4 = polyfitn(time1(6:10)',y4,1);
% y4_fit = p4.Coefficients(1)*t5_line+p4.Coefficients(2);
% y5_spread(1,:) = y4_fit;

% y5 = log2(wt2(1:7)/wt2(1));
% t5_line = 0:160;
% p5 = polyfitn(time2(1:7)',y5,1);
% y5_fit = p5.Coefficients(1)*t5_line+p5.Coefficients(2);

growth_error = [];
growth_mean = [];

for i = 1:6
    growth_mean = [growth_mean mean([log2(wt1(i)/(wt1(1))),log2(delta1(i)/(delta1(1)))])];
    growth_error = [growth_error std([log2(wt1(i)/(wt1(1))),log2(delta1(i)/(delta1(1)))])];
end


% y1_spread(i,:) = normrnd(p1.Coefficients(1),p1.ParameterStd(1))*t1_line+p1.Coefficients(2);
%     y2_spread(i,:) = normrnd(p2.Coefficients(1),p2.ParameterStd(1))*t2_line+p2.Coefficients(2);





figure(1)
scatter(time1(1:6),growth_mean,100,'filled','g')
hold on
errorbar(time1(1:6),growth_mean,growth_error,'-g','LineWidth',2,'LineStyle','none')
hold on
plot(t5_line,y1_fit,'k','LineWidth',1)
hold on
plot(1:89,ones(length(1:89)),'--k')
hold on 
plot(89*ones(length(0:1/89:1)),0:1/89:1,'--k')
title('Growth Curve','FontSize',18,'FontWeight','bold','FontName', 'Arial')
ylabel('{Log_{2}(OD_{600}/OD_{600 t = 0})}','FontSize',18,'FontName', 'Arial')
xlabel('Time (min)','FontSize',18,'FontName', 'Arial')
xt = get(gca, 'XTick');
lgd = legend('Doubling time = 89.05 min');
lgd.FontSize = 10;
lgd.FontName = 'Arial';
lgd.Location = 'northwest';
set(gca,'XLim',[0 125])
set(gca,'FontSize',10,'FontWeight','bold')
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[318,235,556,247])
file1 = strcat([file,'growthCurve']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
print(file1,'-painters','-depsc','-r0')
file1_fig = strcat([file,'growthCurve.fig']);
savefig(gcf,file1_fig)
