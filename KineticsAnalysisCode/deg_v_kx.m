close all
file = '/Users/reyer/Data/SingleCellEpi/degvkx/next/';

b_ms = [.0027929/.0027929,.0173/.0027929,.0188/.0027929,.0481/.0027929,.0139/.0027929];
kx_s = [0.00250/0.00250,.0003225/0.00250,.000680/0.00250,1.4493e-05/0.00250,.00094025/0.00250];
% 
% b_ms_e = [.0002/.0027929,.0018/.0027929,.0042/.0027929,.0081/.0027929,.0006/.0027929];
% kx_s_e = [0.00011/0.00250,.000134/0.00250,.000318/0.00250,1.533e-05/0.00250,0.0005502/0.00250];

b_ms = [.0027929,.0173,.0188,.0481,.0139];
kx_s = [0.002316,.0003225,.000680,1.4493e-05,.00094025];

b_ms_e = [.0002,.0018,.0042,.0081,.0006];
kx_s_e = [0.00011,.000134,.000318,1.533e-05,0.0005502];

% b_ms = [.0173/.0139,.0188/.0139,.0481/.0139,.0139/.0139];
% kx_s = [.0003225/.00094025,.000680/.00094025,1.4493e-05/.00094025,.00094025/.00094025];
% 
% b_ms_e = [.0018/.0139,.0042/.0139,.0081/.0139,.0006/.0139];
% kx_s_e = [.000134/.00094025,.000318/.00094025,1.533e-05/.00094025,0.0005502/.00094025];

% b_ms = [.00027929,.0173,.0188,.0481,.0139,.003029,.01579];
% kx_s = [0.00250,.0003225,.000680,1.4493e-05,.00094025,.0030433,.0009375];
% 
% b_ms_e = [.0002,.0018,.0042,.0081,.0006,.0002,.0013];
% kx_s_e = [0.00011,.000134,.000318,1.533e-05,0.0006502,.00011,.0003288];


fun = @(x,xdata) x(1)*(-1*log(xdata)) - x(2);
%fun = @(x,xdata) x(1)*(-1./(xdata)) - x(2);
x0 = [.001,0];

[x,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(fun,x0,kx_s,b_ms);
ci = nlparci(x,resid,'jacobian',J);
% x(2) = -1;
% ci(2,:) = ci(2,:) +.4982;




% p1 = polyfitn(t1,y1,1);
% y1_fit = p1.Coefficients(1)*t1_line+p1.Coefficients(2);


figure(1)
errorbarxy(kx_s,b_ms,kx_s_e,b_ms_e,'-b','LineWidth',2)
hold on
scatter(kx_s,b_ms,50,'filled','b')
hold on
x_fit = logspace(-5,-2,100);
y1_spread = zeros(50,length(x_fit));
y_fit = x(1)*(-1*log(x_fit))-x(2);
%y_fit = x(1)*(-1./(x_fit))-x(2);
y1_spread(1,:) = y_fit;

for i = 2:50
    
    %y1_spread(i,:) = (ci(1,1)+0.5*((ci(1,2)-ci(1,1))*rand(1)))*(-1*log(x_fit))-(ci(2,1)+0.5*((ci(2,2)-ci(2,1))*rand(1)));
    y1_spread(i,:) = (ci(1,1)+((ci(1,2)-ci(1,1))*rand(1)))*(-1*log(x_fit))-(ci(2,1)+((ci(2,2)-ci(2,1))*rand(1)));
    %y1_spread(i,:) = (ci(1,1)+((ci(1,2)-ci(1,1))*rand(1)))*(-1./(x_fit))-(ci(2,1)+((ci(2,2)-ci(2,1))*rand(1)));
   
end
    
y1_error = zeros(1,length(x_fit));


for i = 1:length(x_fit)
    y1_error(i) = std(y1_spread(:,i));
end

shadedErrorBar(x_fit,y_fit,.5*y1_error,'lineProps', {'Color','r'})
hold on



plot(x_fit,y_fit,'r','LineWidth',1)
hold on

errorbarxy(kx_s,b_ms,kx_s_e,b_ms_e,'-b','LineWidth',2)
hold on
scatter(kx_s,b_ms,50,'filled','b')
hold on 
%errorbarxy(.567*.0023166,(.00350)/.0027929,.042,(.00036/.0027929),'-g','LineWidth',4)
errorbarxy(0.3885*0.002316,0.004315,.042*.0023166,.00036,'-g','LineWidth',4)
hold on
errorbarxy(kx_s,b_ms,kx_s_e,b_ms_e,'-b','LineWidth',2)

title('Degradation vs. Translation Rate','FontSize',12,'FontWeight','bold')
ylabel('Degradation Rate (s^{-1})','FontSize',12,'FontWeight','bold')
xlabel('Translation Rate (s^{-1})','FontSize',12,'FontWeight','bold')
% ylabel('Degradation Rate Fold Change','FontSize',12,'FontWeight','bold')
% xlabel('Translation Rate Reduction','FontSize',12,'FontWeight','bold')
xt = get(gca, 'XTick');
set(gca,'xscale','log')
set(gca,'XLim',[-.1,.01])
set(gca,'yLim',[.0001,.08])
set(gca,'FontSize',10,'FontWeight','bold')
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
set(gca,'XScale','log')
set(gcf,'position',[626,281,248,201])
% set(gca,'YScale','log')
file1 = strcat([file,'kx_v_bm']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
print(file1,'-painters','-depsc','-r0')
file1_fig = strcat([file,'kx_v_bm.fig']);
savefig(gcf,file1_fig)








