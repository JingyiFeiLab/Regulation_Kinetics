function [history,searchdir,xMin] = JustsRNAfit
ydata2 = [29185.23255
743733.1437
1124959.757
1901744.254
2233843.056];

ydata2 = (ydata2 + 5.2879e+04)/7.3302e+03;
%time1 = [0 6 12 24 30];
time2 = [0 6 12 18 24];
time2 = time2*60;
a_s = 10;
b_s = 0.001022905;
%b_s = 1E-4;
%b_s = .0017;
%x0 = [a_s,b_s];
x0 = a_s;
molecules0 = ydata2(1);

funcToFit = @(x0, time) JustsRNA_fit(x0,time2,b_s,molecules0);
subindex = @(A, r, c) A(r, c); 
%funcWeightedLeastSquares = @(x) sum([sum(((subindex(funcToFit(x, time1), ':', 1) - ydata1(:,1)).^2)),sum(((subindex(funcToFit(x, time2), ':', 1) - ydata2(:,1)).^2))]);
funcWeightedLeastSquares = @(x) sum(((subindex(funcToFit(x, time2), ':', 1) - ydata2(:,1)).^2));

lb = [1];
ub = [1000];
A = [];
b = [];
Aeq = [];
beq = [];

history.x = [];
history.fval = [];
searchdir = [];

options = optimoptions(@fminunc,'OutputFcn',@outfun,'Display','iter','Algorithm','quasi-newton');
%options = optimoptions('patternsearch','Display','iter','PlotFcn',@psplotbestf);
[xMin,fval,exitflag,output,grad,hessian] = fminunc(funcWeightedLeastSquares,x0,options);
%xMin = patternsearch(funcWeightedLeastSquares,x0,A,b,Aeq,beq,lb,ub,[],options);
%xMin = patternsearch(funcWeightedLeastSquares,x0,A,b,Aeq,beq,lb,ub,[]);

    function stop = outfun(x,optimValues,state)
        stop = false;
        
        switch state
            case 'init'
                hold on
            case 'iter'
                % Concatenate current point and objective function
                % value with history. x must be a row vector.
                history.fval = [history.fval; optimValues.fval];
                history.x = [history.x; x];
                % Concatenate current search direction with
                % searchdir.
                searchdir = [searchdir;...
                    optimValues.searchdirection'];
                %                     plot(x(1),x(2), x(3), 'o');
                %                     % Label points with iteration number and add title.
                %                     % Add .15 to x(1) to separate label from plotted 'o'
                %                     text(x(1)+.15,x(2),x(3),...
                %                         num2str(optimValues.iteration));
                %                     title('Sequence of Points Computed by fmincon');
            case 'done'
                hold off
            otherwise
        end
    end

pbest = xMin;
tv = linspace(min(time2), max(time2));
Cfit = JustsRNA_fit(xMin,tv,b_s,molecules0);
sqrt(diag(inv(hessian)))

figure(1)
plot(tv,Cfit(:,1),'--r','LineWidth',3)
hold on
scatter(time2,ydata2,'or','LineWidth',3)
title(strcat(['WT SgrS, , a_s = ',num2str(xMin), ', B_m = ',num2str(b_s)]))
end



