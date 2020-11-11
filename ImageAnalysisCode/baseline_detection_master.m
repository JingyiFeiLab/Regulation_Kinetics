function [baseline_curve,myfit] = baseline_detection_master(cell_spot_negative)
%Get spots per cell first for negative control. Spots per cell is correct
%for fitting


baseline_curve = zeros(length(cell_spot_negative),2);
l = length(cell_spot_negative);

for i=1:l
    baseline_curve(i,1)=length(find(cell_spot_negative(:,3)<=cell_spot_negative(i,3)))/l; % = m
    baseline_curve(i,2)=cell_spot_negative(i,3); %Number of spots per cell
end


myfittype = fittype('a1 + a2*exp(a3*x-a4) + a5*exp(a6*x-a7)',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a1','a2','a3','a4','a5','a6','a7',});
myfit = fit(baseline_curve(:,1),baseline_curve(:,2),myfittype);

end
