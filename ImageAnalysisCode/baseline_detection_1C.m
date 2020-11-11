function baseline_curve = baseline_detection_1C(filename);
N=load(filename);
%cor=zeros(length(N), 2);

baseline_curve = zeros(length(N),2);
% for A647

for i=1:length(N)
    m(i,1)=length(find(N(:,1)<=N(i,1)))/length(N);
end

for i=1:length(N)
    baseline_curve(i,1)=m(i,1);
    baseline_curve(i,2)=N(i,1);
end

% for i=1:length(N)
%         if N(i,1)==0
%         cor(i,1)=0;
%         else
%         cor(i,1)=round(N(i,1)-(0.003*exp((m(i,1)-0.566)/0.045)+34.94*exp((m(i,1)-0.566)/0.449)-18.91));
%         end
% end


% for A568
% for i=1:length(N)
%     m(i,2)=length(find(N(:,2)<=N(i,2)))/length(N);
% end
% for i=1:length(N)
%         if N(i,2)==0
%         cor(i,2)=0;
%         else
%         cor(i,2)=round(N(i,2)-(4.81e-5*exp((m(i,2)-0.776)/0.0143)+62.56*exp((m(i,2)-0.776)/0.18)-10.1));
%         end
% end