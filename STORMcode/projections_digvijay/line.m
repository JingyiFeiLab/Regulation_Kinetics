function out= line(guesses)
m=guesses(1);
c=guesses(2);

load('Xset.mat');
load('Yset.mat');
value=0;
for i=1:numel(Xset);
    value=value+(Yset(i)-m*Xset(i)-c).^2;
end
out=value;
end