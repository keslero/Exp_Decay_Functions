function [tau_D,m] = find_slope(t,TKE)


ft = fittype('a-1/b*x');
app = fit(t',log10(TKE)',ft);%,'Startpoint',[t(1),TKE(1)]);
tau_D = app.b;
m = app.a;

end