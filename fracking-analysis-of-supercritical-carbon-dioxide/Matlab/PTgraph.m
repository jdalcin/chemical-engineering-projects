R=8.314;
i=1;
b=4.28*10^-5;
a=.366;
E=[218.2:2:304.2];
PT=zeros(1,length(E));
for T = 218.2:2:304.2
    V = 5*10^-5:0.005*10^-4:5*10^-5;
    P = @(V,T) (((R*T)./(V-b)) - (a./(V.^2)));
    F_spinodal = [R*T -2*a 4*a*b -2*a*b^2];
    V_ans(:,i) = roots(F_spinodal);
    Vg(i) = V_ans(1,i);
    Vl(i) = V_ans(2,i);
    pt=(R*T)/(V_ans(1,i)-b)-a/V_ans(1,i)^2;
    PT(i)=pt;
    i = i+1;
end
plot(E,PT,'-k')
hold on
T = [298 390 437 461.8 475 498 518 533 544 556 566 571 585 592 599 608 613.5 619.3 629 637.8];
Ttrue = [298 409.6 467.2 497.6 513.79 542 566.4 587.4 598.2 612.8 625 631.2 648.3 656.8 665.4 676.4 683.16 690.25 702.1 709.5];
P= [138000 0.5*10^6 0.9*10^6 1180000 1.4*10^6 1.8*10^6 2222000 2.6*10^6 2.9*10^6 3264000 3.6*10^6 3.8*10^6 4306000 4.6*10^6 4.9*10^6 5348000 5.6*10^6 5.9*10^6 6390000 6.89*10^6];
Tfinal = [709.5 298.15];
Pfinal = [6.89*10^6 6.89*10^6];
plot(T,P,'-k',Ttrue,P,'-b')
hold on
plot(Tfinal,Pfinal,'-r');
xlabel('Temperature(K)');
ylabel('Pressure(Pa)');
legend('reversible path-adiabatic process','irreversible path-adiabatic process','isobaric cooling');
xlim([218 710])
