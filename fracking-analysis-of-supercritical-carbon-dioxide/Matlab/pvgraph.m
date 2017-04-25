function pvgraph_Jeremy
a = 0.3649;
b = 4.276*10^-5;
R = 8.314;
Tc = 304.25;
 for T = [277 290 304.25 315 331]
    V = 5*10^-5:0.005*10^-4:2*10^-1;
    P = @(V,T) (((R*T)./(V-b)) - (a./(V.^2)));
    dpdv = @(V,T) ((2*a./V.^3) - (R*T./(V-b).^3));
    semilogx(V,P(V,T),'-r');
    axis([5*10^-5 2*10^-1 0 15*10^6]);
    hold on
 end
legend('From Bottom to Topmost Line:','T1 = 277 K','T2 = 290 K','Tc= 304.25 K','T4 = 315 K','T5 = 331 K');
title('Pressure vs. Volume of CO2')
ylabel('Pressure(Pa)')
xlabel('Volume(m^3/mol)')
xlim([6*10^-5 1*10^-1])
Pc = 7391000;
Tc = 304.12;
Fnew = [Pc -(Pc*b+R*Tc) a -a*b];
Vsp_ans = abs(roots(Fnew));
T = 304.12;
%plot(Vsp_ans(1,1),P(Vsp_ans(1,1),T),'o');
%find Psat at different T
T1 = [277 283 290 295 300 304];
a = 0.3649;
b = 4.276*10^-5;
R = 8.314;
P1_guess = 6*10^6;
for i = 1:length(T1)
    P_bin(i) = fsolve(@(P) myfun(P,T1(i)),P1_guess);
    Fs = [P_bin(i) -(P_bin(i)*b+R*T1(i)) a -a*b];
    Vs = (roots(Fs));
    Vs = abs(Vs);
    Vl(i) = min(Vs)
    Vg(i) = max(Vs)
    semilogx(Vl(i),P_bin(i),'-b',Vg(i),P(Vg(i),T1(i)),'-b')
    end
    hold on
    for i = 1:(length(T1))
    x(i) = Vl(i);
    y(i) = P_bin(i);
    end
    j = 1;
    for i = ((length(T1))+1):(2*length(T1))
    x(i) = Vg((length(T1))-j+1);
    y(i) = P_bin((length(T1))-j+1);
    j = j+1;
    end
    semilogx(x,y)
hold on
P_bin;
i = 1;
TS1=[277:.25:304];
for T = 277:.25:304
    V = 5*10^-5:0.005*10^-4:5*10^-5;
    P = @(V,T) (((R*T)./(V-b)) - (a./(V.^2)));
    F_spinodal = [R*T -2*a 4*a*b -2*a*b^2];
    V_ans(:,i) = roots(F_spinodal);
    Vg(i) = V_ans(1,i);
    Vl(i) = V_ans(2,i);
    PS(i)=(R*T)/(V_ans(1,i)-b)-a/V_ans(1,i)^2;
   semilogx(V_ans(1,i),P(V_ans(1,i),T),'-k',V_ans(2,i),P(V_ans(2,i),T),'-k');
   hold on
    i = i+1;
end
l1=integral(@(V)(R*290)/(V-b)-a/V^2,V_ans(3,52),V_ans(2,52),'Arrayvalued',true)
l2=integral(@(V)(R*290)/(V-b)-a/V^2,V_ans(2,52),V_ans(1,52),'Arrayvalued',true)
% plot([7.92*10^-5 2.83*10^-4],[5.03*10^6 5.03*10^6],'-m')
% plot([8.88*10^-5 2.16*10^-4],[6.09*10^6 6.09*10^6],'-m')
semilogx(7.38*10^6, 1.28*10^-4,'-k')
% VS=[Vl 1.28*10^-4 flip(Vg)];
% PS=[PS 7.38*10^6 flip(PS)];
% TS=[TS1 304.25 flip(TS1)];
% plot(VS,PS,'-k')
V1_pb = 0.0178;
P1_pb = 0.138*10^6;

V2_pb = 7.4434*10^-4;
P2_pb = 6.39*10^6;

V3_pb = 1.5265*10^-4;
P3_pb = P2_pb;

V_partb = [V1_pb V2_pb V3_pb];
P_partb = [P1_pb P2_pb P3_pb];
% %adiabatic path graph of CO2
% syms pf vf
% T=[298:10:640];
% a=.366;
% b=4.28*10^-5;
% R=8.314;
% V=zeros(1,length(T)+1);
% V(1)=.0179;
% P=zeros(1,length(T)+1);
% P(1)=1.38*10^5;
% for i=1:length(T)-1
% p=(R*T(i+1))/(vf-b)-a/vf^2-pf; %van der waals
% S=R*(4.457*log(T(i+1)/T(i))+1.045*10^-3*(T(i+1)-T(i))+5.785*10^4*(T(i+1)^-2-T(i)^-2))+R*log((vf-b)/(V(i)-b)); %entropy
% [pfr,vfr]=solve(p==0,S==0);
% V(i+1)=double(vfr);
% P(i+1)=double(pfr);
% end
% V(length(T)+1)=8.4*10^-4;
% P(length(T)+1)=6.89*10^6;
% T=[T 709.5];
% %isobaric path graph of CO2
% pc=6.89*10^6;
% T1=[709.5 708:-10:298];
% P1=6.89*10^6*ones(1,length(T1));
% VF=zeros(1,length(T1));
% for i=1:length(T1)
% vfr=solve(pc==(R*T1(i))/(vf-b)-a/vf^2);
% VF(i)=double(vfr(1));
% end
T = [298 390 437 461.8 475 498 518 533 544 556 566 571 585 592 599 608 613.5 619.3 629 637.8];
Ttrue = [298 409.6 467.2 497.6 513.79 542 566.4 587.4 598.2 612.8 625 631.2 648.3 656.8 665.4 676.4 683.16 690.25 702.1 709.5];
P= [138000 0.5*10^6 0.9*10^6 1180000 1.4*10^6 1.8*10^6 2222000 2.6*10^6 2.9*10^6 3264000 3.6*10^6 3.8*10^6 4306000 4.6*10^6 4.9*10^6 5348000 5.6*10^6 5.9*10^6 6390000 6.89*10^6];
syms k
for i=1:length(T)
s=double(solve(P(i)==(R*T(i))/(k-b)-a/k^2));
K(i)=max(s);
end
for i=1:length(Ttrue)
s=double(solve(P(i)==(R*Ttrue(i))/(k-b)-a/k^2));
Ktrue(i)=max(s);
end
 plot(K,P,'-r',Ktrue,P,'-b',[8.4*10^-4 9.58*10^-5],[6.89*10^6 6.89*10^6],'-k')
xlabel('Volume(m^3/mol)')
ylabel('Pressure(Pa)')
hold off
end

function F = myfun(P,T)
a = 0.3649;
b = 4.276*10^-5;
R = 8.314;
Fs = [P -(P*b+R*T) a -a*b];
Vs = (roots(Fs));
Vl = min(Vs);
Vg = max(Vs);
F = R*T*log((Vg-b)/(Vl-b)) + a*((1/Vg) - (1/Vl)) - P*(Vg - Vl);
end

function chemp(x,y)
a = 0.3649;
b = 4.276*10^-5;
R = 8.314;
for T = 274.12:15:334.12
p = 0.1*10^4:0.001*10^4:2*10^6;
u = ((R*T.*p*b)./(1-p.*b))- 2*a*p + R*T.*log((p*R*T)./(1-p.*b));
semilogx(p,u);
hold on
end
hold off
end