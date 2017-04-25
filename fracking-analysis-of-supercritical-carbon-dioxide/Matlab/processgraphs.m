%adiabatic path graph of CO2
syms pf vf
T=[298:10:640];
a=.366;
b=4.28*10^-5;
R=8.314;
V=zeros(1,length(T)+1);
V(1)=.0179;
P=zeros(1,length(T)+1);
P(1)=1.38*10^5;
for i=1:length(T)-1
p=(R*T(i+1))/(vf-b)-a/vf^2-pf; %van der waals
S=R*(4.457*log(T(i+1)/T(i))+1.045*10^-3*(T(i+1)-T(i))+5.785*10^4*(T(i+1)^-2-T(i)^-2))+R*log((vf-b)/(V(i)-b)); %entropy
[pfr,vfr]=solve(p==0,S==0);
V(i+1)=double(vfr);
P(i+1)=double(pfr);
end
V(length(T)+1)=8.4*10^-4;
P(length(T)+1)=6.89*10^6;
T=[T 709.5];
%isobaric path graph of CO2
pc=6.89*10^6;
T1=[709.5 708:-10:298];
P1=6.89*10^6*ones(1,length(T1));
VF=zeros(1,length(T1));
for i=1:length(T1)
vfr=solve(pc==(R*T1(i))/(vf-b)-a/vf^2);
VF(i)=double(vfr(1));
end
plot(V,P,'-k',VF,P1,'-b')%pvprocessgraph
title('Steps 1 & 2')
xlabel('Volume(m^3/mol)')
ylabel('Pressure(Pa)')
figure
plot(VF,P1,'-b') %isochoric graph
title('Step 2')
xlabel('Volume(m^3/mol)')
ylabel('Pressure(Pa)')
figure
plot(T,V,'-k')%vtprocessgraph
title('step 1')
xlabel('Temperature(K)')
ylabel('Volume(m^3)')
