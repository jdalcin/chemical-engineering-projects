function [Dpmatrix,Vfmatrix]=solver_gravityscale(Dp)
syms D2 vf2
u1=1.4*10^-5;
u2=1.711*10^-5;
rhop1=999.7;
rhof1=.4;
rhof2=1.2;
vf1=250;
% D=linspace(2*10^-6,80*10^-6,100);
D=Dp;
for i=1:length(D)
    [Dp2,Vf2]=solve(u1/(rhof1*vf1*D(i))==u2/(rhof2*vf2*D2),D(i)/vf1^2==D2/vf2^2,D2,vf2);
    Dpmatrix(i,1)=double(Dp2(1));
    Vfmatrix(i,1)=double(Vf2(1));
end
Vfmatrix=Vfmatrix(1);
end