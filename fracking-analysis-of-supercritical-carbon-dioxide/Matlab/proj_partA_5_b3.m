
clc;
%close all;
clear all;

R = 8.3144621; 

Tc = 31.1+273.15;
Pc = 7.38e6;

Z = 0.272;
%vc = Z*R*Tc/Pc;

a = 27/64*R*R*Tc*Tc/Pc;
b = 1/8*R*Tc/Pc;

W = 44.01*1e-3;

%GRI-3.0 Mech

%H2
%a_high_k(1,:) = [3.33727920 -4.94024731E-5   4.99456778E-07 -1.79566394E-10  2.00255376E-14 -9.50158922E+02 -3.20502331E+00];
%a_low_k(1,:)  = [2.34433112  7.98052075E-03 -1.94781510E-05  2.01572094E-08 -7.37611761E-12 -9.17935173E+02  6.83010238E-01];

a_high_k(1,:) = [3.85746029 4.41437026E-03 -2.21481404E-06 5.23490188E-10  -4.72084164E-14 -4.87591660E+04 2.27163806];
a_low_k(1,:)  = [2.35677352 8.98459677E-03 -7.12356269E-06 2.45919022E-09  -1.43699548E-13 -4.83719697E+04 9.90105222];



G_func_varV = @(T0,v0)( R*T0*b/(v0-b) -R*T0*log(v0-b)  - 2*a/v0);      

A_func_varV = @(T0,v0)(-R*T0*log(v0-b)  - a/v0);      



%a = 3.640*0.101325;
%b = 0.04267*1e-3;

%td = [-0.10 -0.06 -0.03 0.00 0.03 0.06 0.10];

%td = [0.00];
%td = [-0.1];

td = [-0.10];
%td = [-0.08];

T_arr = 277;

%v_arr = [2e-6:1e-6:11e-5];


%v_arr = [6e-5:1e-6:9e-4];

%rho_arr = [W/9e-4:0.1:W/6e-5];
rho_arr = [1.0/5e-3:1:1.0/6e-5];



figure(1);clf

figure(2);clf

for i=1:numel(T_arr);
    T0 = T_arr(i);
    
    
    
    k =1;
    %Rsp = R/W;
     Hig = R*T0*(a_low_k(k,1) + T0*(a_low_k(k,2)/2.0 + T0*(a_low_k(k,3)/3.0 +  ...
                       T0*( a_low_k(k,4)/4.0 +  T0*a_low_k(k,5)/5.0 ))) + a_low_k(k,6)/T0 );
                   
     Cpig = R*( a_low_k(k,1) + T0*(a_low_k(k,2) + T0*(a_low_k(k,3)  ...
                             + T0*(a_low_k(k,4) +  T0*a_low_k(k,5) ) ) ) );
                         
     Sig = R*( a_low_k(k,1)*log(T0) + T0*(a_low_k(k,2) + T0*(a_low_k(k,3)/2.0  ...
                             + T0*(a_low_k(k,4)/3.0 +  T0*a_low_k(k,5)/4.0 ) ) ) + a_low_k(k,7) );                    
                         
     %cvig = R*( a_low_k(k,1) - 1.0 + Tmin*(a_low_k(k,2) + Tmin*(a_low_k(k,3)  ...
     %                        + Tmin*(a_low_k(k,4) +  Tmin*a_low_k(k,5) ) ) ) );
    mu_arr(1) = 0.0;
    for j=1:numel(rho_arr)
       %v00 = 1.0/rho_arr(j-1);
       %P00 = R*T0/(v00-b) - a/(v00*v00); 
                
       rho0 = rho_arr(j);
       v0 = 1.0/rho0;
       v_arr(j) = v0;
       P0 = R*T0/(v0-b) - a/(v0*v0);
       %Hres = R*T0*(b/(v0-b) - 2.0*a/(v0*R*T0));
       %Sres = R*log(P0*(v0-b)/(R*T0));
       %Hres = Rsp*T0*(rho0*b/(W-rho0*b) - 2.0*a*rho0/(W*R*T0));
       %Sres = Rsp*log(P0*(W-rho0*b)/(rho0*R*T0));
       
       %Gres = R*T0*(b/(v0-b) - 2.0*a/(v0*R*T0) - log(P0*(v0-b)/(R*T0)));
       
       P_arr(j) = P0;
       %mu_arr(j) = mu_arr(j-1) + (P0*v0-P00*v00) + (R*T0*log(v0-b) + a/v0) - (R*T0*log(v00-b) + a/v00);
       
       mu_arr(j) = G_func_varV(T0,v0);
       
       %helm_arr(j) = rho0*A_func_varV(T0,v0);
       helm_arr(j) = rho0*mu_arr(j) - P0;
       
       
    end
    
    %plotyy(rho_arr,P_arr,rho_arr,mu_arr,'plot','plot');
    %plotyy(v_arr,P_arr,v_arr,mu_arr,'semilogy','plot');
    
    %loglog(v_arr,P_arr);
    
    figure(1);
    hold on;
    plot(rho_arr,mu_arr);
    hold off;
    
    
    figure(2);
    hold on;
    plot(rho_arr,helm_arr);
    %xlim([175.0 500.0]);
    %plot(v_arr,helm_arr);
    hold off;
    
    %plot(rho_arr,mu_arr);
    %hold on;
    
    
end




P_arr2 = [min(P_arr):1000:max(P_arr)];

%T_arr2 = [Tc-0.4*Tc:0.5:Tc];

T_arr2 = T_arr;

ind0 = 1;
ind =1;

indpt = 1;

for i=1:numel(T_arr2)
    
    T0 = T_arr2(i);
    
    vg =0.0;
    vl = 0.0;
    binodal =0;
for j=1:numel(P_arr2)
        
    P0 = P_arr2(j);
    %T0 = T_arr(i);
    %T0 = Tc - 0.10*Tc;
    a10 = -(b + R*T0/P0);
    a20 = a/P0;
    a30 = -a*b/P0;
    
    [x10,x20,x30] = cubsolve(a10,a20,a30);
    
    if(isreal(x10) && isreal(x20) && isreal(x30))
       vg = max(max(x10,x20),x30);
       vl = min(min(x10,x20),x30);
       
       lhs = P0*(vg-vl);
       
       rhs = (R*T0*log(vg-b) + a/vg) - (R*T0*log(vl-b) + a/vl);
       diff = lhs - rhs
       if(abs(diff)<0.1)
           %binP(ind) = P0;
           %loglog(vg,P0,'r','MarkerSize', 4*72);
           %loglog(vl,P0,'r','MarkerSize', 4*72);
           
           figure(1);
           hold on;
           
           %h2 = scatter(1.0/vl*1e-6,10.0*fun_G(T0,vl),'b','filled'); set(h2,'SizeData',96);
           
           
           xbarr(ind0) = vg;
           Gbarr(ind0) = G_func_varV(T0,vg);
           h1 = scatter(1.0/vg,Gbarr(ind0),'b','filled'); set(h1,'SizeData',96);
           ind0 = ind0+1;
           
           xbarr(ind0) = vl;
           Gbarr(ind0) = G_func_varV(T0,vl);
           h2 = scatter(1.0/vl,Gbarr(ind0),'b','filled'); set(h2,'SizeData',96);
           ind0 = ind0+1;
           
           hold off;
           binodal = 1;
           break;
       end
       
  
    end
     
   
       
end

    if(binodal ==1)
    a1s = -2.0*a/(R*T0);
    a2s = 4.0*a*b/(R*T0);
    a3s = -2.0*a*b*b/(R*T0);
    
    [x1s,x2s,x3s] = cubsolve(a1s,a2s,a3s);
    
    
    if(isreal(x1s) && x1s>vl && x1s <vg)
        
       P0 = R*T0/(x1s-b) - a/(x1s*x1s); 
       G0 = G_func_varV(T0,x1s);
       figure(1);
       hold on;
       h1 = scatter(1.0/x1s,G0,'b','filled'); set(h1,'SizeData',96);
       hold off;
       
       figure(2);
       hold on;
       h1 = scatter(1.0/x1s,G0/x1s-P0,'b','filled'); set(h1,'SizeData',96);
       hold off;
       
       
       %hs = scatter(1.0/x1s*1e-6,10.0*G0,'b','filled'); set(hs,'SizeData',96);
       xsarr(ind) = x1s;
       Gsarr(ind) = G0;
       ind = ind+1;
       
       
       ptxar(indpt) = 1/x1s;
       ptyar(indpt) = G0/x1s-P0;
       indpt = indpt + 1;
    end
    if(isreal(x2s) && x2s>vl && x2s <vg)
        P0 = R*T0/(x2s-b) - a/(x2s*x2s); 
       G0 = G_func_varV(T0,x2s);
       figure(1);
       hold on;
       h1 = scatter(1.0/x2s,G0,'b','filled'); set(h1,'SizeData',96);
       hold off;
       
       
       figure(2);
       hold on;
       h1 = scatter(1.0/x2s,G0/x2s-P0,'b','filled'); set(h1,'SizeData',96);
       hold off;
       
       xsarr(ind) = x2s;
       Gsarr(ind) = G0;
       ind = ind+1;
       
       ptxar(indpt) = 1/x2s;
       ptyar(indpt) = G0/x2s-P0;
       indpt = indpt + 1;
       
       
    end
    if(isreal(x3s) && x3s>vl && x3s <vg)
        P0 = R*T0/(x3s-b) - a/(x3s*x3s); 
       G0 = G_func_varV(T0,x3s);
       figure(1);
       hold on;
       h1 = scatter(1.0/x3s,G0,'b','filled'); set(h1,'SizeData',96);
       hold off;
       
       figure(2);
       hold on;
       h1 = scatter(1.0/x3s,G0/x3s-P0,'b','filled'); set(h1,'SizeData',96);
       hold off;
       
       
       xsarr(ind) = x3s;
       Gsarr(ind) = G0;
       ind = ind+1;
       
       ptxar(indpt) = 1/x3s;
       ptyar(indpt) = G0/x3s-P0;
       indpt = indpt + 1;
       
    end
    end
    

end


figure(1);
%hc = scatter(1.0/vc*1e-6,10.0*fun_G(Tc,vc),'b','filled'); set(hc,'SizeData',96);


vc = MolarVolumeVDWEOS(Tc,Pc,a,b);

xbarr(ind0) = vc;
Gbarr(ind0) = G_func_varV(Tc,vc);


figure(1);
hold on;
h1 = scatter(1.0/xbarr(ind0),Gbarr(ind0),'b','filled'); set(h1,'SizeData',96);
hold off;

cx = mean(xbarr);
cy = mean(Gbarr);

axx = atan2(Gbarr - cy, xbarr - cx);
[~, order] = sort(axx);
xbarr = xbarr(order);
Gbarr = Gbarr(order);

pol1 = polyfit(xbarr,Gbarr,2)

t1 = [4.1692e-5:1e-6:7.4669e-4];
%t1 = v_arr;
y1 = polyval(pol1,t1);

%plot(xsarr,Psarr);

figure(1);
hold on;
plot(1./t1,y1);
hold off;

%plot(1.0./t1*1e-6,10.0*y1,'-k');
%plot(t1*1e6,y1/1e5,'-k');


%1.0./v_arr*1e-6,10.0*mu_arr


xsarr(ind) = vc;
Gsarr(ind) = G_func_varV(Tc,vc);

cx = mean(xsarr);
cy = mean(Gsarr);

axx = atan2(Gsarr - cy, xsarr - cx);
[~, order] = sort(axx);
xsarr = xsarr(order);
Gsarr = Gsarr(order);

pol2 = polyfit(xsarr,Gsarr,2)

t2 = [8.1692e-5:1e-6:5.4669e-4];
%t2 = v_arr;


y2 = polyval(pol2,t2);

figure(1);
hold on;
plot(1./t2,y2);
hold off;


slope = (ptyar(2) - ptyar(1))/(ptxar(2) - ptxar(1));


for i =1:numel(rho_arr)
    y0 = slope*(rho_arr(i) - ptxar(1)) + ptyar(1);
    diffarr(i) =  helm_arr(i) - y0;
end


figure(3);
plot(rho_arr,diffarr,'-r')
title('Helmholtz vs. Density')
xlabel('Density(mol/m^3)')
ylabel('Helmholtz Energy(J/mol)')


figure(1);
xlim([0 18000]);
ylim([1.4e4 2.4e4]);
hold off;