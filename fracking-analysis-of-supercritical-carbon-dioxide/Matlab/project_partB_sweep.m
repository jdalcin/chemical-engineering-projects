function project_partB_sweep
clc;

clear all;


R = 8.3144621; 

Tc = 31.1+273.15;
Pc = 7.38e6;
W = 44.01*1e-3;

delG_298 = -394359; %J/mol
delH_298 = -393509; %J/mol

plotsyms ={'-r','-g','-b','.-','-k','-m','--'};

a = 27/64*R*R*Tc*Tc/Pc;
b = 1/8*R*Tc/Pc;



Tref = 298;
Pref = 101325;


a_high_k(1,:) = [3.85746029 4.41437026E-03 -2.21481404E-06 5.23490188E-10  -4.72084164E-14 -4.87591660E+04 2.27163806];
a_low_k(1,:)  = [2.35677352 8.98459677E-03 -7.12356269E-06 2.45919022E-09  -1.43699548E-13 -4.83719697E+04 9.90105222];


S_func = @(T0)(R*( (a_low_k(1,1)-1.0)*log(T0) + T0*(a_low_k(1,2) + T0*(a_low_k(1,3)/2.0  ...
                             + T0*(a_low_k(1,4)/3.0 +  T0*a_low_k(1,5)/4.0 ) ) )  ) );

                         
%H_func_varT = @(T0,v0)(R*T0*(a_low_k(1,1) - 1.0 + T0*(a_low_k(1,2)/2.0 + T0*(a_low_k(1,3)/3.0 +  ...
%                  T0*( a_low_k(1,4)/4.0 +  T0*a_low_k(1,5)/5.0 )))) + R*T0*v0/(v0-b) );

%H_func_varV = @(T0,v0)( R*T0*(1.0/(v0-b) + log(v0-b)) - 2.0*a/v0);     




%H_func_varT = @(T0,v0)(R*T0*(a_low_k(1,1) - 1.0 + T0*(a_low_k(1,2)/2.0 + T0*(a_low_k(1,3)/3.0 +  ...
%                  T0*( a_low_k(1,4)/4.0 +  T0*a_low_k(1,5)/5.0 )))) + R*T0*v0/(v0-b) -a/v0 );

%H_func_varV = @(T0,v0)( R*T0*log(v0-b));       





H_func_varT = @(T0,v0)(R*T0*(a_low_k(1,1) - 1.0 + T0*(a_low_k(1,2)/2.0 + T0*(a_low_k(1,3)/3.0 +  ...
                  T0*( a_low_k(1,4)/4.0 +  T0*a_low_k(1,5)/5.0 )))) + R*T0*v0/(v0-b)  );

H_func_varV = @(T0,v0)( R*T0*b/(v0-b) -2*a/v0);      


%delH_func = @(T10,v10,T20,v20,vmid0)( (H_func_varT(T20,vmid0) - H_func_varT(T10,v10) )+ ( H_func_varV(T20,v20) - H_func_varV(T20,vmid0) ) );

%delH_func = @(T10,v10,T20,v20,vmid0)( (H_func_varT(T20,v10) - H_func_varT(T10,v10) )+ ( H_func_varV(T20,v20) - H_func_varT(T20,v10) ) );

delH_func = @(T10,v10,T20,v20,vmid0)(( H_func_varV(T10,v20) - H_func_varV(T10,v10) ) + (H_func_varT(T20,v20) - H_func_varT(T10,v20) ));

T1 = 25 + 273.15;
P1 = 0.138e6; 

v1 = MolarVolumeVDWEOS(T1,P1,a,b);
rho1 = W/v1;

S_1 = S_func(T1);
rhs_1 = R*log(v1-b);

P2_f = 6.89e6;

%P2_arr = [1e6:5e5:7e6];

dp = (P2_f - P1)/20.0;
P2_arr = [P1:dp:P2_f];

T2_1(1) = T1;
T2_2(1) = T1;
v2_1(1) = v1;
for zi = 2:numel(P2_arr)

%P2 = 6.89e6;
P2 = P2_arr(zi);



T2g_arr = [T1:0.1:T1+400];

prev = 1.0;
T2_1(zi) = 1.0;

for i=1:numel(T2g_arr)
   
    T2g = T2g_arr(i);
    
    S_2 = S_func(T2g);
    
    v2 = MolarVolumeVDWEOS(T2g,P2,a,b);
    
    rhs_2 = R*log(v2-b);
    
    lhs = S_2 - S_1;
    %rhs = rhs_2 - rhs_1;
    rhs = rhs_1 - rhs_2;
    
    %plot(T2g,(lhs-rhs));
    %plot(T2g,P2);
    %hold on;
    
    if(i>1 && prev*(lhs-rhs)<0.0)
        
         T2_1(zi) = T2g_arr(i-1) + ( T2g_arr(i)-T2g_arr(i-1) )/( (lhs-rhs) - prev )*( 0.0-prev ) ;
         
         
         %rho2 = W/v2;
         %as2_arr(i) = 1.0/sqrt(rho2*ks(j));
         break;
     end
     
     prev = lhs - rhs;
    
end

% 
% figure(2);clf
% plot(Tc,Pc,'-ro');
% hold on;
% plot(T1,P1,'-ko');
% plot(T2_1,P2,'-o');
% hold off;



%disp(T2_1(zi));

T2 = T2_1(zi);

v2 = MolarVolumeVDWEOS(T2,P2,a,b);

v2_1(zi) = v2;

vmid = MolarVolumeVDWEOS(T2,P1,a,b);

%delH_T = H_func_varT(T2,vmid) - H_func_varT(T1,v1);
%delH_V = H_func_varV(T2,v2) - H_func_varV(T2,vmid);

%delH_T = (H_func_varT(T2,v1) - H_func_varT(T1,v1) );
%delH_V = ( H_func_varV(T2,v2) - H_func_varV(T2,v1) );


delH_V = ( H_func_varV(T1,v2) - H_func_varV(T1,v1) );
delH_T = (H_func_varT(T2,v2) - H_func_varT(T1,v2) );

delH_isentropic = delH_T + delH_V;


eta = 0.80;
%eta = 1.00;

W_12 = delH_isentropic/eta;




T2g_arr = [T2-150:0.1:T2+450];

%T2g_arr = [T1:0.1:T1+400];

prev = 1.0;
T2_2(zi) = T2;


%figure(1);clf

for i=1:numel(T2g_arr)
   
    T2g = T2g_arr(i);
    
    %S_2 = S_func(T2g);
    
    
    lhs = W_12;
    
    v2 = MolarVolumeVDWEOS(T2g,P2,a,b);
    
    vmid = MolarVolumeVDWEOS(T2g,P1,a,b);
    
    rhs = delH_func(T1,v1,T2g,v2,vmid);
    
    %plot(T2g,(lhs-rhs));
    %plot(T2g,P2);
    %hold on;
    
    if(i>1 && prev*(lhs-rhs)<0.0)
        
         T2_2(zi) = T2g_arr(i-1) + ( T2g_arr(i)-T2g_arr(i-1) )/( (lhs-rhs) - prev )*( 0.0-prev ) ;
         
         
         %rho2 = W/v2;
         %as2_arr(i) = 1.0/sqrt(rho2*ks(j));
         break;
     end
     
     prev = lhs - rhs;
    
end

%disp(T2_2(zi));


end

figure(1);clf
plot(T2_1,P2_arr,'-o');
hold on;
plot(T2_2,P2_arr,'-k');

figure(2);clf
plot(v2_1,P2_arr);
end

function vm = MolarVolumeVDWEOS(T0,P0,Am0,Bm0)

R = 8.3144621; 
a10 = -(Bm0 + R*T0/P0);
a20 =  Am0/P0;
a30 = - Am0*Bm0/P0;

[x10,x20,x30] = cubsolve(a10,a20,a30);

vm = 0.0;

if(isreal(x10))
   if(vm<x10)
       vm = x10;
   end
end
if(isreal(x20))
   if(vm<x20)
       vm = x20;
   end
end
if(isreal(x30))
   if(vm<x30)
       vm = x30;
   end
end

% 
% if(isreal(x10) && x10>0.0)
%    if(vm>x10)
%        vm = x10;
%    end
% end
% if(isreal(x20) && x20>0.0)
%    if(vm>x20)
%        vm = x20;
%    end
% end
% if(isreal(x30) && x30>0.0)
%    if(vm>x30)
%        vm = x30;
%    end
% end


%vm = realmax(realmax(x10,x20),x30);
end
function [x1,x2,x3]=cubsolve(a,b,c)

Q=(a*a-3.0*b)/9.0;
R=(2.0*a*a*a-9.0*a*b+27.0*c)/54.0;
R2=R*R;
Q3=Q*Q*Q;
if (R2<Q3) %Three real roots
thita=acos(R/sqrt(Q3));
SQ=sqrt(Q);
x1=-2.0*SQ*cos(thita/3)-a/3.0;
x2=-2.0*SQ*cos((thita+2*pi)/3)-a/3.0;
x3=-2.0*SQ*cos((thita-2*pi)/3)-a/3.0;
return
end
A=-(abs(R)+sqrt(R2-Q3))^(1/3);
if (R < 0)
A=-A; %If we used A=A*sign(R) then if R=0 => sign(R)=0 => A=0 which is wrong
end
if (A==0)
B=0;
else
B=Q/A;
end
AB=A+B;
x1=AB-a/3.0;
comp=1i*sqrt(3)*(A-B)/2.0;
x2=-0.5*AB-a/3.0+comp;
x3=-0.5*AB-a/3.0-comp;

end


