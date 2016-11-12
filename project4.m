clear all
close all
M_a         =   28.976; %mass of air
Kappa       =   1.4;
Mach_no     =   1.2; %mach number
CFL         =   0.4;
P0          =   100000; %given pressure
T0          =   300;
T           =   0.77639751*T0; 
R           =   8314.472; %universal gas constant
Im          =   300;
J           =   1:1:Im;
Tst         =   1000;
At          =   1.0022;
c_sound     =   sqrt(Kappa*R/M_a*T);
dx          =   2/Im;
dt          =   CFL*dx/c_sound;
A           =   ones(1,Im);

for i=1:Im;
 xdL        =   2*i/Im;
 A(i)       =   At*(1+0.2223*(xdL)^2);
 end
 p          =   ones(Tst,Im)*P0;
 rho        =   p/R/T*M_a;
 u          =  (ones(Tst,Im)*Mach_no*c_sound);
 e          =   ones(Tst,Im)*(p(1,1)/(Kappa-1)/rho(1,1)+0.5*(u(1,1))^2);
 for n=1:(Tst-1)
  for i=2:Im
 rho(n+1,i)     =   -dt/(dx*A(i))*((rho(n,i)*u(n,i)*A(i))-(rho(n,i-1)...,
                    *u(n,i-1)*A(i-1)))+rho(n,i);
 c              =   ((rho(n,i)*(u(n,i))^2 +p(n,i))*A(i))-((rho(n,i-1)...,
                    *(u(n,i-1))^2+p(n,i-1))*A(i-1));
 u(n+1,i)       =   1/rho(n+1,i)*(-dt/(dx*A(i))*(c-p(n,i)*(A(i)-A(i-1)))...,
                    +rho(n,i)*u(n,i));
 d              =   (e(n,i)+p(n,i))*u(n,i)*A(i)-(e(n,i-1)+p(n,i-1))*...,
                    u(n,i-1)*A(i-1);
 e(n+1,i)       =   (-dt/(dx*A(i))*d+rho(n,i)*e(n,i))/rho(n+1,i);
 p(n+1,i)       =   rho(n+1,i)*R*T/M_a;
 end
end

figure(1);plot(J,A,'k',J,-A,'k');title({'Area'});
figure(2);plot(u(Tst,:));title({'Velocity (m/s)'});
figure(3);plot(p(Tst,:));title({'Pressure (Pa)'});
figure(4);plot(e(Tst,:));title({'Energy (J)'});
figure(5);plot(rho(Tst,:));title({'Density(kg/m3'});

