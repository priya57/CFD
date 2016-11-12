clear all
close all
%Matlab code for Explicit
Up      = 0.1;           %Velocity of Upper Plate
Dg      = 1280.84;       %Density of glycerin
Vi      = 0.8943;        %viscosity of glycerin
dp      = 0.01;          %Distance between grid points
d       = 0.01;          %Distance between plates
l       = 1;             %Length of the plate
D       = 1:1:101;       
Ymin    = 1;             %Minimun y coordiante point
Ymax    = 100;           %Maximum Y coordinate point
M       = zeros(Ymax+1); %Assigning null matrix
Mi      = zeros(Ymax+1);
M(1)    = 0;         %Initial speed at bottom of plate for first time step
M(101)  = 0.1;       %Value of speed at top of the plate for initial Time step
Mi(101) = 0.1;       %Initial speed at bottom of plate for Second time step
Re      = (Up*Dg*d)/Vi;     %formula of Reynolds Number
Ts      = Re*(dp^2)/2;       %calculating timestep
AV      = Ts/(Re*dp^2);      %calculating alpha value
SN      = 100;               %Stopping Number
CT      = 0;                 %Loop counter
        
for T = 1:4500 % Time step for iteration
for i = Ymin+1:Ymax %Loop for 2nd order Equation
Mi(i) = M(i)+AV*(M(i+1)-2*M(i)+M(i-1)); M(i) = Mi(i);
end
if (T == SN)
figure(1)
hold on
plot(Mi,D); % velocity verses height
CT = CT+1;
SN=SN+100;
end
end
title('graph using Explicit')
ylabel('height(m)');
xlabel('Velocity(m/s)')
ylim([0 100]);