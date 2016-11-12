clear all
close all
clc
l    =  20;       %Wall Length 
h    =  5;        %Wall Height 
Uint =  10;      %Initial Velocity
 
Psi_old = zeros(51,201);% Defining a Matrix M whose values will be replaced  
                      % by stream function (SF) values    
u       =   zeros(51,201);   % Matrix u for velocity in x-direction
v       =   zeros(51,201);   % Matrix v for velocity in y-direction
V       =   zeros(51,201);   % Matrix V for Total velocity 
Cp      =   zeros(51,201);   % Matrix for Pressure coefficient values.
P       =   zeros(51,201);   % Matrix for Pressure value at each point
 
for i=2:51         %Assigning inlet Psi values 
Psi_old(i,1)=Psi_old(i-1,1)+10*0.1;
end
 
%Assigning initial boundary conditions
Psi_old(1,2:90)             =   0; 
Psi_old(11,91:110)          =   Psi_old(11,1);
Psi_old(1,111:201)          =   0;
Psi_old(51,2:90)            =   Psi_old(51,1);
Psi_old(31,91:130)          =   Psi_old(31,1);
Psi_old(51,131:201)         =   Psi_old(51,1);
 
Psi_new  =      Psi_old; %Assigning the Psi values to the new Psi matrix
Abserr   =      1;       %Initial error value 
Abstolr  =      1e-7;    %Max. Iteration deviation allowed
 
%The given diagram is divided into 4 sections.
while max(abs(Abserr)) > Abstolr

% For Streamlines of area 1
for i=2:50
for j=2:90
Psi_new(i,j)=0.25*(Psi_old(i-1,j)+Psi_old(i+1,j)+Psi_old(i,j-1)+Psi_old(i,j+1));
end
end

% For Streamlines of area 2
for i=11:30
for j=91:100
Psi_new(i,j)=0.25*(Psi_old(i-1,j)+Psi_old(i+1,j)+Psi_old(i,j-1)+Psi_old(i,j+1));
end
end

% for Streamlines of area 3
for i=2:30
for j=101:130
Psi_new(i,j)=0.25*(Psi_old(i-1,j)+Psi_old(i+1,j)+Psi_old(i,j-1)+Psi_old(i,j+1));
end
end
 
% For Streamlines of area 4
for i=2:50
for j=131:200
Psi_new(i,j)=0.25*(Psi_old(i-1,j)+Psi_old(i+1,j)+Psi_old(i,j-1)+Psi_old(i,j+1));
end
end
 
Abserr = Psi_new(:)-Psi_old(:); % Error calculation
Psi_old=Psi_new;
end
 
%Velocity and coefficient of pressure calculations
for i=2:50
for j=2:90
u(i,j)       =   (Psi_new(i,j)-Psi_new(i-1,j))/0.1;
v(i,j)       =   (Psi_new(i,j+1)-Psi_new(i,j))/0.1;
V(i,j)       =   (((u(i,j))^2)+(v(i,j))^2)^0.5;
Cp(i,j)      =   1-((((u(i,j))^2)+(v(i,j))^2)/100);
P(i,j)       =   0.5*1000*V(i,j)^2;
end
end
 
for i=11:30
for j=91:100
u(i,j)       =   (Psi_new(i,j)-Psi_new(i-1,j))/0.1;
v(i,j)       =   (Psi_new(i,j+1)-Psi_new(i,j))/0.1;
V(i,j)       =   (((u(i,j))^2)+(v(i,j))^2)^0.5;
Cp(i,j)      =   1-((((u(i,j))^2)+(v(i,j))^2)/100);
P(i,j)       =   0.5*1000*V(i,j)^2;
end
end
 
for i=2:30
for j=101:130
u(i,j)       =   (Psi_new(i,j)-Psi_new(i-1,j))/0.1;
v(i,j)       =   (Psi_new(i,j+1)-Psi_new(i,j))/0.1;
V(i,j)       =   (((u(i,j))^2)+(v(i,j))^2)^0.5;
Cp(i,j)      =   1-((((u(i,j))^2)+(v(i,j))^2)/100);
P(i,j)       =   0.5*1000*V(i,j)^2;
end
end
 
for i=2:50
for j=131:200
u(i,j)       =   (Psi_new(i,j)-Psi_new(i-1,j))/0.1;
v(i,j)       =   (Psi_new(i,j+1)-Psi_new(i,j))/0.1;
V(i,j)       =   (((u(i,j))^2)+(v(i,j))^2)^0.5;
Cp(i,j)      =   1-((((u(i,j))^2)+(v(i,j))^2)/100);
P(i,j)       =   0.5*1000*V(i,j)^2;
end
end

% Plotting the graph with MESH
figure(1);
mesh(Psi_new);
xlabel('X Position');
ylabel('Y Position');
zlabel('Z Position');

% Plotting the graph with contour streamlines
figure(2)
contour(Psi_new,5000);
xlabel('Length ');
ylabel('Height');

% Plotting graph contour of resultant velocity
figure(3)
hold on
contour(V,5000);
colormap jet
xlabel('Length');
ylabel('Height')

% Plotting graph with contour of pressure coefficient
figure(4)
hold on
contour(Cp,5000);
xlabel('Length');
ylabel('Height');


