clc
close all
clear all

syms x y thetax thetay l
syms dx dy dthetax dthetay dl 
syms ux uy ul

s=[dthetax; dthetay; dx; dy; dl];

%%State X=[theta_X,theta_Y,X,Y,L,dtheta_X,dtheta_Y,dx,dy,dl];

% Absolute position of the mass 
% 
% xm= x + l*sin(thetax)*cos(thetay);
% ym= y + l*sin(thetay);
% zm= -l*cos(thetax)*cos(thetay);


% mass components
m=10; % load [kg]
mc=10; % mass of the crane [kg]
mw=20 ;% mass of the structure [kg]

Mx= mc + mw;
My= mc;
Ml=m;
g=9.81; 
%friction coeff
mu= 0.001; % dimensionless
cp=0.7;




% Inertial Matrix
%syms m l thetax thetay My Mx Ml

m11=[m*(l^2)*cos(thetay)^2 0 ; 0 m*l^2];

m12=[m*l*cos(thetax)*cos(thetay) 0 0;
    -m*l*sin(thetax)*sin(thetay) m*l*cos(thetay) 0];

m21=[m*l*cos(thetax)*cos(thetay) -m*l*sin(thetax)*sin(thetay);
    0 m*l*cos(thetay);
    0 0];
m22=[Mx + m 0 m*sin(thetax)*cos(thetay);
    0 My+m m*l*sin(thetay);
    m*sin(thetax)*cos(thetay) m*sin(thetay) Ml+m];
M=[m11 m12 ;
    m21 m22];
%Minv=inv(M);



c11=[ -m*(l^2)*dthetay*sin(thetay)*cos(thetay) -m*(l^2)*dthetax*sin(thetay)*cos(thetay);
    m*(l^2)*dthetax*sin(thetay)*cos(thetay) m*l*dl];

c12= [0 0 m*l*dthetax*cos(thetay)^2;
    0 0 m*l*dthetay];

c21= [m*dl*cos(thetax)*cos(thetay)-m*l*dthetax*sin(thetax)*cos(thetay)-m*l*dthetay*cos(thetax)*sin(thetay) -m*dl*sin(thetax)*sin(thetay)-m*l*dthetax*cos(thetax)*sin(thetay)-m*l*dthetay*sin(thetax)*cos(thetay);
    0 m*dl*cos(thetay)-m*l*dthetay*sin(thetay);
    -m*l*dthetax*cos(thetay)^2 -m*l*dthetay];

c22= [0 0 m*dthetax*cos(thetax)*cos(thetay)-m*dthetay*sin(thetax)*sin(thetay);
    0 0 m*dthetay*cos(thetay);
    0 0 0];

C=[ c11 c12; % (5x5)
    c21 c22];


% Gravity Terms

G1= [m*g*l*sin(thetax)*cos(thetay);
    m*g*l*cos(thetax)*sin(thetay)];
G2=[ 0;
    0;
    -m*g*cos(thetax)*cos(thetay)];

G= [ G1;
    G2];



% Driving force vector 

F= [ 0 ; 0 ; ux ; uy ; ul ];

% friction forces 
%%Remarque : friction solid est propotionelle à dx et à dy
T=[
  cp*(dthetax^2)*(l^3)*cos(thetay)^3 ; cp*(dthetay^2)*(l^3)*sin(thetax)^3; mu*s(3,1) ;mu*s(4,1);cp*dl];
%%%%%%% Dynamic Model %%%%%%%%%%


o =M\(F-C*s -G -T); %output

eq=[dthetax; dthetay; dx; dy; dl; o(1,:); o(2,:); o(3,:); o(4,:); o(5,:)];
variable=[thetax;thetay; x; y; l; dthetax; dthetay; dx; dy; dl];
u=[ux;uy;ul];

%%Computation of the matrix in analytical way

A=jacobian(eq,variable);
B=jacobian(eq,u);
% 
% syst_non=ss((A-BK),BN,Cnew,zeros(size(N,1)));
% H=tf(syst_non);

save('SystemMatrix','A','B');
