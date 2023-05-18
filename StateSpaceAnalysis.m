close all;
clear all;
clc;
% model variables 

%% Motor
% Resistance
Rm = 8.4;
% Current-torque (N-m/A)
kt = 0.042;
% Back-emf constant (V-s/rad)
km = 0.042;
%
%% Rotary Arm
% Mass (kg)
mr = 0.095;
% Total length (m)
r = 0.085;
% Moment of inertia about pivot (kg-m^2)
Jr = mr*r^2/3;
% Equivalent Viscous Damping Coefficient (N-m-s/rad)
br = 1e-3; % damping tuned heuristically to match QUBE-Sero 2 response
%
%% Pendulum Link
% Mass (kg)
mp = 0.024;
% Total length (m)
Lp = 0.129;
% Pendulum center of mass (m)
l = Lp/2;
% Moment of inertia about pivot (kg-m^2)
Jp = mp*Lp^2/3;
% Equivalent Viscous Damping Coefficient (N-m-s/rad)
bp = 5e-5; % damping tuned heuristically to match QUBE-Sero 2 response
% Gravity Constant
g = 9.81;
Jt=Jp*Jr-mp^2*l^2*r^2;
% state variables
syms ddtheta dtheta theta ddalpha dalpha alpha;
%Control voltage and the torque
syms Tau vm;
Tau=km/Rm*(vm-km*dtheta);
%%
% state vectors
u = Tau;
x = [theta, alpha, dtheta, dalpha]; 
y = [theta, alpha];
%%
%matrix forms 
A=[0  0                   1                 0; 
   0 0                    0                 1; 
   0 (1/Jt)*mp^2*l^2*r*g -(1/Jt)*Jp*br     (1/Jt)*mp*l*r*bp;
   0 -(1/Jt)*mp*g*l*Jr   (1/Jt)*mp*l*r*br  -(1/Jt)*Jp*bp];

B=[0; 0; Jp/Jt; -mp*r*l/Jt];
C=[1 0 0 0;
   0 1 0 0];
D=[0 ; 0];

%%
% linearized system
linsys = ss(A, B, C, D);
%%
%observability
ob=obsv(linsys);
rank(ob);
%%
%controllability
cb=ctrb(linsys);
rank(cb);

%%stability
syms s;
eig(A) % asimptotically stable 
%%transfer function

H=C*(inv(s*eye(4)-A))*B;
%%
%impulse response function
h=ilaplace(H);
%%
%Discretization
deltat=0.005; %time period of sampling
dsys=c2d(linsys,deltat);
