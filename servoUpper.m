clear all;
close all;
clc;
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

%%
%matrix forms 
Jt=Jp*Jr-mp^2*l^2*r^2;
A =[0,                                       0,                                     1,                                     0;
0,                                       0,                                     0,                                     1;
0, (g*l^2*mp^2*r)/(- l^2*mp^2*r^2 + Jp*Jr),     -(Jp*br)/(- l^2*mp^2*r^2 + Jp*Jr), -(bp*l*mp*r)/(- l^2*mp^2*r^2 + Jp*Jr);
0,    (Jr*g*l*mp)/(- l^2*mp^2*r^2 + Jp*Jr), -(br*l*mp*r)/(- l^2*mp^2*r^2 + Jp*Jr),     -(Jr*bp)/(- l^2*mp^2*r^2 + Jp*Jr)];

B = [0; 0; Jp/Jt; mp*l*r/Jt];
C=[1 0 0 0;
   0 1 0 0];
D=[0 ; 0];

%%Pole placement
syms s
K=place(A,B,[-11,-12,-13,-14]);

%%State observer
A2=transpose(A);
B2=transpose(C);
L=transpose(place(A2,B2,[-15,-16,-17,-18]));
F=A-L*C;
H=B;

%LQR controller
Q=[0.75 0   0   0;
   0    20  0   0;
   0    0   22  0;
   0    0   0   15];

R=[20];

Klqr=lqr(A,B,Q,R);
