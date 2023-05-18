close all;
clear all;
clc;
%Linearization of the system
%Variables
syms theta alpha dtheta dalpha ddtheta ddalpha Tau;
syms alpha2 theta2
syms Jr Jp mp l r br g bp;

%Nonlinear system 
eqn1=dtheta==theta2;
eqn2=dalpha==alpha2;
eqn3=(Jr+Jp*sin(alpha)^2)*ddtheta+mp*l*r*cos(alpha)*ddalpha+2*Jp*sin(alpha)*cos(alpha)*theta2*alpha2-mp*l*r*sin(alpha)*alpha2^2==Tau-br*theta2;
eqn4=Jp*ddalpha+mp*l*r*cos(alpha)*ddtheta-Jp*sin(alpha)*cos(alpha)*theta2^2+mp*g*l*sin(alpha)==-bp*alpha2;

%parametric solution of the model for ddtheta and ddalpha
[dtheta2,dalpha2]=solve(eqn3,eqn4,ddtheta,ddalpha);
solution1=subs(dalpha2,[alpha,alpha2,theta,theta2,Tau],[pi,0,0,0,0]);
solution2=subs(dtheta2,[alpha,alpha2,theta,theta2,Tau],[pi,0,0,0,0]);
%Jacobian matrix
J=jacobian([theta2,alpha2,dtheta2, dalpha2],[theta,alpha,theta2,alpha2]);

%Substituing the values of the equilibrium point
%The A matrix of the state space model:
%If alpha=0 it is the stable equilibrium point if alpha=pi, it is the
%unstable equilibrium point
%%
A=subs(J,[theta,alpha,theta2,alpha2,Tau],[0,pi,0,0,0]);

