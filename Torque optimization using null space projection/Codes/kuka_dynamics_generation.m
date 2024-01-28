clear all; close all; clc;
syms t th1(t) th2(t) th3(t) th4(t) th5(t) th6(t) th7(t)
syms th1_d(t) th2_d(t) th3_d(t) th4_d(t) th5_d(t) th6_d(t) th7_d(t)
syms th1dd(t) th2dd(t) th3dd(t) th4dd(t) th5dd(t) th6dd(t) th7dd(t)
g = 9.81;

[Jv1,Jw1,R1,y1] = Jacobian_kuka_1(th1(t));
[Jv2,Jw2,R2,y2] = Jacobian_kuka_2(th1(t),th2(t));
[Jv3,Jw3,R3,y3] = Jacobian_kuka_3(th1(t),th2(t),th3(t));
[Jv4,Jw4,R4,y4] = Jacobian_kuka_4(th1(t),th2(t),th3(t),th4(t));
[Jv5,Jw5,R5,y5] = Jacobian_kuka_5(th1(t),th2(t),th3(t),th4(t),th5(t));
[Jv6,Jw6,R6,y6] = Jacobian_kuka_6(th1(t),th2(t),th3(t),th4(t),th5(t),th6(t));
[Jv7,Jw7,R7,y7] = Jacobian_kuka_7(th1(t),th2(t),th3(t),th4(t),th5(t),th6(t),th7(t));

%% 

m1 = 3.94781;
m2 = 4.50275;
m3 = 2.45520;
m4 = 2.61155;
m5 = 3.41000;
m6 = 3.38795;
m7 = 0.35432;

%Inertia matrix
I1 = [0.00455,0,0;0,0.00454,0.00001;0,0.00001,0.00029];
I2 = [0.00032,0,0;0,0.0001,0;0,0,0.00042];
I3 = [0.00223,-0.00005,0.00007;-0.00005,0.00219,0.00007;0.00007,0.00007,0.00073];
I4 = [0.03844,0.00088,-0.00112;0.00088,0.01144,-0.00111;-0.00112,-0.00111,0.04988];
I5 = [0.00277,-0.00001,0.00001;-0.00001,0.00284,0;0.00001,0,0.00012];
I6 = [0.0005,-0.00005,-0.00003;-0.00005,0.00281,-0.00004;-0.00003,0.00004,0.00232];
I7 = [0.00795,0.00022,-0.00029;0.00022,0.01089,-0.00029;-0.00029,-0.00029,0.00294];

%Kinetic energy
k1 = 0.5*diff(th1(t), t)*((Jw1.'*R1*I1*R1.'*Jw1) + m1*(Jv1.'*Jv1))*[diff(th1(t), t)];
k2 = 0.5*[diff(th1(t), t),diff(th2(t), t)]*((Jw2.'*R2*I2*R2.'*Jw2) + m2*(Jv2.'*Jv2))*[diff(th1(t), t);diff(th2(t), t)];
k3 = 0.5*[diff(th1(t), t),diff(th2(t), t),diff(th3(t), t)]*((Jw3.'*R3*I3*R3.'*Jw3) + m3*(Jv3.'*Jv3))*[diff(th1(t), t);diff(th2(t), t);diff(th3(t), t)];
k4 = 0.5*[diff(th1(t), t),diff(th2(t), t),diff(th3(t), t),diff(th4(t), t)]*((Jw4.'*R4*I4*R4.'*Jw4) + m4*(Jv4.'*Jv4))*[diff(th1(t), t);diff(th2(t), t);diff(th3(t), t);diff(th4(t), t)];
k5 = 0.5*[diff(th1(t), t),diff(th2(t), t),diff(th3(t), t),diff(th4(t), t),diff(th5(t), t)]*((Jw5.'*R5*I5*R5.'*Jw5) + m5*(Jv5.'*Jv5))*[diff(th1(t), t);diff(th2(t), t);diff(th3(t), t);diff(th4(t), t);diff(th5(t), t)];
k6 = 0.5*[diff(th1(t), t),diff(th2(t), t),diff(th3(t), t),diff(th4(t), t),diff(th5(t), t),diff(th6(t), t)]*((Jw6.'*R6*I6*R6.'*Jw6) + m6*(Jv6.'*Jv6))*[diff(th1(t), t);diff(th2(t), t);diff(th3(t), t);diff(th4(t), t);diff(th5(t), t);diff(th6(t), t)];
k7 = 0.5*[diff(th1(t), t),diff(th2(t), t),diff(th3(t), t),diff(th4(t), t),diff(th5(t), t),diff(th6(t), t),diff(th7(t), t)]*((Jw7.'*R7*I7*R7.'*Jw7) + m7*(Jv7.'*Jv7))*[diff(th1(t), t);diff(th2(t), t);diff(th3(t), t);diff(th4(t), t);diff(th5(t), t);diff(th6(t), t);diff(th7(t), t)];

%% 
% 
KE_T = k1 + k2 + k3 + k4 + k5 + k6 +k7;
P1=m1*g*y1; P2=m2*g*y2; P3=m3*g*y3; P4=m4*g*y4; P5=m5*g*y5; P6=m6*g*y6; P7 = m7*g*y7;
PE_T = P1 + P2 + P3 + P4 + P5 + P6 + P7;
Lg = KE_T - PE_T;
%% 

dLbdth1=(diff(Lg,diff(th1(t), t)));dLbdth2=real(diff(Lg,diff(th2(t), t)));dLbdth3=real(diff(Lg,diff(th3(t), t)));dLbdth4=real(diff(Lg,diff(th4(t), t))); dLbdth5=real(diff(Lg,diff(th5(t), t)));dLbdth6=real(diff(Lg,diff(th6(t), t)));dLbdth7=real(diff(Lg,diff(th7(t), t)));
dbdtodLbdth1=diff(dLbdth1, t);dbdtodLbdth2=diff(dLbdth2, t);dbdtodLbdth3=diff(dLbdth3, t);dbdtodLbdth4=diff(dLbdth4, t);dbdtodLbdth5=diff(dLbdth5, t);dbdtodLbdth6=diff(dLbdth6, t);dbdtodLbdth7=diff(dLbdth7, t);
dLbth1=diff(Lg,th1(t));dLbth2=diff(Lg,th2(t));dLbth3=diff(Lg,th3(t));dLbth4=diff(Lg,th4(t));dLbth5=diff(Lg,th5(t));dLbth6=diff(Lg,th6(t));dLbth7=diff(Lg,th7(t));
%% 

T1 = dbdtodLbdth1-dLbth1;
% M1 = subs(T1,{diff(th1(t),t,t),diff(th2(t),t,t),diff(th3(t),t,t),diff(th4(t),t,t),diff(th5(t),t,t),diff(th6(t),t,t),diff(th7(t),t,t)},{th1dd(t),th2dd(t),th3dd(t),th4dd(t),th5dd(t),th6dd(t),th7dd(t)});
% M1_1 = subs(M1,{diff(th1(t),t),diff(th2(t),t),diff(th3(t),t),diff(th4(t),t),diff(th5(t),t),diff(th6(t),t),diff(th7(t),t),g},{0,0,0,0,0,0,0,0});
% G1 = subs(T1,{diff(th1(t),t),diff(th2(t),t),diff(th3(t),t),diff(th4(t),t),diff(th5(t),t),diff(th6(t),t),diff(th7(t),t)},{0,0,0,0,0,0,0});
% H1 = subs(T1,{diff(th1(t),t,t),diff(th2(t),t,t),diff(th3(t),t,t),diff(th4(t),t,t),diff(th5(t),t,t),diff(th6(t),t,t),diff(th7(t),t,t),g},{0,0,0,0,0,0,0,0});
% %% 
% 
T2 = dbdtodLbdth2-dLbth2;
% M2 = subs(T2,{diff(th1(t),t,t),diff(th2(t),t,t),diff(th3(t),t,t),diff(th4(t),t,t),diff(th5(t),t,t),diff(th6(t),t,t),diff(th7(t),t,t)},{th1dd(t),th2dd(t),th3dd(t),th4dd(t),th5dd(t),th6dd(t),th7dd(t)});
% M2_2 = subs(M2,{diff(th1(t),t),diff(th2(t),t),diff(th3(t),t),diff(th4(t),t),diff(th5(t),t),diff(th6(t),t),diff(th7(t),t),g},{0,0,0,0,0,0,0,0});
% G2 = subs(T2,{diff(th1(t),t),diff(th2(t),t),diff(th3(t),t),diff(th4(t),t),diff(th5(t),t),diff(th6(t),t),diff(th7(t),t)},{0,0,0,0,0,0,0});
% H2 = subs(T2,{diff(th1(t),t,t),diff(th2(t),t,t),diff(th3(t),t,t),diff(th4(t),t,t),diff(th5(t),t,t),diff(th6(t),t,t),diff(th7(t),t,t),g},{0,0,0,0,0,0,0,0});
% 
T3 = dbdtodLbdth3-dLbth3;
% M3 = subs(T3,{diff(th1(t),t,t),diff(th2(t),t,t),diff(th3(t),t,t),diff(th4(t),t,t),diff(th5(t),t,t),diff(th6(t),t,t),diff(th7(t),t,t)},{th1dd(t),th2dd(t),th3dd(t),th4dd(t),th5dd(t),th6dd(t),th7dd(t)});
% M3_3 = subs(M3,{diff(th1(t),t),diff(th2(t),t),diff(th3(t),t),diff(th4(t),t),diff(th5(t),t),diff(th6(t),t),diff(th7(t),t),g},{0,0,0,0,0,0,0,0});
% G3 = subs(T3,{diff(th1(t),t),diff(th2(t),t),diff(th3(t),t),diff(th4(t),t),diff(th5(t),t),diff(th6(t),t),diff(th7(t),t)},{0,0,0,0,0,0,0});
% H3 = subs(T3,{diff(th1(t),t,t),diff(th2(t),t,t),diff(th3(t),t,t),diff(th4(t),t,t),diff(th5(t),t,t),diff(th6(t),t,t),diff(th7(t),t,t),g},{0,0,0,0,0,0,0,0});
% 
T4 = dbdtodLbdth4-dLbth4;
% M4 = subs(T4,{diff(th1(t),t,t),diff(th2(t),t,t),diff(th3(t),t,t),diff(th4(t),t,t),diff(th5(t),t,t),diff(th6(t),t,t),diff(th7(t),t,t)},{th1dd(t),th2dd(t),th3dd(t),th4dd(t),th5dd(t),th6dd(t),th7dd(t)});
% M4_4 = subs(M4,{diff(th1(t),t),diff(th2(t),t),diff(th3(t),t),diff(th4(t),t),diff(th5(t),t),diff(th6(t),t),diff(th7(t),t),g},{0,0,0,0,0,0,0,0});
% G4 = subs(T4,{diff(th1(t),t),diff(th2(t),t),diff(th3(t),t),diff(th4(t),t),diff(th5(t),t),diff(th6(t),t),diff(th7(t),t)},{0,0,0,0,0,0,0});
% H4 = subs(T4,{diff(th1(t),t,t),diff(th2(t),t,t),diff(th3(t),t,t),diff(th4(t),t,t),diff(th5(t),t,t),diff(th6(t),t,t),diff(th7(t),t,t),g},{0,0,0,0,0,0,0,0});
% 
T5 = dbdtodLbdth5-dLbth5;
% M5 = subs(T5,{diff(th1(t),t,t),diff(th2(t),t,t),diff(th3(t),t,t),diff(th4(t),t,t),diff(th5(t),t,t),diff(th6(t),t,t),diff(th7(t),t,t)},{th1dd(t),th2dd(t),th3dd(t),th4dd(t),th5dd(t),th6dd(t),th7dd(t)});
% M5_5 = subs(M5,{diff(th1(t),t),diff(th2(t),t),diff(th3(t),t),diff(th4(t),t),diff(th5(t),t),diff(th6(t),t),diff(th7(t),t),g},{0,0,0,0,0,0,0,0});
% G5 = subs(T5,{diff(th1(t),t),diff(th2(t),t),diff(th3(t),t),diff(th4(t),t),diff(th5(t),t),diff(th6(t),t),diff(th7(t),t)},{0,0,0,0,0,0,0});
% H5 = subs(T5,{diff(th1(t),t,t),diff(th2(t),t,t),diff(th3(t),t,t),diff(th4(t),t,t),diff(th5(t),t,t),diff(th6(t),t,t),diff(th7(t),t,t),g},{0,0,0,0,0,0,0,0});
% 
T6 = dbdtodLbdth6-dLbth6;
% M6 = subs(T6,{diff(th1(t),t,t),diff(th2(t),t,t),diff(th3(t),t,t),diff(th4(t),t,t),diff(th5(t),t,t),diff(th6(t),t,t),diff(th7(t),t,t)},{th1dd(t),th2dd(t),th3dd(t),th4dd(t),th5dd(t),th6dd(t),th7dd(t)});
% M6_6 = subs(M6,{diff(th1(t),t),diff(th2(t),t),diff(th3(t),t),diff(th4(t),t),diff(th5(t),t),diff(th6(t),t),diff(th7(t),t),g},{0,0,0,0,0,0,0,0});
% G6 = subs(T6,{diff(th1(t),t),diff(th2(t),t),diff(th3(t),t),diff(th4(t),t),diff(th5(t),t),diff(th6(t),t),diff(th7(t),t)},{0,0,0,0,0,0,0});
% H6 = subs(T6,{diff(th1(t),t,t),diff(th2(t),t,t),diff(th3(t),t,t),diff(th4(t),t,t),diff(th5(t),t,t),diff(th6(t),t,t),diff(th7(t),t,t),g},{0,0,0,0,0,0,0,0});
% 
T7 = dbdtodLbdth7-dLbth7;
% M7 = subs(T7,{diff(th1(t),t,t),diff(th2(t),t,t),diff(th3(t),t,t),diff(th4(t),t,t),diff(th5(t),t,t),diff(th6(t),t,t),diff(th7(t),t,t)},{th1dd(t),th2dd(t),th3dd(t),th4dd(t),th5dd(t),th6dd(t),th7dd(t)});
% M7_7 = subs(M7,{diff(th1(t),t),diff(th2(t),t),diff(th3(t),t),diff(th4(t),t),diff(th5(t),t),diff(th6(t),t),diff(th7(t),t),g},{0,0,0,0,0,0,0,0});
% G7 = subs(T7,{diff(th1(t),t),diff(th2(t),t),diff(th3(t),t),diff(th4(t),t),diff(th5(t),t),diff(th6(t),t),diff(th7(t),t)},{0,0,0,0,0,0,0});
% H7 = subs(T7,{diff(th1(t),t,t),diff(th2(t),t,t),diff(th3(t),t,t),diff(th4(t),t,t),diff(th5(t),t,t),diff(th6(t),t,t),diff(th7(t),t,t),g},{0,0,0,0,0,0,0,0});
% o=0
%
% [c1,t1] = coeffs(simplify(M1_1),[th1dd(t),th2dd(t),th3dd(t),th4dd(t),th5dd(t),th6dd(t),th7dd(t)]);
% o =1
% [c2,t2] = coeffs(simplify(M2_2),[th1dd(t),th2dd(t),th3dd(t),th4dd(t),th5dd(t),th6dd(t),th7dd(t)]);
% o=2
% [c3,t3] = coeffs(simplify(M3_3),[th1dd(t),th2dd(t),th3dd(t),th4dd(t),th5dd(t),th6dd(t),th7dd(t)]);
% [c4,t4] = coeffs(simplify(M4_4),[th1dd(t),th2dd(t),th3dd(t),th4dd(t),th5dd(t),th6dd(t),th7dd(t)]);
% [c5,t5] = coeffs(simplify(M5_5),[th1dd(t),th2dd(t),th3dd(t),th4dd(t),th5dd(t),th6dd(t),th7dd(t)]);
% o=3
% [c6,t6] = coeffs(simplify(M6_6),[th1dd(t),th2dd(t),th3dd(t),th4dd(t),th5dd(t),th6dd(t),th7dd(t)]);
% [c7,t7] = coeffs(simplify(M7_7),[th1dd(t),th2dd(t),th3dd(t),th4dd(t),th5dd(t),th6dd(t),th7dd(t)]);
% o=4
% 

function [Jv,Jw,R,Y] = Jacobian_kuka_1(theta1)
syms theta alph d a;
M = [cos(theta), -sin(theta)*cos(alph), sin(theta)*sin(alph), a*cos(theta);sin(theta),cos(theta)*cos(alph),-cos(theta)*sin(alph), a*sin(theta); 0,sin(alph),cos(alph),d; 0, 0, 0, 1];
A1 = subs(M,{a,alph,d,theta},{0,-pi/2,0.36,theta1});
T01 = A1;
T02 = T01;
T03 = T02;
T04 = T03;
T05 = T04;
T06 = T05;
T07 = T06;
z0 = [0; 0; 1];
Jv = [diff(T07(1,4),theta1);
    diff(T07(2,4),theta1);
    diff(T07(3,4),theta1)];
Jw = [z0];
R = [T07(1,1),T07(1,2),T07(1,3);
    T07(2,1),T07(2,2),T07(2,3);
    T07(3,1),T07(3,2),T07(3,3)];
Y = T07(3,4);
end

function [Jv,Jw,R,Y] = Jacobian_kuka_2(theta1,theta2)
syms theta alph d a;
M = [cos(theta), -sin(theta)*cos(alph), sin(theta)*sin(alph), a*cos(theta);sin(theta),cos(theta)*cos(alph),-cos(theta)*sin(alph), a*sin(theta); 0,sin(alph),cos(alph),d; 0, 0, 0, 1];
A1 = subs(M,{a,alph,d,theta},{0,-pi/2,0.36,theta1});
A2 = subs(M,{a,alph,d,theta},{0,pi/2,0,theta2});
T01 = A1;
T02 = T01*A2;
T03 = T02;
T04 = T03;
T05 = T04;
T06 = T05;
T07 = T06;
z0 = [0; 0; 1];
z1 = [T01(1,3);T01(2,3);T01(3,3)];
Jv = [diff(T07(1,4),theta1),diff(T07(1,4),theta2);
    diff(T07(2,4),theta1),diff(T07(2,4),theta2);
    diff(T07(3,4),theta1),diff(T07(3,4),theta2)];
Jw = [z0 z1];
R = [T07(1,1),T07(1,2),T07(1,3);
    T07(2,1),T07(2,2),T07(2,3);
    T07(3,1),T07(3,2),T07(3,3)];
Y = T07(3,4);
end

function [Jv,Jw,R,Y] = Jacobian_kuka_3(theta1,theta2,theta3)
syms theta alph d a;
M = [cos(theta), -sin(theta)*cos(alph), sin(theta)*sin(alph), a*cos(theta);sin(theta),cos(theta)*cos(alph),-cos(theta)*sin(alph), a*sin(theta); 0,sin(alph),cos(alph),d; 0, 0, 0, 1];
A1 = subs(M,{a,alph,d,theta},{0,-pi/2,0.36,theta1});
A2 = subs(M,{a,alph,d,theta},{0,pi/2,0,theta2});
A3 = subs(M,{a,alph,d,theta},{0,pi/2,0.42,theta3});
T01 = A1;
T02 = T01*A2;
T03 = T02*A3;
T04 = T03;
T05 = T04;
T06 = T05;
T07 = T06;
z0 = [0; 0; 1];
z1 = [T01(1,3);T01(2,3);T01(3,3)];
z2 = [T02(1,3);T02(2,3);T02(3,3)];
Jv = [diff(T07(1,4),theta1),diff(T07(1,4),theta2),diff(T07(1,4),theta3);
    diff(T07(2,4),theta1),diff(T07(2,4),theta2),diff(T07(2,4),theta3);
    diff(T07(3,4),theta1),diff(T07(3,4),theta2),diff(T07(3,4),theta3)];
Jw = [z0 z1 z2];
R = [T07(1,1),T07(1,2),T07(1,3);
    T07(2,1),T07(2,2),T07(2,3);
    T07(3,1),T07(3,2),T07(3,3)];
Y = T07(3,4);
end

function [Jv,Jw,R,Y] = Jacobian_kuka_4(theta1,theta2,theta3,theta4)
syms theta alph d a;
M = [cos(theta), -sin(theta)*cos(alph), sin(theta)*sin(alph), a*cos(theta);
    sin(theta),cos(theta)*cos(alph),-cos(theta)*sin(alph), a*sin(theta); 0,sin(alph),cos(alph),d; 0, 0, 0, 1];
A1 = subs(M,{a,alph,d,theta},{0,-pi/2,0.36,theta1});
A2 = subs(M,{a,alph,d,theta},{0,pi/2,0,theta2});
A3 = subs(M,{a,alph,d,theta},{0,pi/2,0.42,theta3});
A4 = subs(M,{a,alph,d,theta},{0,-pi/2,0,theta4});
T01 = A1;
T02 = T01*A2;
T03 = T02*A3;
T04 = T03*A4;
T05 = T04;
T06 = T05;
T07 = T06;
z0 = [0; 0; 1];
z1 = [T01(1,3);T01(2,3);T01(3,3)];
z2 = [T02(1,3);T02(2,3);T02(3,3)];
z3 = [T03(1,3);T03(2,3);T03(3,3)];
Jv = [diff(T07(1,4),theta1),diff(T07(1,4),theta2),diff(T07(1,4),theta3),diff(T07(1,4),theta4);
    diff(T07(2,4),theta1),diff(T07(2,4),theta2),diff(T07(2,4),theta3),diff(T07(2,4),theta4);
    diff(T07(3,4),theta1),diff(T07(3,4),theta2),diff(T07(3,4),theta3),diff(T07(3,4),theta4)];
Jw = [z0 z1 z2 z3];
R = [T07(1,1),T07(1,2),T07(1,3);
    T07(2,1),T07(2,2),T07(2,3);
    T07(3,1),T07(3,2),T07(3,3)];
Y = T07(3,4);
end

function [Jv,Jw,R,Y] = Jacobian_kuka_5(theta1,theta2,theta3,theta4,theta5)
syms theta alph d a;
M = [cos(theta), -sin(theta)*cos(alph), sin(theta)*sin(alph), a*cos(theta);sin(theta),cos(theta)*cos(alph),-cos(theta)*sin(alph), a*sin(theta); 0,sin(alph),cos(alph),d; 0, 0, 0, 1];
A1 = subs(M,{a,alph,d,theta},{0,-pi/2,0.36,theta1});
A2 = subs(M,{a,alph,d,theta},{0,pi/2,0,theta2});
A3 = subs(M,{a,alph,d,theta},{0,pi/2,0.42,theta3});
A4 = subs(M,{a,alph,d,theta},{0,-pi/2,0,theta4});
A5 = subs(M,{a,alph,d,theta},{0,-pi/2,0.4,theta5});
T01 = A1;
T02 = T01*A2;
T03 = T02*A3;
T04 = T03*A4;
T05 = T04*A5;
T06 = T05;
T07 = T06;
z0 = [0; 0; 1];
z1 = [T01(1,3);T01(2,3);T01(3,3)];
z2 = [T02(1,3);T02(2,3);T02(3,3)];
z3 = [T03(1,3);T03(2,3);T03(3,3)];
z4 = [T04(1,3);T04(2,3);T04(3,3)];
Jv = [diff(T07(1,4),theta1),diff(T07(1,4),theta2),diff(T07(1,4),theta3),diff(T07(1,4),theta4),diff(T07(1,4),theta5);
    diff(T07(2,4),theta1),diff(T07(2,4),theta2),diff(T07(2,4),theta3),diff(T07(2,4),theta4),diff(T07(2,4),theta5);
    diff(T07(3,4),theta1),diff(T07(3,4),theta2),diff(T07(3,4),theta3),diff(T07(3,4),theta4),diff(T07(3,4),theta5)];
Jw = [z0 z1 z2 z3 z4];
R = [T07(1,1),T07(1,2),T07(1,3);
    T07(2,1),T07(2,2),T07(2,3);
    T07(3,1),T07(3,2),T07(3,3)];
Y = T07(3,4);
end

function [Jv,Jw,R,Y] = Jacobian_kuka_6(theta1,theta2,theta3,theta4,theta5,theta6)
syms theta alph d a;
M = [cos(theta), -sin(theta)*cos(alph), sin(theta)*sin(alph), a*cos(theta);sin(theta),cos(theta)*cos(alph),-cos(theta)*sin(alph), a*sin(theta); 0,sin(alph),cos(alph),d; 0, 0, 0, 1];
A1 = subs(M,{a,alph,d,theta},{0,-pi/2,0.36,theta1});
A2 = subs(M,{a,alph,d,theta},{0,pi/2,0,theta2});
A3 = subs(M,{a,alph,d,theta},{0,pi/2,0.42,theta3});
A4 = subs(M,{a,alph,d,theta},{0,-pi/2,0,theta4});
A5 = subs(M,{a,alph,d,theta},{0,-pi/2,0.4,theta5});
A6 = subs(M,{a,alph,d,theta},{0,pi/2,0,theta6});
T01 = A1;
T02 = T01*A2;
T03 = T02*A3;
T04 = T03*A4;
T05 = T04*A5;
T06 = T05*A6;
T07 = T06;
z0 = [0; 0; 1];
z1 = [T01(1,3);T01(2,3);T01(3,3)];
z2 = [T02(1,3);T02(2,3);T02(3,3)];
z3 = [T03(1,3);T03(2,3);T03(3,3)];
z4 = [T04(1,3);T04(2,3);T04(3,3)];
z5 = [T05(1,3);T05(2,3);T05(3,3)];

Jv = [diff(T07(1,4),theta1),diff(T07(1,4),theta2),diff(T07(1,4),theta3),diff(T07(1,4),theta4),diff(T07(1,4),theta5),diff(T07(1,4),theta6);
    diff(T07(2,4),theta1),diff(T07(2,4),theta2),diff(T07(2,4),theta3),diff(T07(2,4),theta4),diff(T07(2,4),theta5),diff(T07(2,4),theta6);
    diff(T07(3,4),theta1),diff(T07(3,4),theta2),diff(T07(3,4),theta3),diff(T07(3,4),theta4),diff(T07(3,4),theta5),diff(T07(3,4),theta6)];
Jw = [z0 z1 z2 z3 z4 z5];
R = [T07(1,1),T07(1,2),T07(1,3);
    T07(2,1),T07(2,2),T07(2,3);
    T07(3,1),T07(3,2),T07(3,3)];
Y = T07(3,4);
end

function [Jv,Jw,R,Y] = Jacobian_kuka_7(theta1,theta2,theta3,theta4,theta5,theta6,theta7)
syms theta alph d a;
M = [cos(theta), -sin(theta)*cos(alph), sin(theta)*sin(alph), a*cos(theta);sin(theta),cos(theta)*cos(alph),-cos(theta)*sin(alph), a*sin(theta); 0,sin(alph),cos(alph),d; 0, 0, 0, 1];
A1 = subs(M,{a,alph,d,theta},{0,-pi/2,0.36,theta1});
A2 = subs(M,{a,alph,d,theta},{0,pi/2,0,theta2});
A3 = subs(M,{a,alph,d,theta},{0,pi/2,0.42,theta3});
A4 = subs(M,{a,alph,d,theta},{0,-pi/2,0,theta4});
A5 = subs(M,{a,alph,d,theta},{0,-pi/2,0.4,theta5});
A6 = subs(M,{a,alph,d,theta},{0,pi/2,0,theta6});
A7 = subs(M,{a,alph,d,theta},{0,0,0.126,theta7});
T01 = A1;
T02 = T01*A2;
T03 = T02*A3;
T04 = T03*A4;
T05 = T04*A5;
T06 = T05*A6;
T07 = T06*A7;
z0 = [0; 0; 1];
z1 = [T01(1,3);T01(2,3);T01(3,3)];
z2 = [T02(1,3);T02(2,3);T02(3,3)];
z3 = [T03(1,3);T03(2,3);T03(3,3)];
z4 = [T04(1,3);T04(2,3);T04(3,3)];
z5 = [T05(1,3);T05(2,3);T05(3,3)];
z6 = [T06(1,3);T06(2,3);T06(3,3)];
Jv = [diff(T07(1,4),theta1),diff(T07(1,4),theta2),diff(T07(1,4),theta3),diff(T07(1,4),theta4),diff(T07(1,4),theta5),diff(T07(1,4),theta6),diff(T07(1,4),theta7);
    diff(T07(2,4),theta1),diff(T07(2,4),theta2),diff(T07(2,4),theta3),diff(T07(2,4),theta4),diff(T07(2,4),theta5),diff(T07(2,4),theta6),diff(T07(2,4),theta7);
    diff(T07(3,4),theta1),diff(T07(3,4),theta2),diff(T07(3,4),theta3),diff(T07(3,4),theta4),diff(T07(3,4),theta5),diff(T07(3,4),theta6),diff(T07(3,4),theta7)];
Jw = [z0 z1 z2 z3 z4 z5 z6];
R = [T07(1,1),T07(1,2),T07(1,3);
    T07(2,1),T07(2,2),T07(2,3);
    T07(3,1),T07(3,2),T07(3,3)];
Y = T07(3,4);
end

