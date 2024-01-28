clear
close all
clc
%%
syms th1 th2 th3 th4 th5 th6 th7 fx fy fz mx my mz
K=diag([fx,fy,fz,mx,my,mz]);
th=[th1 th2 th3 th4 th5 th6 th7];
J = Jacobian_kuka_iiwa(th);
A=J*J.';
B=K;
A=A/trace(A);
B=B/trace(B);
W=sqrt(trace((A-B)*(A-B).'));
qo=[diff(W,th1);diff(W,th2);diff(W,th3);diff(W,th4);diff(W,th5);...
    diff(W,th6);diff(W,th7)];
qof=symfun(qo,[th1 th2 th3 th4 th5 th6 th7 fx fy fz mx my mz]);
matlabFunction(qof,'File','nullspce_cotribution_qo')






function J = Jacobian_kuka_iiwa(q)
syms theta alph d a theta1 theta2 theta3 theta4 theta5 theta6 theta7;
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
J = [Jv;Jw];
J = subs(J,{theta1,theta2,theta3,theta4,theta5,theta6,theta7},{q(1),q(2),q(3),q(4),q(5),q(6),q(7)});

end
