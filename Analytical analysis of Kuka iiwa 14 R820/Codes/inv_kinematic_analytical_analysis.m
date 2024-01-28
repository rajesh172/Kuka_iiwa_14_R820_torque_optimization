clear all;clc;close all;
%T07 = [1 0 0 0; 0 1 0 0; 0 0 1 1.3060; 0 0 0 1]; for all thetas= 0
T07 = [0 0 -1 -0.5260; 0 1 0 0; 1 0 0 0.78; 0 0 0 1]; % for th4=pi/2
%T07 = transformation_T07(0,0,0,pi/2,0,0,0);
%T07 = transformation_T07(0,0,0,0,0,0,0);
GC2 =1;GC4 =1;GC6 =1;% changing these as the angles are changing

dbs = 0.36;% distance from base to shoulder
dse = 0.42;% distance from shoulder to elbow
dew= 0.4 ; % distance from elbow to wrist
dwf = 0.126;% sitance from wrist to flange

p02 = [0 0 dbs]';
p24 = [0 dse 0]';
p46 = [0 0 dew]';
p67 = [0 0 dwf]';
p07 = [T07(1,4) T07(2,4) T07(3,4)]';
R07 = [T07(1,1) T07(1,2) T07(1,3);
    T07(2,1) T07(2,2) T07(2,3);
    T07(3,1) T07(3,2) T07(3,3)];

p06 = p07-R07*p67;
p26 = p07-p02-R07*p67;
th4_v = real(GC4*acos(((norm(p26))^2-dse^2-dew^2)/(2*dse*dew)));
phi = acos((dse^2 + (norm(p26))^2 -dew^2)/(2*dse*norm(p26)));
th3_v =0;

if norm(cross(p26,[0;0;1]))>0
    th1_v = atan2(p26(2,1),abs(p26(1,1)));
else
    th1_v=0;
end
th2_v = -atan2(sqrt(p26(1,1)^2 + p26(2,1)^2),p26(3,1))+GC4*phi;
%th2_v = -(atan(sqrt(p26(1,1)^2 + p26(2,1)^2)/(p26(3,1)))+GC4*phi);

[T04,T03] = transformation_T04(th1_v,th2_v,th3_v,th4_v); % calling tranformation function to find T04 and T03
psi = 0;

R03 = [T03(1,1) T03(1,2) T03(1,3);
    T03(2,1) T03(2,2) T03(2,3);
    T03(3,1) T03(3,2) T03(3,3)];

R04 = [T04(1,1) T04(1,2) T04(1,3);
    T04(2,1) T04(2,2) T04(2,3);
    T04(3,1) T04(3,2) T04(3,3)];

R0psi = eye(3);

p26_u = p26/(norm(p26)); % unit vector of p26

Skew_u = [0 -p26_u(3,1) p26_u(2,1);
    p26_u(3,1) 0 -p26_u(1,1);
    -p26_u(2,1) p26_u(1,1) 0];

% all auxilary parameters
As = Skew_u*R03;
Bs = -Skew_u*Skew_u*R03;
Cs = p26_u*p26_u'*R03;

R34 = pinv(R03)*R04;
R47 = pinv(R34)*pinv(R03)*R07;
Aw  = (R34'*As'*R07);
Bw  = (R34'*Bs'*R07);
Cw  = (R34'*Cs'*R07);

th1 = th1_v;
th2 = th2_v;
th3 = th3_v;
th4 = th4_v;

th5_num = (Aw(2,3)*sin(psi)+Bw(2,3)*cos(psi)+ Cw(2,3))*GC6;
th5_den = (Aw(1,3)*sin(psi)+Bw(1,3)*cos(psi)+ Cw(1,3))*GC6;
th5 = atan2(-th5_num,-th5_den);

th6= GC6*(acos((Aw(3,3)*sin(psi)+Bw(3,3)*cos(psi)+Cw(3,3))));

th7_num = (Aw(3,2)*sin(psi)+Bw(3,2)*cos(psi)+ Cw(3,2))*GC6;
th7_den = (-Aw(3,1)*sin(psi)-Bw(3,1)*cos(psi)- Cw(3,1))*GC6;
th7 = atan2(-th7_num,-th7_den);

thv = [th1;th2;th3;th4;th5;th6;th7] % all theta


