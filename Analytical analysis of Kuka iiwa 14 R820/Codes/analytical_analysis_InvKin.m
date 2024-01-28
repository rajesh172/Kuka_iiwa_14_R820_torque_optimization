function q = analytical_analysis_InvKin(th1,th2,th3,th4,th5,th6,th7)

%T07 = [0 0 1 0.5260; 0 1 0 0; -1 0 0 0.78; 0 0 0 1];
%T07 = [0 0 -1 -0.5260; 0 1 0 0; 1 0 0 0.78; 0 0 0 1];
%T07 = [0 0 -1 0.2940; 0 1 0 0; 1 0 0 -0.04; 0 0 0 1];
T07 = transformation_T07(th1,th2,th3,th4,th5,th6,th7);
GC2 =1;GC4 =1;GC6 =1;
dbs = 0.36;
dse = 0.42;
dew= 0.4 ;
dwf = 0.126;

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
[T04,T03] = transformation_T04(th1_v,th2_v,th3_v,th4_v);
psi = 0;
R03 = [T03(1,1) T03(1,2) T03(1,3);
    T03(2,1) T03(2,2) T03(2,3);
    T03(3,1) T03(3,2) T03(3,3)];

R04 = [T04(1,1) T04(1,2) T04(1,3);
    T04(2,1) T04(2,2) T04(2,3);
    T04(3,1) T04(3,2) T04(3,3)];
R0psi = eye(3);
p26_u = p26/(norm(p26));
Skew_u = [0 -p26_u(3,1) p26_u(2,1);
    p26_u(3,1) 0 -p26_u(1,1);
    -p26_u(2,1) p26_u(1,1) 0];

As = Skew_u*R03;
Bs = -Skew_u*Skew_u*R03;
Cs = p26_u*p26_u'*R03;
R34 = pinv(R03)*R04;
R47 = pinv(R34)*pinv(R03)*R07;
Aw  = (R34'*As'*R07);
Bw  = (R34'*Bs'*R07);
Cw  = (R34'*Cs'*R07);

% th1_num = (As(2,2)*sin(psi)+Bs(2,2)*cos(psi)+ Cs(2,2))*GC2;
% th1_den = (As(1,2)*sin(psi)+Bs(1,2)*cos(psi)+ Cs(1,2))*GC2;
% th1 = atan2(th1_num,th1_den);
% 
% th2 = GC2*(acos((As(3,2)*sin(psi)+Bs(3,2)*cos(psi)+Cs(3,2))));
% 
% th3_num = (-As(3,3)*sin(psi)-Bs(3,3)*cos(psi)- Cs(3,3))*GC2;
% th3_den = (-As(3,1)*sin(psi)-Bs(3,1)*cos(psi)-Cs(3,1))*GC2;
% th3 = atan2(th3_num,th3_den);
% th4 = th4_v;

th5_v_num = (Aw(2,3)*sin(psi)+Bw(2,3)*cos(psi)+ Cw(2,3))*GC6;
th5_v_den = (Aw(1,3)*sin(psi)+Bw(1,3)*cos(psi)+ Cw(1,3))*GC6;
th5_v = atan2(-th5_v_num,-th5_v_den);

th6_v = GC6*(acos((Aw(3,3)*sin(psi)+Bw(3,3)*cos(psi)+Cw(3,3))));

th7_v_num = (Aw(3,2)*sin(psi)+Bw(3,2)*cos(psi)+ Cw(3,2))*GC6;
th7_v_den = (-Aw(3,1)*sin(psi)-Bw(3,1)*cos(psi)- Cw(3,1))*GC6;
th7_v = atan2(-th7_v_num,-th7_v_den);
thv = [th1_v;th2_v;th3_v;th4_v;th5_v;th6_v;th7_v];
q = thv(:);
% th = [th1;th2;th3;th4;th5_v;th6_v;th7_v]



end

