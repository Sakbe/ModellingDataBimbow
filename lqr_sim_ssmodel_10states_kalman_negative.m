%% LQR WITH SS MODEL OF FLAT TOP
close all
clear
clc

%%import data
load('C:\Users\Il_Ciancio\Desktop\ISSTOK\MATLAB\ISTTOK_Data_collecting-master\shot_46628.mat')
data_3 = load('C:\Users\Il_Ciancio\Desktop\ISSTOK\MATLAB\ISTTOK_Data_collecting-master\shot_45993.mat');
temp = 350;
time_sec = (1e-6*data.time(11:end-temp));
data_2 = data_3.data;
%% 46628
Ip_magn = double(data.Ip_magn(11:end-temp));
Rc = double(data.r0_corr(11:end-temp));
Zc = double(data.z0_corr(11:end-temp));
Primary = double(data.prim(11:end-temp));
Horizontal = double(data.hor(11:end-temp));
Vertical = double(data.vert(11:end-temp));
Ts =1e-4;
% plot(time_sec, Ip_magn);
ts1 = find(time_sec >= 0.0762);
te1 = find(time_sec >= 0.0976);

Ip_magn_1 = Ip_magn(ts1(1):te1(1));
Rc_1 = Rc(ts1(1):te1(1));
Zc_1 = Zc(ts1(1):te1(1));
Primary_1 = Primary(ts1(1):te1(1));
Horizontal_1 = Horizontal(ts1(1):te1(1));
Vertical_1 = Vertical(ts1(1):te1(1));

dt1 = time_sec(1:length(ts1(1):te1(1)));

%% 45993
time_sec_2 = (1e-6*data_2.time(11:end-temp));
Ip_magn_2 = double(data_2.Ip_magn(11:end-temp));
Rc_2 = double(data_2.r0_corr(11:end-temp));
Zc_2 = double(data_2.z0_corr(11:end-temp));
Primary_2 = double(data_2.prim(11:end-temp));
Horizontal_2 = double(data_2.hor(11:end-temp));
Vertical_2 = double(data_2.vert(11:end-temp));
ts7 = find(time_sec >= 0.3552);
te7 = find(time_sec >= 0.3788);

Ip_magn_7 = Ip_magn_2(ts7(1):te7(1));
Rc_7 = Rc_2(ts7(1):te7(1));
Zc_7 = Zc_2(ts7(1):te7(1));
Primary_7 = Primary_2(ts7(1):te7(1));
Horizontal_7 = Horizontal_2(ts7(1):te7(1));
Vertical_7 = Vertical_2(ts7(1):te7(1));
dt7 = time_sec_2(1:length(ts7(1):te7(1)));
%% NEGATIVE model without D
states = 10;
data = load('ss10_states_negative.mat');
ss_model = data.ss10_states_negative;
% ss_model = data.ss_10states_contnuos;
input_error1 = 0.1*ones(states,1);
input_error2 = 0.1*ones(states,1);
B_error = [ss_model.B, input_error1, input_error2];
C_eye = eye(states);
C_mod = [ss_model.C; C_eye];
D_mod = zeros(states+2,2+2);
x0 = ss_model.x0;
%% CONTROL AND OBSERVER
fprintf('normal system\n');
if(rank(ctrb(ss_model.A, ss_model.B)) == length(ss_model.A))
    fprintf('fully controllable system\n');
end

if(rank(obsv(ss_model.A, ss_model.C)) == length(ss_model.A))
    fprintf('fully osservable system\n');
end


%% LQR positive
%old
% R = (ss_model.C*ss_model.C');
% R1 = 1e2*eye(2);
% R1(1,1) = 0.7e-4; 
% R1(2,2) = 0.34e-3;
% Q=zeros(states,states);
% for i = 1: 1: 10 
%     Q(i,i) = 1e4; 
% end
% Q2=diag([7000,%magin in z at 100
%          20, 
%          0.001, 
%          100, 
%          800, 
%          1, 
%          1, 
%          1, 
%          1, 
%          1]);
% end old

R = (ss_model.C*ss_model.C');

R1 = 1e2*eye(2);
R1(1,1) = 1e-6; 
R1(2,2) = 1e-6;
%slow control
R2 = 10*R1;
R3 = 100*R1;
%end slow control
Q=zeros(states,states);
for i = 1: 1: 10 
    Q(i,i) = 1e4; 
end
Q2=diag([0,%magin in z at 100
         0, 
         0.001, 
         100, 
         1, 
         1, 
         1, 
         0, 
         0, 
         0]);
[K1_LQR,S1_LQR,E1_LQR] = dlqr(ss_model.A,ss_model.B,Q2, R1);
[K2_LQR,S2_LQR,E2_LQR] = dlqr(ss_model.A,ss_model.B,Q2, R2);
[K3_LQR,S3_LQR,E3_LQR] = dlqr(ss_model.A,ss_model.B,Q2, R3);
w_var = 1e-3;
v_var = 1e-4;
bar = rscale_corona(ss_model,K1_LQR);
bar_2 = rscale_corona(ss_model,K2_LQR);
bar_3 = rscale_corona(ss_model,K3_LQR);
time_simulation = 0.03;
r_sp = 0.025;
z_sp = 0.0;
% bar = bar;
SYS = ss(ss_model.A, B_error, ss_model.C, zeros(2,4), Ts);
[kest, L] = kalman(SYS, w_var*eye(2), v_var*eye(2));
[KEST,L,P,Mx,Z,My]  = kalman(SYS, w_var*eye(2), v_var*eye(2));
Cf = eye(size(ss_model.A,1));
Df = zeros(size(ss_model.A,1),size(ss_model.C,1)+2);

sim('lqr_ssmodel_10states_kalmanretro_negative.slx')