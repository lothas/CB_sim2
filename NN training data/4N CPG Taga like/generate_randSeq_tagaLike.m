function [seq] = generate_randSeq_tagaLike(N,seqOrder)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% % % % CPG parameters:

tau_min = 0.02;     tau_max = 0.25;
tau = (tau_max-tau_min).*rand(1,N) + tau_min;

% b_min = 0.2;     b_max = 10;
b_min = 0.2;     b_max = 2.5;
b = (b_max-b_min).*rand(1,N) + b_min;

c_hip_min = 0;     c_hip_max = 8;
c_hip = (c_hip_max-c_hip_min).*rand(2,N) + c_hip_min;
c_ankle_min = 0;     c_ankle_max = 20;
c_ankle = (c_ankle_max-c_ankle_min).*rand(2,N) + c_ankle_min;
c = [c_ankle;c_hip];

% W_min = 0;     W_max = 10;
W_min = 0;     W_max = 5;
W = (W_max-W_min).*rand(4,N) + W_min;

ks_tau_min = -10;     ks_tau_max = 10;
ks_tau = (ks_tau_max-ks_tau_min).*rand(1,N) + ks_tau_min;

ks_c_hip_min = -0.8;     ks_c_hip_max = 0.8;
ks_c_hip = (ks_c_hip_max-ks_c_hip_min).*rand(2,N) + ks_c_hip_min;
ks_c_ankle_min = -2;     ks_c_ankle_max = 2;
ks_c_ankle = (ks_c_ankle_max-ks_c_ankle_min).*rand(2,N) + ks_c_ankle_min;
ks_c = [ks_c_ankle;ks_c_hip];

seq = [tau;b;c;W;ks_tau;ks_c];

end

