function EV = SM_Period_1(Data,ii,A,B,G,m,mh)
syms  th1_td th2_td dth1_td  dth2_td 
syms  th1 th2 dth1_td dth2_td phi
q_td_gal = [th1_td; th2_td; dth1_td; dth2_td;phi]; 

 
omega = 1/Data.LCt{ii}(end);%seq(1);%
% OMEGA(ii)=omega;
% Phi_impact= Data.LCx{1,ii}(end,end);  
t0 = 0;%(1-Phi_impact)/omega;

Tr_an = diff(Data.LCtorques{1,ii}(:,1)); 
IND_Torq_an = find(Tr_an);
Tr_hip = diff(Data.LCtorques{1,ii}(:,2)); 
IND_Torq_hip = find(Tr_hip);
Torq_vec = [[IND_Torq_an+1;IND_Torq_hip+1] [Data.LCtorques{1,ii}(IND_Torq_an+1,1);Data.LCtorques{1,ii}(IND_Torq_hip+1,2)] [ones(length(IND_Torq_an),1);2*ones(length(IND_Torq_hip),1)]];
IND_event_T = sortrows(Torq_vec);
Time = t0+[Data.LCt{ii}(IND_event_T(:,1));Data.LCt{ii}(end)];
State = zeros(length(IND_event_T(:,1))+1,5); 
State(1:end-1,:) = Data.LCx{1,ii}((IND_event_T(:,1)),:);
State(end,:) = Data.LCx{1,ii}(end,:);
U_vec = zeros(3,length(IND_event_T(:,1))+1);
U_vec(:,1) = [Data.LCtorques{1,ii}(1,1);Data.LCtorques{1,ii}(1,2);omega];

for pp = 2:1:length(IND_event_T(:,1))+1
    U_vec(:,pp) = [Data.LCtorques{1,ii}(IND_event_T(pp-1,1),1);Data.LCtorques{1,ii}(IND_event_T(pp-1,1),2);omega]; 
end

% Saltation Matrix for Torques events
H_k_x = [0,0,0,0,1];
G_x_k_xgal = eye(5);
F_k_gal = zeros(5,length(State)-1);
F_k_kova = zeros(5,length(State)-1);


for ss = 1:1:length(State)-1
    F_k_gal(:,ss) =  A*State(ss,:)'+B*U_vec(:,ss);
    F_k_kova(:,ss) = A*State(ss,:)'+B*U_vec(:,ss+1);
end

Lk = zeros(5,5,length(State)-1);
for kk = 1:1:length(State)-1
    Lk(:,:,kk) = G_x_k_xgal + ( (F_k_kova(:,kk)-F_k_gal(:,kk))*H_k_x )/(H_k_x*F_k_gal(:,kk));
end

L_torq = 1;
for tt=1:1:length(Time)-1
  L = expm( A*(Time(tt+1)-Time(tt)) )*Lk(:,:,tt);
  L_torq = L*L_torq;
end 
%%%%%%%%%%%%%%%%%% For PWL %%%%%%%%%%%%%%%%%%%%%%
% Saltation Matrix for impact events
% q = [th1; th2; dth1_td ;dth2_td]; 
% 
% M = [ 5*l^2*m/4 + l^2*mh,    -l^2*m*cos(th1-th2)/2  ;
%      -l^2*m*cos(th1-th2)/2   ,  l^2*m/4            ];
% H = [  0,          -m*l^2/2*sin(q(1)-q(2))*q(4);
%      m*l^2/2*sin(q(1)-q(2))*q(3),    0];
% Gln = [ -g*l*sin(th1)*((3*m)/2 + mh) ;
%        g*l*m*sin(th2)/2          ];
% G_tot = -inv(M)*(Gln+H*[dth1_td ;dth2_td]); 
% B_tot = M\E;
% Bf = [  0,      0,      0;
%         0,      0,      0;
% B_tot(1,1), B_tot(1,2), 0;
% B_tot(2,1), B_tot(2,2), 0;
%         0,      0,      1];
% Be = double(subs(Bf,[th1,th2,dth1_td,dth2_td],[ State(end,1),State(end,2),State(end,3),State(end,4)])); 
% Ff = [dth1_td;dth2_td;G_tot(1);G_tot(2);0]; 
% Fe = double(subs(Ff,[th1,th2,dth1_td,dth2_td],[ State(end,1),State(end,2),State(end,3),State(end,4)])); 
%Fe+Be*U_vec(:,end);%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Saltation Matrix for Impact events
F_td = A*State(end,:)'+B*U_vec(:,end);
H_td_x = [1,1,0,0,0];
HF_td = H_td_x*F_td;
G_td = G*q_td_gal;

G_x_td_gal = jacobian(G_td,q_td_gal);
L_hit = G_x_td_gal*( eye(5)-F_td*H_td_x/HF_td );
S_hit = double(subs(L_hit,[th1_td,th2_td,dth1_td,dth2_td,phi],[ State(end,1),State(end,2),State(end,3),State(end,4),State(end,5)]));

P = S_hit*L_torq*expm( A*(Time(1)) );
EV = eig(P);
