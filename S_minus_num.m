load('GA_11_24_10_25.mat');
load('Gen4CL');
ID = 4;
ii=230;
Slope = Data.Slopes(ii);
Sim = deepcopy(GA.Sim);
Sim.Con.FBType = 2;
Sim.Env = Terrain(0,Slope,Slope);
Sim = GA.Gen.Decode(Sim, GA.Seqs(ID,:,GA.Progress));

PMeps = 5e-6;
IC = Data.IC(:,ii);
for tr = 1:1:5
     vec=[0;0;0;0;0];
     vec(tr)=1;
     delta0 = vec*5e-6;
     ICerr = IC + delta0;%IC + delta0';
     [X,State_OUT,f] = take_step(Sim,ICerr,Slope/180*pi);
     X
     
     T(:,tr)= State_OUT.T(end,:)-Data.LCt{1,ii}(end,:);
     tstep = 1e-4; % smaller simulation time
     if T(:,tr)>0
         % Run more accurate simulation until t = T* (normal period)
         t_bef = find(State_OUT.T<Data.LCt{1,ii}(end,:),1,'last');
         X0 = State_OUT.X(t_bef,:);
         T0 = State_OUT.T(t_bef,:);
         Sim.Mod.LegShift = 0;
     else
         % Run more accurate simulation until t = T* (normal period, without impact!)
         X0 = State_OUT.X(end,:);
         T0 = State_OUT.T(end,:);
         Sim.Mod.LegShift = 2*Sim.Mod.Clearance; % shorter leg to avoid the ground
     end
     simT = {T0,tstep,Data.LCt{1,ii}(end,:)+tstep};
     [X,State_OUT,f] = take_step(Sim,X0,Slope/180*pi,simT,1);
     % Still need to post-process to find the desired point in State_OUT
     
     Del_gal(:,tr) = State_OUT.X(end,:)-Data.LCx{1,ii}(end,:)
     Del_plus(:,tr)= X-IC;
     P_xnum(:,tr)= X-IC;
end
eig(P_xnum/(5e-6))
Data.EigV(:,ii)

C1 = [1,0,0,0,0;0,0,1,0,0;0,0,0,1,0;0,0,0,0,1];
H1 = [1,0,0,0;-1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1];
C1*(Del_plus/(5e-6))*H1
C1*(Del_gal/(5e-6))*H1
S_hit_num = C1*Del_plus/(5e-6)*H1/(C1*Del_gal/(5e-6)*H1)

load('EIGVal_S_minus_Gen_4');

eig(C1*Del_plus/(5e-6)*H1/(C1*Del_gal/(5e-6)*H1))

EIGVal_S_minus(:,ii)

eig(Del_plus/(5e-6)/(Del_gal/(5e-6)))

delta0 = PMeps*[1,0,0,0,0];
