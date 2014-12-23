syms  th1_td th2_td dth1_td  dth2_td phi%  m mh l g betta eps
syms  th1 th2 dth1_td dth2_td


%    load('Gen144');seq = 1.4746932;
filename = 'Gen4CL';
fnext = [filename,'.mat'];
load(fnext);

sm_filename = [filename,'SaltMat.mat'];
if exist(sm_filename,'file') == 2
    % Just plot
    Plot_EV_CB (AnEigV,NumEigV,Slopes,Alpha)
else
    % Calculate save and plot
    q_td_gal = [th1_td; th2_td; dth1_td; dth2_td;phi]; 

    % Model 
    m = 3; % leg mass
    mh = 10; % hip mass
    l = 1; % leg length
    g = 9.81;
    betta = mh/m;
    I = 0;%1/3*m*l^2;

    M0 = [ 5/4*m*l^2+mh*l^2+I ,      -m*l^2/2;
              -m*l^2/2 ,              1/4*m*l^2+I];
    G0 = [-(mh+m*3/2)*l*g,     0   ;
                   0,          0.5*m*l*g ];
    E = [1,-1;
         0, 1] ;        
    Gtot = -inv(M0)*(G0);
    Bb = M0\E;
    A = [    0,       0,      1, 0, 0;
             0,       0,      0, 1, 0;
         Gtot(1,1), Gtot(1,2),0, 0, 0;
         Gtot(2,1), Gtot(2,2),0, 0, 0;
             0,       0,      0, 0, 0];
    B = [  0,      0,     0;
           0,      0,     0;
        Bb(1,1), Bb(1,2), 0;
        Bb(2,1), Bb(2,2), 0;
           0,      0,     1];
    % Impact 
    P = [         1-4*I/m/l^2,                  0 ;
          1-(4+4*betta)*cos(th1_td-th2_td)-4*I/m/l^2,   1-4*I/m/l^2];
    Q = [2*cos(th1_td-th2_td), -1-4*I/m/l^2;
         2*cos(th1_td-th2_td)-5-4*betta-4*I/m/l^2, 2*cos(th1_td-th2_td)-1-4*I/m/l^2];

    Z = Q\P;

    %Impact Law
    G = [0,1,   0,     0,     0;
         1,0,   0,     0,     0;
         0,0, Z(1,1),Z(1,2),  0;
         0,0, Z(2,1),Z(2,2),  0;
         0,0,   0,     0,     1];


    for ii = 1:1:length(Data.Slopes)%6:1:200%
        if Data.Period(ii,1)==1
           EIGVal(:,ii) = SM_Period_1(Data,ii,A,B,G,m,mh);
    %     elseif Data.Period(ii,1)==2
    %        EIGVal(:,ii) = SM_Period_2(Data);
        end
    end
    %%

    AnEigV = ConnectVectors(EIGVal(:,1:1:end));
    NumEigV = ConnectVectors(Data.EigV);
    Alpha = (Data.IC(1,:)-Data.IC(2,:))/2;
    Slopes = Data.Slopes;

    save(sm_filename,'AnEigV','NumEigV','Alpha','Slopes');

    Plot_EV_CB (AnEigV,NumEigV,Slopes,Alpha)
end


