function [ K ] = MimicGains( NC, NCorig, Slope )
%MIMICGAINS Calculates the gains needed to mimic NCorig behavior
%   MimicGains updates the controller's gains to mimic the torque signals
%   of NCorig at a certain slope. First the frequency is copied and then
%   the amplitudes are adjusted.
%   NOTE: Slope should be provided in degrees.

Phi = Slope*pi/180;
NC.FBType = 2;
NC.lastPhi = Phi;
NCorig.FBType = 2;
NCorig.lastPhi = Phi;
NCorig = NCorig.Adaptation(Phi);

% Mimic the CPG frequency
kOmega = (NCorig.omega-NC.omega0)/Phi;
if Slope>0
    NC.kOmega_u = kOmega;
else
    NC.kOmega_d = kOmega;
end

% Mimic the output signal
% Try different approaches
K = zeros(6,length(NC.kTorques_u));
fval = zeros(6,1);
for i = 1:3
    for s = 0:1
        [K(2*i+s-1,:),fval(2*i+s-1)] = FindMin(i,s,NC,NCorig,Phi);
    end
end

% Try with nearest neighbor
KNN = zeros(1,length(NC.kTorques_u));
for j = 1:size(NC.OutM,1)
    NNIDs = NearestN(NC,NCorig,j);
    
    for n = 1:size(NNIDs,2)
        KNN(NNIDs(1,n)) = (NCorig.Amp(NNIDs(2,n))-NC.Amp0(NNIDs(1,n)))/Phi;
    end
end
fvalNN = FindGains(KNN,NC,NCorig,Phi);
K = [K;KNN];
fval = [fval;fvalNN];
                
% Output gains:
opt = find(fval==min(fval),1,'first');
K = [kOmega, K(opt,:)];

%%%%%%%%%%%%%%%%%%%%%% Auxiliary Functions %%%%%%%%%%%%%%%%%%%%%%

% Find minimum difference between signals by optimizing K
% with different optimization parameters
    function [K,fval] = FindMin(initg,sep,NC,NCorig,Phi)
        switch initg
            case 1 % Initial guess 0
                K = zeros(size(NC.kTorques_u));
            case 2 % Initial guess: current gains
                if Phi>0
                    K = NC.kTorques_u;
                else
                    K = NC.kTorques_d;
                end
            case 3 % Initial guess: mimic orig gains
                if Phi>0
                    K = NCorig.kTorques_u;
                else
                    K = NCorig.kTorques_d;
                end
        end
        
        options = optimset('MaxIter',2000,'TolFun',1e-3); ...
            % 'PlotFcns',@optimplotfval);
        if sep == 1
            fvals = zeros(1,size(NC.OutM,1));
            for jo = 1:size(NC.OutM,1)
                IDs = find(NC.OutM(jo,:)==1);
                [K(IDs),fvals(jo)] = fminsearch(@FindJointGains,K(IDs),...
                    options,NC,NCorig,Phi,jo,IDs);
            end
            fval = sum(fvals);
        else
            [K,fval] = fminsearch(@FindGains,K,options,NC,NCorig,Phi);
        end
            
    end

% Optimize separately
    function diff = FindJointGains(K,NC,NC2,Phi,jID,kIDs)
        if Phi>0
            NC.kTorques_u(kIDs) = K;
        else
            NC.kTorques_d(kIDs) = K;
        end
        NC = NC.Adaptation(Phi);
        diffs = CompareSigs(NC,NC2,0.001);
        
        diff = diffs(jID);
    end

% Optimize all together
    function diff = FindGains(K,NC,NC2,Phi)
        if Phi>0
            NC.kTorques_u = K;
        else
            NC.kTorques_d = K;
        end
        NC = NC.Adaptation(Phi);
        diffs = CompareSigs(NC,NC2,0.001);
        
        diff = sum(diffs);
    end

% Find nearest neighbors
    function NNIDs = NearestN(NC,NCorig,j)
        kIDs = find(NC.OutM(j,:)==1);
        Np = length(kIDs);
        Dist = zeros(Np);
        for op = 1:Np
            oStart = NC.Offset(kIDs(op));
            oEnd = oStart + NC.Duration(kIDs(op));
            for mp = 1:length(kIDs)
                mStart = NCorig.Offset(kIDs(mp));
                mEnd = mStart + NCorig.Duration(kIDs(mp));
                Dist(op,mp) = abs(oStart-mStart) + abs(oEnd-mEnd);
            end
        end

        % Check all possible connections
        PosC = perms(1:Np);
        totDist = zeros(1,size(PosC,1));
        IDbase = 0:Np:Np*(Np-1);
        for td = 1:size(PosC,1)
            pIDs = IDbase + PosC(td,:);
            totDist(td) = sum(Dist(pIDs));
        end
        NN = find(totDist==min(totDist),1,'first');
        NNIDs = [kIDs;
                 kIDs(PosC(NN,:))];
    end
end

