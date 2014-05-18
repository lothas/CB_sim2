function [ Genome ] = CheckGenome( GA, Genome )
% Checks a gene sequence to verify that it satisfies
% all the min/max conditions

    % Check the controller
    Con = GA.Sim.Con;

    % Make sure that a torque pulse doesn't end after the next neuron spike
    for ge = 1:Con.nPulses
        if any(ge == Con.ExtPulses)
            % External pulses con finish at any point
            continue;
        end
        
        End=Genome(3*ge+1)+Genome(3*ge+2);
        if End>=1
            % Shorten both the offset and the duration
            % so that End will be slightly smaller than 1
            Genome(3*ge+1)=Genome(3*ge+1)/(End*1.001);
            Genome(3*ge+2)=Genome(3*ge+2)/(End*1.001);
        end
    end
    
    % Make sure that the max sum of torques remains below the max torque
    CPG=Controller();
    CPG.NumActJoints=NumActJoints;
    CPG=CPG.LoadParameters(Genome);
    [ThisTime,ThisTorques]=CPG.GetTorqueSig(0.005); %#ok<ASGLU>
    N=size(ThisTorques,2);
    Torque=zeros(2,N);
    switch NumActJoints
        case 1
            % Use only hip joint
            for ti=1:NumTorques
                Torque(2,:)=Torque(2,:)+ThisTorques(ti,:);
            end
            
            MaxT2=max(abs(Torque(2,:)));
            MaxAllowed=2*max(abs(TorqueMin(1)),abs(TorqueMax(1)));
            if MaxT2>MaxAllowed
                for ge=1:NumTorques
                    Genome(3*ge)=Genome(3*ge)/MaxT2*MaxAllowed;
                end
            end
        case 2
            % Use ankle and hip joints
            for ti=1:NumTorques/2
                Torque(1,:)=Torque(1,:)+ThisTorques(ti,:);
                Torque(2,:)=Torque(2,:)+ThisTorques(NumTorques/2+ti,:);
            end
            
            MaxT1=max(abs(Torque(1,:)));
            MaxT2=max(abs(Torque(2,:)));
            MaxAllowed=2*max(abs(TorqueMin(1)),abs(TorqueMax(1)));
            
            if MaxT1>MaxAllowed
                for ge=1:NumTorques/2
                    Genome(3*ge)=Genome(3*ge)/MaxT1*MaxAllowed;
                end
            end
                
            if MaxT2>MaxAllowed
                for ge=NumTorques/2+1:NumTorques
                    Genome(3*ge)=Genome(3*ge)/MaxT2*MaxAllowed;
                end
            end
    end
end

