function [ ] = NumAnalysis( )
% This function runs the numerical analysis for the system.
% This includes obtaining steady state initial conditions for
% a range of slopes, generating the corresponding limit cycle
% data for each slope and calculating the numerical Poincare
% map (external m file)
% Version 0.6 - 9/06/2012
clc

filename='NumericalData.mat';

% Initial conditions for start slope of 0
% InitCond0=[0.18425, -0.18425, -0.72367, -0.59289, 0.88958];
InitCond0=[0.1710746, -0.1710746, -0.6820206, -0.5715793, 0.8916005];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Load existing data or set-up simulations %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist(filename,'file')==2
    load(filename);
    
    samples=length(Slopes); %#ok<NODEF>
else
    Progress=0;
    
%     % Define number of samples (number of slopes that the robot will walk on)
%     samples=401;
% 
%     % Define the slope range
%     Slopes=linspace(-10,10,samples);
    Slopes=[-10:0.01:-9.5 -9.45:0.05:-8.55 -8.54:0.01:-8.36 -8.35:0.05:-3.8 -3.79:0.01:-3.61 ...
        -3.6:0.05:0.8 0.81:0.01:1 1.05:0.05:1.35 1.36:0.01:1.59 1.6:0.05:4.95 4.96:0.01:5.09 ...
        5.1:0.05:5.6 5.61:0.01:5.69 5.7:0.05:9.8 9.81:0.01:10];
    samples=length(Slopes);

    % Arrays to collect simulation values
    EndGame=zeros(samples,1);           % Simulation end result
    LimitCycle_Period=zeros(samples,1); % Limit cycle period
    LimitCycle_ICb=zeros(5*4,samples);  % Steady state initial conditions 
    LimitCycle_Xb=cell(samples,1);       % State-space values for the limit cycle
    LimitCycle_Tb=cell(samples,1);       % Time values for the limit cycle
    LimitCycle_Torquesb=cell(samples,1); % Torques values for the limit cycle
    LimitCycle_eigb=zeros(5*4,samples); % Eigenvalues of the linearized Poincare Map

    % Define the slopes that need to be re-done
    ReDoSlopeID=zeros(samples,1);

    save(filename,'Slopes','EndGame','LimitCycle_Period','LimitCycle_ICb','Progress');
    save(filename,'LimitCycle_Xb','LimitCycle_Tb','LimitCycle_eigb','ReDoSlopeID','-append');
end

if matlabpool('size')==0 && Progress<4
    matlabpool open % Work in parallel to finish faster
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Run simulations for every slope %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Progress==0
    disp('Running first batch of simulations...');
    parfor i=1:samples
        if EndGame(i)~=0
            % Already processed, skip
        else
            [EndGame(i), EndText, LimitCycle_Period(i), IC]=Simulation( InitCond0, 3, 0, Slopes(i), 2, 1, 0 );
            disp(EndText);

            switch LimitCycle_Period(i)
                case 0
                    % Simulation failed, do nothing
                case 1
                    % All initial conditions are the same
                    LimitCycle_ICb(:,i)=[IC;IC;IC;IC];
                case 2
                    % Initial conditions repeat once
                    LimitCycle_ICb(:,i)=[IC(:,1);IC(:,2);IC(:,1);IC(:,2)];
                case 3
                    % Repeat only the first set of IC
                    LimitCycle_ICb(:,i)=[IC(:,1);IC(:,2);IC(:,3);IC(:,1)];
                case 4
                    % All initial conditions are different
                    LimitCycle_ICb(:,i)=[IC(:,1);IC(:,2);IC(:,3);IC(:,4)];
            end

            if EndGame(i)~=4
                % Simulation didn't converge, add it to re-do list
                ReDoSlopeID(i)=i;
            end
        end
    end

    Progress=1;
    save(filename,'EndGame','LimitCycle_Period','LimitCycle_ICb','ReDoSlopeID','Progress','-append');
    disp([10,'']);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Re-do simulations that failed %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Progress<=1
    ReIndex=find(ReDoSlopeID~=0);
%     ReIndex=[337];
    ReSlopes=Slopes(ReIndex);
    N=length(ReSlopes);
    disp(['Running ',num2str(N),' failed simulation(s)...']);

    for i=1:N
        % Find closest slope that didn't fail
        Index=ReIndex(i);
        Index2=Index+find(EndGame(Index+1:end)==4,1,'first');
%         Index2=510;
        if Index2
            % Try re-running the simulation with the closest bigger slope
            InitCond=LimitCycle_ICb(1:5,Index2);
            InitSlope=Slopes(Index2);

            [EndGame(Index), EndText, LimitCycle_Period(Index), IC]=Simulation( InitCond, 3, InitSlope, Slopes(Index), 2, 1, 0 );
            disp(EndText);

            switch LimitCycle_Period(Index)
                case 0
                    if EndGame(Index)==0
                        % All initial conditions are the same
                        LimitCycle_ICb(:,Index)=[IC;IC;IC;IC];
                    else
                        % Simulation failed, do nothing
                    end
                case 1
                    % All initial conditions are the same
                    LimitCycle_ICb(:,Index)=[IC;IC;IC;IC];
                case 2
                    % Initial conditions repeat once
                    LimitCycle_ICb(:,Index)=[IC(:,1);IC(:,2);IC(:,1);IC(:,2)];
                case 3
                    % Repeat only the first set of IC
                    LimitCycle_ICb(:,Index)=[IC(:,1);IC(:,2);IC(:,3);IC(:,1)];
                case 4
                    % All initial conditions are different
                    LimitCycle_ICb(:,Index)=[IC(:,1);IC(:,2);IC(:,3);IC(:,4)];
            end
        end

        if (EndGame(Index)~=4 && EndGame(Index)~=0) && Index>1
            Index2=find(EndGame(1:Index-1)==4,1,'last');
            if Index2
                % Try re-running the simulation with the closest smaller slope
                InitCond=LimitCycle_ICb(1:5,Index2);
                InitSlope=Slopes(Index2);

                [EndGame(Index), EndText, LimitCycle_Period(Index), IC]=Simulation( InitCond, 3, InitSlope, Slopes(Index), 2, 1, 0 );
                disp(EndText);

                switch LimitCycle_Period(Index)
                    case 0
                        if EndGame(Index)==0
                            % All initial conditions are the same
                            LimitCycle_ICb(:,Index)=[IC;IC;IC;IC];
                        else
                            % Simulation failed, do nothing
                        end
                    case 1
                        % All initial conditions are the same
                        LimitCycle_ICb(:,Index)=[IC;IC;IC;IC];
                    case 2
                        % Initial conditions repeat once
                        LimitCycle_ICb(:,Index)=[IC(:,1);IC(:,2);IC(:,1);IC(:,2)];
                    case 3
                        % Repeat only the first set of IC
                        LimitCycle_ICb(:,Index)=[IC(:,1);IC(:,2);IC(:,3);IC(:,1)];
                    case 4
                        % All initial conditions are different
                        LimitCycle_ICb(:,Index)=[IC(:,1);IC(:,2);IC(:,3);IC(:,4)];
                end
            end
        end

        if EndGame(Index)==4
            ReDoSlopeID(Index)=0;
            save(filename,'EndGame','LimitCycle_Period','LimitCycle_ICb','ReDoSlopeID','-append');
        end
    end
    
    if isempty(find(ReDoSlopeID~=0, 1))
        Progress=2;
        save(filename,'Progress','-append');
    else
        disp('Some slopes are still failing');
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Sort out numerical data to avoid "jumping" lines %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Progress<=2
    % Ask user for input to decide whether to continue
    % Show results for ICs and periods
    figure()
    plot(Slopes,LimitCycle_ICb);
    hold on
    plot(Slopes,1+LimitCycle_Period/4);

    choice = questdlg(['All the steady state initial conditions have been calculated.',...
                        10,10,'Would you like to continue?'],...
                        'IC processing done',...
                        'Yes','No','Yes');
    switch choice
        case 'Yes'
            % continue
        case 'No'
            return;
        otherwise
            return;
    end

    disp('Sorting out numerical data to avoid "jumping" lines...');

    [LimitCycle_IC] = SortLines(LimitCycle_ICb);
        
    % By now each row in LimitCycle_IC should hold a single "curve" of
    % initial conditions.
    % Check that the sets of initial conditions 1:5 and 6:10 intersect
    % the poincare section, i.e. they satisfy the geometric condition
    % for impact.
    GeoCond=LimitCycle_IC(1,:)+LimitCycle_IC(2,:)-2*Slopes*pi/180;
    if ~isempty(find(abs(GeoCond)>1e-5,1))
        % This probably means that the sets are mixed.
        % Switching between row 2 and 7 should solve the issue
        Temp=LimitCycle_IC(2,:);
        LimitCycle_IC(2,:)=LimitCycle_IC(7,:);
        LimitCycle_IC(7,:)=Temp;
        disp([10,10,'Warning: Switched lines',10,10]);
    end
        
    Progress=3;
    save(filename,'Slopes','Progress','LimitCycle_IC','LimitCycle_ICb','-append');
    disp([10,'']);
end

function [Lines] = SortLines( Lines )
   % Works great as long as lines don't intersect
   N=size(Lines,2);
   
   % We'll only look at lines 1 and 6 and switch the whole sets
   % 1:5 and 6:10 accordingly. Let row one be always "above" row 6
   for c=1:N
       if Lines(1,c)<Lines(6,c)
           % Need to switch
           Temp=Lines(1:5,c:end);
           Lines(1:5,c:end)=Lines(6:10,c:end);
           Lines(6:10,c:end)=Temp;
       end
   end
   
   % Clear the rest of the rows if there's really no period 4
   if isempty(find(LimitCycle_Period==4,1))
       Lines=Lines(1:10,:);
   end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Generate each slope's limit cycle %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Progress<=3
    disp('Processing Limit Cycles...');
    parfor i=1:samples
        if isempty(LimitCycle_Tb{i})==1
            % Run the simulation for the required number of steps
            % based on the gait's period
            [ Result, ResText, LimitCycle_Tb{i}, LimitCycle_Xb{i}, LimitCycle_Torquesb{i} ]=Simulation( LimitCycle_IC(1:5,i), 0, Slopes(i), Slopes(i), 2, [2,LimitCycle_Period(i)], 0 ); %#ok<ASGLU>
            if Result~=5
                disp(['Failed to process limit cycle for slope: ',num2str(Slopes(i)),char(176)]);
            else
                disp(['Succesfully processed limit cycle for slope: ',num2str(Slopes(i)),char(176)]);
            end
        end
    end
    
    save(filename,'LimitCycle_Xb','LimitCycle_Tb','LimitCycle_Torquesb','-append');
    disp([10,'']);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Add one point to limit cycle data to close the plot  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Progress<=3
    LimitCycle_X=cell(samples,1);       % State-space values for the limit cycle
    LimitCycle_T=cell(samples,1);       % Time values for the limit cycle
    LimitCycle_Torques=cell(samples,1); % Torque values for the limit cycle
    
    for i=1:samples
        if LimitCycle_Period(i)==1
            LastLine=LimitCycle_Xb{i}(1,:);
            LLOrdered=[LastLine(1:2),LastLine(4),LastLine(3),LastLine(6),LastLine(5),LastLine(7)];
        else
            LLOrdered=LimitCycle_Xb{i}(1,:);
        end
        LimitCycle_X{i}=[LimitCycle_Xb{i};LLOrdered];
        LimitCycle_T{i}=[LimitCycle_Tb{i};LimitCycle_Tb{i}(end)];
        LimitCycle_Torques{i}=[LimitCycle_Torquesb{i}';LimitCycle_Torquesb{i}(:,end)'];
    end

    save(filename,'LimitCycle_X','LimitCycle_T','LimitCycle_Torques','-append');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate linearized Poincare Map Matrix and obtain its eigenvalues %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Progress<=3
    disp('Calculating the numerical Poincare map...');
    
    LC_eig1=zeros(5,samples);
    LC_eig2=zeros(5,samples);
    LC_eig3=zeros(5,samples);
    LC_eig4=zeros(5,samples);
    
    parfor i=1:samples
        switch LimitCycle_Period(i)
            case 1
                % Run a single simulation for the only
                % set of initial conditions
                [Result,LC_eig1(:,i)]=Poincare5(Slopes(i), LimitCycle_IC(1:5,i), 1, 2);
                if Result==1
                    disp(['Succesfully processed eigenvalues for ',num2str(Slopes(i)),char(186),' slope']);
                else
                    disp(['Failed to process eigenvalues for ',num2str(Slopes(i)),char(186),' slope']);
                end
                
            case 2
                % Run one simulation for each set of initial conditions
                [Result,LC_eig1(:,i)]=Poincare5(Slopes(i), LimitCycle_IC(1:5,i), 2, 2);
                if Result==1
                    disp(['Succesfully processed eigenvalues for ',num2str(Slopes(i)),char(186),' slope (1)']);
                else
                    disp(['Failed to process eigenvalues for ',num2str(Slopes(i)),char(186),' slope (1)']);
                end
                
                [Result,LC_eig2(:,i)]=Poincare5(Slopes(i), LimitCycle_IC(6:10,i), 2, 2);
                if Result==1
                    disp(['Succesfully processed eigenvalues for ',num2str(Slopes(i)),char(186),' slope (2)']);
                else
                    disp(['Failed to process eigenvalues for ',num2str(Slopes(i)),char(186),' slope (2)']);
                end
                
            case 4
                % Run one simulation for each set of initial conditions
                [Result,LC_eig1(:,i)]=Poincare5(Slopes(i), LimitCycle_IC(1:5,i), 4, 2);
                if Result==1
                    disp(['Succesfully processed eigenvalues for ',num2str(Slopes(i)),char(186),' slope (1)']);
                else
                    disp(['Failed to process eigenvalues for ',num2str(Slopes(i)),char(186),' slope (1)']);
                end
                
                [Result,LC_eig2(:,i)]=Poincare5(Slopes(i), LimitCycle_IC(6:10,i), 4, 2);
                if Result==1
                    disp(['Succesfully processed eigenvalues for ',num2str(Slopes(i)),char(186),' slope (2)']);
                else
                    disp(['Failed to process eigenvalues for ',num2str(Slopes(i)),char(186),' slope (2)']);
                end
                
                [Result,LC_eig3(:,i)]=Poincare5(Slopes(i), LimitCycle_IC(11:15,i), 4, 2);
                if Result==1
                    disp(['Succesfully processed eigenvalues for ',num2str(Slopes(i)),char(186),' slope (3)']);
                else
                    disp(['Failed to process eigenvalues for ',num2str(Slopes(i)),char(186),' slope (3)']);
                end
                
                [Result,LC_eig4(:,i)]=Poincare5(Slopes(i), LimitCycle_IC(16:20,i), 4, 2); %#ok<PFBNS>
                if Result==1
                    disp(['Succesfully processed eigenvalues for ',num2str(Slopes(i)),char(186),' slope (4)']);
                else
                    disp(['Failed to process eigenvalues for ',num2str(Slopes(i)),char(186),' slope (4)']);
                end
        end
    end
    
    LimitCycle_eigb=[LC_eig1;LC_eig2;LC_eig3;LC_eig4]; %#ok<*NASGU>
    
    Progress=4;
    save(filename,'LimitCycle_eigb','Progress','-append');
    disp([10,'']);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Sort out numerical data to avoid "jumping" lines %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Progress<=4
    disp('Sorting out numerical data to avoid "jumping" lines...');
    
    [LimitCycle_eig(1:5,:)] = ConnectVectors(LimitCycle_eigb(1:5,:));
    [LimitCycle_eig(6:10,:)] = ConnectVectors(LimitCycle_eigb(6:10,:));
    
    Progress=5;
    save(filename,'Progress','LimitCycle_eig','-append');
end

function [Lines] = ConnectVectors( Lines )
    N=size(Lines,2);
    
    for d=1:N-1
        L=length(find(Lines(:,d)~=0));
        if L==10
            % To shorten the computation time of calculating
            % all the permutations of a 10 elements vector,
            % we remove the 2 eigenvalues that are closest
            % to zero (there's a zero eigenvalue in the PM)
            
            CurVec=[Lines(1:L,d), (1:L)']; % Vector that holds the current points
            NextVec=[Lines(1:L,d+1), (1:L)']; % Vector that holds the next set of points
            
            CurVecS=sortrows(CurVec,1);
            NextVecS=sortrows(NextVec,1);
            
            % Grab the last 8
            CurVec=CurVecS(3:end,1);
            NextVec=NextVecS(3:end,1);
            
            L=8;
            if length(CurVec)~=length(NextVec)
                A=1;
            end
        else
            CurVec=Lines(1:L,d); % Vector that holds the current points
            NextVec=Lines(1:L,d+1); % Vector that holds the next set of points
        end

        % Calculate all possible permutations of the next points
        % This takes a LONG LONG LONG while when M is large
        PermVec=perms(NextVec);

        % Calculate distance from current vector to all permutations
        Dist=abs(PermVec-repmat(CurVec',size(PermVec,1),1));

        % Find the most distant point for each permutation
        MaxDist=max(Dist,[],2);

        % Find the permutation where the max. distance is the smallest
        MinPerm=find(MaxDist==min(MaxDist));

        MP=length(MinPerm);
        if MP==1
            NextVec=PermVec(MinPerm,:)';
            if L~=8
                Lines(1:L,d+1)=NextVec;
            else
                Lines(CurVecS(1:2,2),d+1)=Lines(NextVecS(1:2,2),d+1);
                Lines(CurVecS(3:end,2),d+1)=NextVec;
            end
        else            
            % Too many best permutations found
            % We'll select the one with the minimal overall distance
            BestPerms=PermVec(MinPerm,:);

            % Calculate absoulte distance from CurVec to each candidate
            Dist=abs(BestPerms-repmat(CurVec',MP,1));
            AbsDist=diag(Dist*Dist');

            MinPerm=find(AbsDist==min(AbsDist));
            if length(MinPerm)>1
                % screw it, choose the first one
                NextVec=BestPerms(MinPerm(1),:)';
            else
                NextVec=BestPerms(MinPerm,:)';
            end
            
            if L~=8
                Lines(1:L,d+1)=NextVec;
            else
                Lines(CurVecS(1:2,2),d+1)=Lines(NextVecS(1:2,2),d+1);
                Lines(CurVecS(3:end,2),d+1)=NextVec;
            end
        end
    end
end
        

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Find the convergence rate for each slope  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Progress<=5
    Epsilon=0.001;
    ConvRate=zeros(samples,5);
    NumSteps=20;

    disp('Calculating convergence rates...');
    parfor i=1:samples
        InitCond=LimitCycle_IC(1:5,i);
        InitCond1=InitCond+[Epsilon,0,0,0,0]';
        InitCond2=InitCond+[0,Epsilon,0,0,0]';
        InitCond3=InitCond+[0,0,Epsilon,0,0]';
        InitCond4=InitCond+[0,0,0,Epsilon,0]';
        InitCond5=InitCond+[0,0,0,0,Epsilon]';

        % Run simulations with different disturbances
        [IC_Store1,~,~]=Simulation(InitCond1,0,Slopes(i),Slopes(i),2,[2,NumSteps],0);
        [IC_Store2,~,~]=Simulation(InitCond2,0,Slopes(i),Slopes(i),2,[2,NumSteps],0);
        [IC_Store3,~,~]=Simulation(InitCond3,0,Slopes(i),Slopes(i),2,[2,NumSteps],0);
        [IC_Store4,~,~]=Simulation(InitCond4,0,Slopes(i),Slopes(i),2,[2,NumSteps],0);
        [IC_Store5,~,~]=Simulation(InitCond5,0,Slopes(i),Slopes(i),2,[2,NumSteps],0);

        % Calculate change in initial conditions for each step
        % with respect to steady state 
        ICdiff1=IC_Store1-repmat(InitCond,1,10);
        ICdiff2=IC_Store2-repmat(InitCond,1,10);
        ICdiff3=IC_Store3-repmat(InitCond,1,10);
        ICdiff4=IC_Store4-repmat(InitCond,1,10);
        ICdiff5=IC_Store5-repmat(InitCond,1,10);

        % Calculate change in the perturbed coordinate
        nDiff1=ICdiff1(1,:)/Epsilon;
        nDiff2=ICdiff2(2,:)/Epsilon;
        nDiff3=ICdiff3(3,:)/Epsilon;
        nDiff4=ICdiff4(4,:)/Epsilon;
        nDiff5=ICdiff5(5,:)/Epsilon;

        % Obtain rate of change
        Lambda1=abs(nDiff1.^(1./(NumSteps:-1:NumSteps-9)));
        Lambda2=abs(nDiff2.^(1./(NumSteps:-1:NumSteps-9)));
        Lambda3=abs(nDiff3.^(1./(NumSteps:-1:NumSteps-9)));
        Lambda4=abs(nDiff4.^(1./(NumSteps:-1:NumSteps-9)));
        Lambda5=abs(nDiff5.^(1./(NumSteps:-1:NumSteps-9)));

        ConvRate(i,:)=[mean(Lambda1);mean(Lambda2);mean(Lambda3);mean(Lambda4);mean(Lambda5)];
    end
    
    Progress=6;
    save(filename,'Progress','ConvRate','-append');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Calculate Potential, Kinetic, and other energies  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Robot = CompassBiped();
CPG = Controller();
CPG.Adaptive=2;
    
if Progress<=6
    disp('Calculating energy consumption...');
    StepLength=zeros(samples,1);
    dKineticEIm=zeros(samples,1);
    dKineticEFr=zeros(samples,1);
    dPotentialE=zeros(samples,1);
    ControllerE=zeros(samples,1);
    AbsControllerE=zeros(samples,1);

    % Arrays for double period data
    numDP=length(find(LimitCycle_Period==2));
    DPSlopes=zeros(numDP,1);
    DPStepLength=zeros(numDP,2);
    DPdKineticEIm=zeros(numDP,2);
    DPdKineticEFr=zeros(numDP,2);
    DPdPotentialE=zeros(numDP,2);
    DPControllerE=zeros(numDP,2);
    DPAbsControllerE=zeros(numDP,2);
    j=0;

    for i=1:samples
        if LimitCycle_Period(i)==1
            StepLength(i)=Robot.GetStepLength(LimitCycle_X{i}(1,3:end));
            dKineticEIm(i)=dKineticEnergyImpact(LimitCycle_X{i}(end-1,3:6));
            dKineticEFr(i)=dKineticEnergyFriction(LimitCycle_T{i}(1:end-1), LimitCycle_X{i}(1:end-1,3:6));
            dPotentialE(i)=dPotentialEnergy(LimitCycle_X{i}(1,3:end),Slopes(i));
            ControllerE(i)=ControllerEnergy(LimitCycle_T{i}, LimitCycle_X{i}(:,3:7), LimitCycle_Torques{i});
            AbsControllerE(i)=AbsControllerEnergy(LimitCycle_T{i}, LimitCycle_X{i}(:,3:7), LimitCycle_Torques{i});
        else
            j=j+1;
            DPSlopes(j)=Slopes(i);
            % Find the spot where impact occurs
            diffLC=diff(LimitCycle_X{i}(:,5));
            Idx=find(diffLC>mean(diffLC)+3*std(diffLC));
            DPStepLength(j,1)=Robot.GetStepLength(LimitCycle_X{i}(1,3:end));
            DPStepLength(j,2)=Robot.GetStepLength(LimitCycle_X{i}(Idx(1)+1,3:end));
            DPdKineticEIm(j,1)=dKineticEnergyImpact(LimitCycle_X{i}(Idx(1),3:6),LimitCycle_X{i}(1,3:6));
            DPdKineticEIm(j,2)=dKineticEnergyImpact(LimitCycle_X{i}(Idx(2),3:6),LimitCycle_X{i}(Idx(1)+1,3:6));
            DPdKineticEFr(j,1)=dKineticEnergyFriction(LimitCycle_T{i}(1:Idx(1)), LimitCycle_X{i}(1:Idx(1),3:6));
            DPdKineticEFr(j,2)=dKineticEnergyFriction(LimitCycle_T{i}(Idx(1)+1:end-1), LimitCycle_X{i}(Idx(1)+1:end-1,3:6));
            DPdPotentialE(j,1)=dPotentialEnergy(LimitCycle_X{i}(1,3:end),Slopes(i));
            DPdPotentialE(j,2)=dPotentialEnergy(LimitCycle_X{i}(Idx(1)+1,3:end),Slopes(i));
            DPControllerE(j,1)=ControllerEnergy(LimitCycle_T{i}(1:Idx(1)), LimitCycle_X{i}(1:Idx(1),3:7), LimitCycle_Torques{i}(1:Idx(1),:));
            DPControllerE(j,2)=ControllerEnergy(LimitCycle_T{i}(Idx(1)+1:end-1), LimitCycle_X{i}(Idx(1)+1:end-1,3:7), LimitCycle_Torques{i}(Idx(1)+1:end-1,:));
            DPAbsControllerE(j,1)=AbsControllerEnergy(LimitCycle_T{i}(1:Idx(1)), LimitCycle_X{i}(1:Idx(1),3:7), LimitCycle_Torques{i}(1:Idx(1),:));
            DPAbsControllerE(j,2)=AbsControllerEnergy(LimitCycle_T{i}(Idx(1)+1:end-1), LimitCycle_X{i}(Idx(1)+1:end-1,3:7), LimitCycle_Torques{i}(Idx(1)+1:end-1,:));

            StepLength(i)=DPStepLength(j,1)+DPStepLength(j,2);
            dKineticEIm(i)=DPdKineticEIm(j,1)+DPdKineticEIm(j,2);
            dKineticEFr(i)=DPdKineticEFr(j,1)+DPdKineticEFr(j,2);
            dPotentialE(i)=DPdPotentialE(j,1)+DPdPotentialE(j,2);
            ControllerE(i)=DPControllerE(j,1)+DPControllerE(j,2);
            AbsControllerE(i)=DPAbsControllerE(j,1)+DPAbsControllerE(j,2);
        end
    end
    
    Cmt=zeros(1,samples);
    Cet=zeros(1,samples);
    Cmt2=zeros(1,samples);
    Cet2=zeros(1,samples);
    Weight=Robot.GetWeight();
    for i=1:samples
        Cmt(i)=CostOfTransport(AbsControllerE(i),Weight,StepLength(i));
        Cet(i)=CostOfTransport(AbsControllerE(i),dPotentialE(i),Weight,StepLength(i));
        Cmt2(i)=CostOfTransport(ControllerE(i),Weight,StepLength(i));
        Cet2(i)=CostOfTransport(ControllerE(i),dPotentialE(i),Weight,StepLength(i));
    end
    
    Progress=7;
    save(filename,'Progress','StepLength',  'dKineticEIm',  'dKineticEFr',  'dPotentialE',  'ControllerE',  'AbsControllerE',...
                  'DPSlopes','DPStepLength','DPdKineticEIm','DPdKineticEFr','DPdPotentialE','DPControllerE','DPAbsControllerE',...
                  'Cmt','Cet','Cmt2','Cet2','-append');
end

    function [dKEi] = dKineticEnergyImpact(varargin)
        switch nargin
            case 1
                X=varargin{1};
                Period=1;
            case 2
                X0=varargin{1};
                X1=varargin{2};
                Period=2;
        end
        
        if Period==1
            KE0=Robot.GetKineticEnergy(X);
            Xa=Robot.CalcImpact(X);
            KE1=Robot.GetKineticEnergy(Xa);
            dKEi=KE1-KE0;        
        else
            if abs(X0(1)-X1(1)<1e-4)
                KE0=Robot.GetKineticEnergy(X0);
                KE1=Robot.GetKineticEnergy(X1);
                dKEi=KE1-KE0;
            else
                X0Ordered=[X0(2) X0(1) X0(4) X0(3)];
                X1Ordered=[X1(2) X1(1) X1(4) X1(3)];
                KE0=Robot.GetKineticEnergy(X0Ordered);
                KE1=Robot.GetKineticEnergy(X1Ordered);
                dKEi=KE1-KE0;
            end
        end
    end

    function [dKEf] = dKineticEnergyFriction(T,X)
        dKEf=trapz(T,-Robot.dampS*X(:,3).^2)+trapz(T,-Robot.dampNS*X(:,4).^2);
    end

    function [dPE] = dPotentialEnergy(X,Slope)
        Length=Robot.GetStepLength(X);
        deltay=Length*sind(Slope);
        dPE=(2*Robot.m+Robot.mh)*Robot.g*deltay;
    end

    function [CE] = ControllerEnergy(T,X,Control)     
        CE=trapz(T,Control(:,1).*X(:,4))+trapz(T,Control(:,2).*X(:,3));
    end

    function [CE] = AbsControllerEnergy(T,X,Control)     
        CE=trapz(T,abs(Control(:,1).*X(:,4)))+trapz(T,abs(Control(:,2).*X(:,3)));
    end

    function [COT] = CostOfTransport(varargin)
        if nargin==3
            ControlEffort=varargin{1};
            Weight=varargin{2};
            Distance=varargin{3};
            
            COT=ControlEffort/(Weight*Distance);
        end
        if nargin==4
            ControlEffort=varargin{1};
            Potential=varargin{2};
            Weight=varargin{3};
            Distance=varargin{4};
            
            COT=(ControlEffort-Potential)/(Weight*Distance);
        end
        
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Calculate Ground Reaction Forces  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Progress<8
    disp('Calculating ground reaction forces...');
    GRF=cell(1,samples); % Actual GRF (Fx,Fy)
    MeasuredGRF=cell(1,samples); % GRF projection on stance leg
    MeanMeasuredGRF=zeros(1,samples); % Mean of this ^
    MaxMeasuredGRF=zeros(1,samples); % Max of this ^^
    MaxGRFRatio=zeros(1,samples); % Max ratio Fx/Fy (for friction assumption)
    
    parfor s=1:samples
        N=length(LimitCycle_X{s});
        
        GRF{s}=zeros(2,N);
        MeasuredGRF{s}=zeros(1,N);
        MaxRatio=0;
        for i=1:N
            % Get ground reaction force (Fx, Fy) for a whole step
            % for each slope        
            GRF{s}(:,i)=Robot.GetGRF(LimitCycle_X{s}(i,3:end));
            
            % Calculate angle between GRF and leg
            alpha=atan2(GRF{s}(2,i),GRF{s}(1,i))-LimitCycle_X{s}(i,3)-pi/2;
            % Project GRF on stance leg
            MeasuredGRF{s}(i)=sqrt(GRF{s}(2,i)^2+GRF{s}(1,i)^2)*cos(alpha);
            
            % Find max ratio, i.e. friction coefficient
            % required to avoid slipage
            Ratio=abs(GRF{s}(1,i)/GRF{s}(2,i));
            if Ratio>MaxRatio
                MaxRatio=Ratio;
            end
        end
        
        MaxGRFRatio(s)=MaxRatio;
        MeanMeasuredGRF(s)=mean(MeasuredGRF{s});
        MaxMeasuredGRF(s)=max(MeasuredGRF{s});
    end
    
    Progress=8;
    save(filename,'Progress','GRF','MeasuredGRF','MeanMeasuredGRF','MaxMeasuredGRF','MaxGRFRatio','-append');
end
    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Find Max slope reached VS Curvature  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Progress<9
    disp('Calculating relation of max slope VS Terrain curvature...');
    samples2=10;
    MaxUpSlope=zeros(1,samples2); % Max slope reached uphill
    MaxDownSlope=zeros(1,samples2); % Max slope reached downhill
    Curvatures=linspace(0.005,0.1,samples2); % Ground curvature
    
    parfor s=1:samples2
        % Run uphill simulation
        [EndReached, EndText, T, X, TorquesP] = Simulation( InitCond0, 2, 0, Curvatures(s), 2, 1, 0 ); %#ok<ASGLU>
        
        % Obtain max slope from last stance position
        Floor = Terrain(2,0,Curvatures(s));
        MaxUpSlope(s)=Floor.SurfSlope(X(end,1));
        
        % Run downhill simulation
        [EndReached, EndText, T, X, TorquesP] = Simulation( InitCond0, 2, 0, -Curvatures(s), 2, 1, 0 ); %#ok<ASGLU>
        
        % Obtain max slope from last stance position
        Floor = Terrain(2,0,-Curvatures(s));
        MaxDownSlope(s)=Floor.SurfSlope(X(end,1));
    end
    
    Progress=9;
    save(filename,'Progress','MaxUpSlope','MaxDownSlope','Curvatures','-append');
end
    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% ZMP Calculations  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Progress<10
    disp('Performing GRF and ZMP calculations...');
    
    LimitCycle_GRF = cell(samples,1);
    LimitCycle_ZMP = cell(samples,1);
    MaxFriction = zeros(samples,1);
    MaxZMP = zeros(samples,1);
    
    parfor s=1:samples
        ThisX = LimitCycle_X{s};
        ThisT = LimitCycle_T{s};
        ThisTorque = LimitCycle_Torques{s};
        
        NT = length(ThisT);
        ThisGRF = zeros(NT,2);
        ThisZMP = zeros(NT,1);
        
        % Save stance foot position
        SFoot = ThisX(1,1);
        for t=1:NT
            % Check for a step taking place
            if ThisX(t,1) ~= SFoot
                % Switch legs
                ThisTemp = ThisX(t:end,[3,5]);
                ThisX(t:end,[3,5]) = ThisX(t:end,[4,6]);
                ThisX(t:end,[4,6]) = ThisTemp;
                
                ThisTemp = ThisTorque(t:end,1);
                ThisTorque(t:end,1) = ThisTorque(t:end,2);
                ThisTorque(t:end,2) = ThisTemp;
                
                SFoot = ThisX(t,1);
            end
            
            ThisGRF(t,:) = Robot.GetGRF(ThisX(t,3:end))';
            ThisZMP(t) = ThisTorque(t,2)/ThisGRF(t,2);
        end
        
        LimitCycle_GRF{s} = ThisGRF;
        LimitCycle_ZMP{s} = ThisZMP;
        
        % Find max values
        MaxFriction(s) = max(abs(ThisGRF(:,1)./ThisGRF(:,2)));
        Min = min(ThisZMP);
        Max = max(ThisZMP);
        if abs(Min)>abs(Max)
            MaxZMP(s) = Min;
        else
            MaxZMP(s) = Max;
        end
    end
    
    Progress = 10;
    save(filename,'Progress','LimitCycle_GRF','LimitCycle_ZMP','MaxFriction','MaxZMP','-append');
end
    

disp('All done!!!');
end

