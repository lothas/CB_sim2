function GA_Stage1(redo)
    % Version 2.3 - 23/03/14
    
    % Version 2.2 - 14/03/14
    
    % Version 2.1 - 02/06/13
    
    % New version:
    % Once a robot gets a good height grade (i.e. it managed to take a few
    % steps) two new simulations are run to check the slope fitness
    
    % Version 1.0 - 13/12/12
    % This version uses 3 objectives to look for good parameters that will
    % cause the robot to walk in a stable, fast, energy efficient manner.
    % The best solutions are selected using pareto fronts and paired
    % randomly.
    
    outer_tic = tic;
    
    if nargin<1
        redo=0;
    end
    
    if matlabpool('size')==0
        matlabpool open 7 % Work in parallel to finish faster
    end
    
    tend_start=10;      % Duration of simulations for 1st gen (sec)
    tend_increase=1;   % extend duration each generation
    tend_max=60;       % Max duration
    tend=tend_start;    
        
   
    % Use manual solution
%     Genes(:,1,1)=[1.106512566, 0.4579, 0.9, 0.65, -13.6255, 0, 0.1,   95,   80,...
%                                                     4.6281, 0, 0.4, -443, -295,...
%                                                    13.6255, 0, 0.1,   95,   80,...
%                                                          0, 0,   0,    0,    0, 0];
                                 
    % Use genomes from previous run
%     Data=load('GAResults.mat');
%     Genes(:,:,1)=Data.Genes(:,:,end);
    
    Filename='GAResults.mat';
              
    g0=0;     % Last generation where there was no change in fitness
    gcounter=0; % Counter since last gen when there was a change in fitness
    LastMaxFit=zeros(NumObj,1);
    
    % Check that genomes loaded OK before starting
    if ~isempty(find(Genes(1,:,1) == 0, 1))
        disp('ERROR: Some genomes are missing');
        error=1;
    end
    
    for g=InitGen:Generations
        inner_tic = tic;
        disp(['Running Generation ',num2str(g),'     -     ',datestr(now,'HH:MM:SS')]);
        % Run simulation for each individual
        parfor p=StartG:Population
            Genome=Genes(:,p,g);
            [EndReached, EndText, T, X, Torques]=Simulation(Genome,NumActJoints,tend,0); %#ok<ASGLU>
            close all
            % disp(EndText);
            ThisFit=zeros(1,NumObj);
            ThisFit(1:3)=Fitness(T,X,Torques,tend);
            
            % Make sure the controller walks forward
            if ThisFit(2)<0
                % Change the sign of all torques
                for i=1:NumTorques
                    Genome(3*i)=-Genome(3*i);
                end
                
                % Update results
                Genes(:,p,g)=Genome;
                ThisFit(2)=-ThisFit(2);
            end
           
            Fit(:,p,g)=ThisFit;
        end
        
        StartG=TopPop+1;
        
        if g<Generations
            % Select top TopPop fit controllers
            TopIDs = GetTopIDs(Fit(:,:,g),TopPop);
%             disp(sort(TopIDs)');
            
            % Copy top genomes to next generation
            Genes(:,1:TopPop,g+1)=Genes(:,TopIDs,g);
            Fit(:,1:TopPop,g+1)=Fit(:,TopIDs,g);
           
            % Add a mutated copy of top genomes
            for i=1:TopPop
                MutGene = Genes(:,i,g+1);
                Genes(:,TopPop+i,g+1)=Mutate(MutGene,GenesMin,GenesMax,GenProb(SingMutProb,g,g0,gcounter));
            end

            % Pair top controllers randomly
            Pairs=zeros(NumPairs,2);
            NewPop=2*TopPop+1;
            TopIDs_copy=TopIDs;
            for p=1:NumPairs
                pick=ceil(rand(1)*length(TopIDs_copy));
                Pairs(p,1)=TopIDs_copy(pick);

                % Remove pick from top IDs
                TopIDs_copy(pick)=[];

                pick=ceil(rand(1)*length(TopIDs_copy));
                Pairs(p,2)=TopIDs_copy(pick);

                % Remove pick from top IDs
                TopIDs_copy(pick)=[];

                % Crossover and mutate genes
                for of=1:NumOffs
                    Offspring=CrossOver(Genes(:,Pairs(p,1),g),Genes(:,Pairs(p,2),g));
                    
                    if rand(1)<=MutProb
                        Genes(:,NewPop,g+1)=Mutate(Offspring,GenesMin,GenesMax,GenProb(SingMutProb,g,g0,gcounter));
                        % the Mutate function checks the genome after mutation
                    else
                        Genes(:,NewPop,g+1)=Offspring;
                    end 
                    
                    NewPop=NewPop+1;
                end
                
                if NewPop>Population
                    % Stop reproducing if population is full
                    break
                end
            end

            % Add 1 random new individual
            Genes(:,end,g+1)=RandomGenome(GenesMin,GenesMax);
            
            % Increase simulation time
            tend=min(tend_start+g*tend_increase,tend_max);
        end
        
        % Display Generation fitness
        AvgFit=zeros(NumObj,1);
        TopFit=zeros(NumObj,1);
        for i=1:NumObj
            AvgFit(i)=mean(Fit(i,:,g));
            TopFit(i)=max(Fit(i,:,g));
        end
        
        disp(['Top Gen. Fitness: H=',num2str(TopFit(1)),...
            ' D=',num2str(TopFit(2)),' E=',num2str(TopFit(3))]);
        disp(['Average Gen. Fitness: H=',num2str(AvgFit(1)),...
            ' D=',num2str(AvgFit(2)),' E=',num2str(AvgFit(3))]);
        
        if sum(TopFit>LastMaxFit)==0
            % There was no change in fitness
            g0 = g;
            gcount=gcount+1;
            disp(['There hasn''t been an improvement in ',num2str(gcount),' generations']);
        else
            disp(['Change in max fitness: H=',num2str(TopFit(1)-LastMaxFit(1)),...
                ' D=',num2str(TopFit(2)-LastMaxFit(2)),' E=',num2str(TopFit(3)-LastMaxFit(3))]);
            LastMaxFit=TopFit;
            gcount=0;
        end
        
        t_diff = toc(inner_tic);
        minutes = floor(t_diff/60);
        seconds = mod(t_diff,60);
        if minutes>0
            disp(['Time to run this generation: ',num2str(minutes),''' ',num2str(seconds,'%.2f'),'"',10])
        else
            disp(['Time to run this generation: ',num2str(seconds,'%.2f'),'"',10])
        end
        
        InitGen=g; %#ok<NASGU>
            
        save(Filename);
    end
    
    toc(outer_tic)
    
    disp('Preparing Data Analysis');
    DataAnalysis(1);

function [Prob] = GenProb(Base,cur_gen,last_gen,num_no_change)
    Prob = 0.01 + Base * (1-(cur_gen-last_gen)/Generations) + min(0.08*min(num_no_change,10),0.89-Base);
end

function [fit] = Fitness(T,X,Torques,tend)
    fit=zeros(3,1);
    
    Robot=CompassBiped();
    Weight=2*Robot.m+Robot.mh;
    
    steps=length(T);
    HipX=zeros(1,steps);
    HipY=zeros(1,steps);
    StanceX=zeros(1,steps);
    StanceY=zeros(1,steps);
    SwingX=zeros(1,steps);
    SwingY=zeros(1,steps);
    for t=1:steps
        [HipX(t),HipY(t)]=Robot.GetPos(X(t,3:end),'Hip');
        [StanceX(t),StanceY(t)]=Robot.GetPos(X(t,3:end),'S');
        [SwingX(t),SwingY(t)]=Robot.GetPos(X(t,3:end),'NS');
    end
    % Make X position absolute
    StanceX=StanceX+X(:,1)';
    SwingX=SwingX+X(:,1)';
    PosX = min(StanceX,SwingX);
    
    % Give points for keeping the hip above a certain height
    fit(1)=HeightFitness(T,HipY,Robot.L,tend);
    
    % Give points for distance traveled
    fit(2)=DistanceFitness(T,PosX,Robot.L,tend);
    
    % Give points for energy efficiency
    [HipX0,HipY0]=Robot.GetPos(X(1,3:end),'Hip');
    [HipX1,HipY1]=Robot.GetPos(X(end,3:end),'Hip');
    % make hip position absolute
    HipX0=HipX0+X(1,1);
    HipX1=HipX1+X(end,1);
    HipY0=HipY0+X(1,2);
    HipY1=HipY1+X(end,2);
    
    if fit(1)>0.9
        fit(3)=EnergyFitness(T,X,Torques,[HipX0,HipY0],[HipX1,HipY1],Weight,0);
    else
        fit(3)=EnergyFitness(T,X,Torques,[HipX0,HipY0],[HipX1,HipY1],Weight,3*Robot.L);
    end
end

function [Points] = EnergyFitness(T,X,Torques,Hip0,Hip1,Weight,MinDist)
    % Calculate absolute control effort
    ControlEffort=trapz(T,abs((Torques(1,:)+Torques(2,:)).*X(:,5)'))+trapz(T,abs(Torques(2,:).*X(:,6)'));
%     ControlEffort=trapz(T,abs(Torques(1,:).*X(:,5)'))+trapz(T,abs(Torques(2,:).*(X(:,6)'-X(:,5)')));
    
    % Calculate difference in potential energy
    dPotentialE=Weight*9.81*(Hip1(2)-Hip0(2));
    
    % Calculate distance travelled
    DistanceTravelled=abs(Hip1(1)-Hip0(1));
    
    if DistanceTravelled>MinDist
        % Calculate Cost Of Transport
        COT=(ControlEffort-dPotentialE)/(Weight*DistanceTravelled);

        % Low COT (as low as 0) is best so points are given by
        Points=1/(1+5*COT);
        % COT of 0 gives 1
        % COT of 0.03 gives 0.869
        % COT of 0.12 gives 0.625
        % COT of 0.3 gives 0.4
    else
        Points=0;
    end
end
