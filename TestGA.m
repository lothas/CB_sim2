function [  ] = TestGA(  )
GA = MOOGA();
NPop = 2500;
TopPop = NPop/5;
NGen = 20;

Pop = GA.Gen.RandSeq(NPop);
plot3(Pop(:,1),Pop(:,2),Pop(:,3),'o')
    axis([0 1 0 1 0 1])
view(35,35)
grid
title('Initial population');

for g = 1:NGen
    pause(0.1)
    NewPop = zeros(size(Pop));
    TopIDs = GA.GetTopPop(Pop, TopPop); % fitness = genes
    
    % Transfer top IDs to new population
    NewPop(1:TopPop,:) = Pop(TopIDs,:);
    
    % Add a mutated copy of top IDs
    NewPop(TopPop+1:2*TopPop,:) = GA.Gen.Mutate(NewPop(1:TopPop,:));
    
    % Add children (by pairs)
    for c = 1:NPop/2-TopPop
        % Select 2 parents randomly
        IDs = randsample(1:TopPop,2);
        NewPop(2*TopPop+2*c-1:2*TopPop+2*c,:) = ...
            GA.Gen.Crossover(NewPop(IDs(1),:),NewPop(IDs(2),:),'1point');
    end
    
    % Show new population
    Pop = NewPop;
    cla
    plot3(Pop(:,1),Pop(:,2),Pop(:,3),'o')
    axis([0 1 0 1 0 1])
    title(['Generation ',num2str(g)]);    
end

