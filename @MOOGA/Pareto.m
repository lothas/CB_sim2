function [Fronts] = Pareto(GA, Data)
%PARETO Finds Pareto fronts for multi-objective genetic algorithms
%   The algorithm builds each front by "dropping" elements that are
%   dominated until only dominant or weak dominant elements remain.
%   These elements are removed from the set and the process is repeated to
%   find a second pareto front and so forth until all elements are
%   accounted for

% Version 0.4 - 11/05/2014

if length(Data)<1
    nI = 5; nJ = 100;
    % Run sample code
    % Build sample data
    x = zeros(nI*nJ,3);
    for i = 1:nI
        r = 10*i;
        for j = 1:nJ
            phi = pi*rand();
            th = pi*rand();
            x(nJ*(i-1)+j,:) = [...
                r*cos(phi)*cos(th), ...
                r*cos(phi)*sin(th), ...
                r*sin(phi)];
        end
    end
    x = abs(x);
    
    % Plot sample data before running algorithm
    figure()
    plot3(x(:,1),x(:,2),x(:,3),'+');
    
    % Save data to origx
    origx = x;
end

% Round off to 3 decimal places
% x = round(x*1000)/1000;

% Start separating layers
% Get number of samples and objectives
data = size(x,1);
% Give unique ID to each sample
x = [x (1:data)'];

NFronts = 1;
Fronts = {};
out = [];
% Sort by first objective
SortedData = sortrows(x,-1);

while 1
    % Drop samples that are "dominated"
    % i.e. that x1>x2 for each objective
    i = 1;
    Nd = size(SortedData,1);
    if Nd == 0
        break;
    end
    
    while i < Nd
        j = i+1;
        while j <= Nd
            if SortedData(j,1:end-1) <= SortedData(i,1:end-1)
                % Sample j is dominated by sample i
                SortedData(j,:)=[];
                Nd = Nd-1;
            else
                % Check next sample
                j = j+1;
            end
        end
        i = i+1;
    end
    
    % Save "undominated" samples to current front
    Fronts{NFronts} = SortedData(:,end); %#ok<AGROW>
    out = [out; Fronts{NFronts}]; %#ok<AGROW>
    NFronts = NFronts+1;
    
    % Restore samples to original data - samples already on fronts
    SortedData = x; SortedData(out,:) = [];
    SortedData = sortrows(SortedData,-1);
end

if length(Data)<1
    % Run sample code
    close all
    N=length(Fronts);

    figure()
    hold on
    for i=1:N
        if length(Fronts{i})<3
            continue;
        end
        x = origx(Fronts{i},1);
        y = origx(Fronts{i},2);
        z = origx(Fronts{i},3);
        tri = delaunay(x,y);
        trisurf(tri, x, y, z,'FaceColor',[1-i/N; i/2/N; i/N]);
    end
end

end

%% %%%%%%%%%%%%%% Slow method %%%%%%%%%%%%%% %%
% % Get number of samples and objectives
% data = size(x,1);
% % Give unique ID to each sample
% x=[x (1:data)'];
% tempx = x;
% NFronts=1;
% Fronts={};
% out = [];
% while 1
%     i = 1;
%     while i <= data
%         j = 1;
%         while j <= data
%             if i == j
%                 j = j+1;
%                 continue;
%             end
% 
%             if tempx(j,1:end-1) <= tempx(i,1:end-1)
%                 tempx(j,:) = [];
%                 data = data - 1;
%                 if i>j
%                     i = i-1;
%                 end
%             else
%                 j = j+1;
%             end
%         end
%         i = i+1;
%     end
% 
%     Fronts{NFronts} = tempx(:,end); %#ok<AGROW>
%     out = [out; Fronts{NFronts}]; %#ok<AGROW>
%     NFronts = NFronts+1;
% 
%     tempx = x; tempx(out,:) = [];
%     data = size(tempx,1);
%     if data == 0
%         break;
%     end
% end

