classdef plotMOOGA4Paper
    % This class is use to plot the MOOGA results for our paper
    
    properties
        
        % 1x4 cell array contain the names of GA data files: 
        data_names = [];
        data = [];
        
        % make the cluster colors constant for clarity:
        colors = {'r.','b.','g.','m.','c.','k.'};
        colors1 = {'r','b','g','m','c','k'};
        
        % legend for graph:
        Legends = [];
        
        % the order of the fitnesses from the GA run:
        fitnessOrder = {'VelFit','NrgEffFit',...
            'VelRangeFit #1','VelRangeFit #2','VelRangeFit #3',...
            'VelRangeFit #4','VelRangeFit #5','VelRangeFit #6',...
            'VelRangeFit #7','VelRangeFit #8','EigenFit'};
        
        % the parametrs seq order in the Matsuoka CPG:
        seqOrder = [];

        % Matsuoka CPG simulation class (contains the paramets ranges)
        MML = [];
        
        
    end
    
    methods
        function obj = ...
                plotMOOGA4Paper(MML_in,InFiles_names,Legends1,seqorder)
            
            % get MatsuokaML class:
            obj.MML = MML_in;
            
            % get the GA data from the MOGa results files:
            obj.data_names = InFiles_names;
            for i=1:numel(InFiles_names)
                obj.data{1,i} = load(InFiles_names{1,i});
            end
            
            % get the names of the Type of results in the MOGA files:
            obj.Legends = Legends1;
            
            % get the order of the parameters in the results.seq
            obj.seqOrder = seqorder;
        end
        
    end
    
end

