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
        legends = {'MOGA','MOGA + NN',...
            'MOGA + re-scaling','MOGA + NN + re-scaling'};
        % titles additions for different cases
        titleAdd = {'GA only','GA + NN','GA + rescale','everything'};
        
        % the order of the fitnesses from the GA run:
        fitnessOrder = {'VelFit','NrgEffFit',...
            'VelRangeFit #1','VelRangeFit #2','VelRangeFit #3',...
            'VelRangeFit #4','VelRangeFit #5','VelRangeFit #6',...
            'VelRangeFit #7','VelRangeFit #8','EigenFit'};
        
        % the parametrs seq order in the Matsuoka CPG:
        seqOrder = {'tau','b','c_1','c_2','c_3','c_4',...
            'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
            'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}'};
        
        % the parametrs seq order in the Matsuoka CPG: (including the
        % adaption parametrs)
        seqOrder_extend = {'tau','b','c_1','c_2','c_3','c_4',...
                         'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
                         'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}',...
                         'ks_\tau','ks_c1','ks_c2','ks_c3','ks_c4'};

        % Matsuoka CPG simulation class (contains the paramets ranges)
        MML = [];
        
        
    end
    
    methods
        function obj = plotMOOGA4Paper(InFiles_names)
            obj.data_names = InFiles_names;
            for i=1:4
                obj.data{1,i} = load(InFiles_names{1,i});
            end
            
            % define the class for CPG simulation:
            MML = MatsuokaML();
            MML.perLim = [0.68 0.78];
            MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
            MML.tStep = 0.05;
            MML.tEnd = 15;
            MML.nNeurons = 4;
            
            obj.MML = MML;
        end
        
    end
    
end

