classdef class4paper
    % thi class contain all of the function needed for the data analisys of
    % the paper
    
    properties
        
        MML = []; %Matsuoka object for the sims
        seqOrder = {'tau','b','c_1','c_2','c_3','c_4',...
                    'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
                    'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}'}; % define the order in which the seq is coded
            % 'b','tau' - beta and tau in the 4 neuron Matsuoka CPG (by defualt
            %               T=5*tau)
            %        'w_12',...,'w_43' - weights
            %       'c_1',...'c_4' - tonic inputs
  
        results = []; % all off the Matsuoka sim results
        
        ids_osc = []; % only CPGs that converged
        ids_des_period = []; % only CPGs that converged in range
        ids_not_des_period = []; % CPGs that didn't convereged (not osc CPGs)
        % TODO: think what is better. to store the results? or only the ids?
        
        
    end
    
    methods
        % %%% % class constructor % %%% %
        function obj = class4paper(varagin)
            % define the class for CPG simulation:
            obj.MML = MatsuokaML();
            obj.MML.perLim = [0.68 0.78];
            obj.MML.perLimOut = obj.MML.perLim + [-0.08 0.08]; % Desired period range
            obj.MML.tStep = 0.05;
            obj.MML.tEnd = 15;
            obj.MML.nNeurons = 4;

            switch nargin
                case 1
                    obj.results = varagin;
                otherwise
                    error('too many inputs...');
            end
            
            % get osc CPGs:
            obj = obj.filter_ids();
        end
        
    end
    
end

