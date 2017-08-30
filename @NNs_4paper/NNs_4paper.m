classdef NNs_4paper
    % this class contains all of the NN analisys for paper
    
    properties
        
        % name of the file contain the training data
        results_fileName = [];
        MML = []; %Matsuoka object for the sims (contains the parametrs border);
        
        % the order of the parametrs in CPG Sequence:
        seqOrder = {'tau','b','c_1','c_2','c_3','c_4',...
            'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
            'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}'};

        results = [];
        periods = []; % Save the periods
        
        osc_ids = []; % ids of oscillating CPGs
        num_of_osc_ids = [] % the number of oscillating CPGs
        num_of_osc_ids_exluded_param_range = []; % self checking about 
                                                 % the amount of osc ids when i dont
                                                 % care about the parameters range
        
        osc_inRange_ids = []; % ids CPGs which oscillates in range
        num_of_inRange_ids = [] % the number of oscillating CPGs
        
        Inputs_names_train = [];
        Inputs_train = []; % inputs to NN or MoE
        
        % same but with 'period_desired' instead of 'periods'
        Inputs_names_actual = [];
        Inputs_actual = [];
        
        Targets_names = [];
        Targets = []; % targets to NN or MoE
        
        % Neural network object
        NN = [];
        NN_training_data = [];
        
        % "Mixture of Experts" object:
        MoE = [];
        
        % save cahnges in 'seq':
        seq = [];
        tau_rescaled = []; %change in seq due to rescaling of 'tau'
        seq_NN = [];
        
        disp_information = true; % show garphs and stuff while training
        
    end
    
    methods
        function obj = NNs_4paper(results_fileName,MatsuokaSimObject)
            obj.results_fileName = results_fileName;
            
            res = load(results_fileName,'results');
            obj.results = res.results;
            
%             % Show header:
%             hed = load(results_fileName,'header');
%             disp(' The header of the data file:');
%             disp(hed.header);

            obj.MML = MatsuokaSimObject;
            
            obj.seq = vertcat(obj.results(:).seq);
            
            obj = obj.filter_ids_4N_general;
%             obj = obj.filter_ids_2N_Symm;
            
            
        end
        
        function obj = filter_ids_4N_general(obj)
            % this function get the CPGs in 'results' and filer out the osc CPGs and
            % the CPGs with missdefined periods.

            % get CPG periods:
            periods = horzcat(obj.results(:).periods);
            
            % Filter CPG's where not both signals oscillating:
            osc_ids_temp = ~isnan(periods);
            osc_ids_temp = osc_ids_temp(1,:) & osc_ids_temp(2,:);
            disp(['Number of non-osc CPGs: ',num2str(sum(~osc_ids_temp))]);
            
            % Filter CPG's where the is a big difference between hip and ankle:
            periods_ratios = (periods(1,:)./periods(2,:));
            diff_ids = (periods_ratios >  0.85) & (periods_ratios <  1.15); 
            disp(['Number of CPGs with not matching periods: (from the osc ones)',...
                num2str(sum(osc_ids_temp & ~diff_ids))]);
            obj.periods = mean(periods,1);

            % % plot the distribution of the missdefined CPG periods:
            if false
                figure;
                h=histogram(periods_ratios,100); grid minor;
                h.BinLimits = [0,2.5];
                h.BinWidth = 0.1;
                h.Normalization = 'pdf';
                xlabel('P_{hip} / P_{ankle} : the ratio between the two CPG outputs');
                ylabel('probability density');
                title('histogram of the ratio between the two CPG outputs');
                set(gca,'FontSize',10);
                savefig('figure_TBD_Histogram_of_ratio_between_periods_hipAnkle')
            end
            
            % % check that all of the parameters are in the genome range:
            seq = (vertcat(obj.results(:).seq))';
            ids_in_genome_range = true(1,size(seq,2));
            for n=1:obj.MML.Gen.Length
                ids_temp = (seq(n,:) > obj.MML.Gen.Range(1,n)) &...
                    (seq(n,:) < obj.MML.Gen.Range(2,n));
                ids_in_genome_range = ids_in_genome_range & ids_temp;
            end
            disp(['Number of CPGs with parameters not in range: ',...
                num2str(sum(~ids_in_genome_range))]);
            
            obj.num_of_osc_ids_exluded_param_range = sum(osc_ids_temp & diff_ids);
            obj.osc_ids = osc_ids_temp & diff_ids & ids_in_genome_range;

            obj.osc_inRange_ids = obj.osc_ids &...
                ( (periods(1,:) > obj.MML.perLimOut(1,1)) &...
                (periods(1,:) < obj.MML.perLimOut(1,2)) );
            
            obj.num_of_osc_ids = sum(obj.osc_ids);
            obj.num_of_inRange_ids = sum(obj.osc_inRange_ids);
            
            % save the good periods
%             periods_mean = mean(periods,1);
%             obj.periods = periods_mean(1,obj.osc_ids);
        end
        
        function obj = filter_ids_2N_Symm(obj)
            % same as 'filter_ids_4N_general()' but for 2N symmetric CPG
            
            % Change the sequence order to match the new case:
            obj.seqOrder = ...
                {'tau','b','c','NR','a','NR','NR','NR'}; % 'NR'-not relevant
            
            % get CPG periods:
            periods = horzcat(obj.results(:).periods);
            
            % Filter CPG's where not both signals oscillating:
            osc_ids_temp = ~isnan(periods);
            disp(['Number of non-osc CPGs: ',num2str(sum(~osc_ids_temp))]);
            
            % Filter CPG's where the is a big difference between hip and ankle:
            % if we have only 2 Matsuoka neurons,
            %   then assume that the period is correct
            diff_ids = true(1,size(periods,2));
            obj.periods = periods;
            
            % % check that all of the parameters are in the genome range:
            seq = (vertcat(obj.results(:).seq))';
            ids_in_genome_range = true(1,size(seq,2));
            
            for n=1:obj.MML.Gen.Length
                if ~strcmp(obj.seqOrder{1,n},'NR')
                    ids_temp = (seq(n,:) > obj.MML.Gen.Range(1,n)) &...
                        (seq(n,:) < obj.MML.Gen.Range(2,n));
                    ids_in_genome_range = ids_in_genome_range & ids_temp;
                end
            end
            
            disp(['Number of CPGs with parameters not in range: ',...
                num2str(sum(~ids_in_genome_range))]);
            
            obj.num_of_osc_ids_exluded_param_range = sum(osc_ids_temp & diff_ids);
            obj.osc_ids = osc_ids_temp & diff_ids & ids_in_genome_range;

            obj.osc_inRange_ids = obj.osc_ids &...
                ( (periods(1,:) > obj.MML.perLimOut(1,1)) &...
                (periods(1,:) < obj.MML.perLimOut(1,2)) );
            
            obj.num_of_osc_ids = sum(obj.osc_ids);
            obj.num_of_inRange_ids = sum(obj.osc_inRange_ids);
            
        end
    end
    
end

