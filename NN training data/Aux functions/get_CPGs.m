function [results_filtered,period_filtered,seq_filtered,ids] =...
    get_CPGs(results,Type,MML)
%GET_CPGS - take 'results' structure and filter out the right 'Type' of
%CPGs. either Oscillatory ones or not oscillatory ones.
% 'MML' - is the Matsuoka Sim class


% get CPG periods:
periods = horzcat(results(:).periods);

if size(periods,1) == 2 % 4N CPG:
    % Filter CPG's where not both signals oscillating:
    osc_ids_temp = ~isnan(periods);
    osc_ids_temp = osc_ids_temp(1,:) & osc_ids_temp(2,:);
    disp(['Number of non-osc CPGs: ',num2str(sum(~osc_ids_temp))]);

    % Filter CPG's where the is a big difference between hip and ankle:
    periods_ratios = (periods(1,:)./periods(2,:));
    diff_ids = (periods_ratios >  0.85) & (periods_ratios <  1.15); 
    disp(['Number of CPGs with not matching periods: (from the osc ones)',...
        num2str(sum(osc_ids_temp & ~diff_ids))]);
    periods = mean(periods,1);

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

    
    
elseif size(periods,1) == 1  % 2N CPG:
    % Filter CPG's where not both signals oscillating:
    osc_ids_temp = ~isnan(periods);
    disp(['Number of non-osc CPGs: ',num2str(sum(~osc_ids_temp))]);

    diff_ids = ones(size(osc_ids_temp));
    
else
    error('invalid period structure');
    
end

% % check that all of the parameters are in the genome range:
seq = (vertcat(results(:).seq))';
ids_in_genome_range = true(1,size(seq,2));
for n=1:MML.Gen.Length
    ids_temp = (seq(n,:) > MML.Gen.Range(1,n)) &...
        (seq(n,:) < MML.Gen.Range(2,n));
    if ~any(ids_temp)
        disp('found parameter not in genome range:');
        disp(['the param is num: ',num2str(n),...
            '    ',num2str(sum(~ids_temp)),...
            ' param ecceed range']);
    end
    ids_in_genome_range = ids_in_genome_range & ids_temp;
end
disp(['Number of CPGs with parameters not in range: ',...
    num2str(sum(~ids_in_genome_range))]);

disp(['num of osc CPG before range checking: ',...
    num2str(sum(osc_ids_temp & diff_ids))]);

% get the oscillatory ones:
osc_ids = osc_ids_temp & diff_ids & ids_in_genome_range;
disp(['num of oscillating CPGs: ',...
    num2str(sum(osc_ids))]);

% check which ones oscillates in range:
osc_inRange_ids = osc_ids &...
    ( (periods(1,:) > MML.perLimOut(1,1)) &...
    (periods(1,:) < MML.perLimOut(1,2)) );

disp(['num of CPG which oscillate in range: ',...
    num2str(sum(osc_inRange_ids))]);

switch Type
    case {'osc','oscillating'}
        ids = osc_ids;
        results_filtered = results(ids);
        period_filtered = periods(ids);
        seq_filtered = seq(:,ids);
    case {'osc_in_per_range'}
        ids = osc_inRange_ids;
        results_filtered = results(ids);
        period_filtered = periods(ids);
        seq_filtered = seq(:,ids);
    case {'n-osc','not oscillatory'}
        % we want oscillatory CPGs with genome in range:
        ids = ~osc_ids_temp & ids_in_genome_range;
        results_filtered = results(ids);
        period_filtered = periods(ids);
        seq_filtered = seq(:,ids);
    otherwise
        error('invalid type: can be either "osc" or "n-osc"');
end
end

