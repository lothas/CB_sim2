function [obj] = filter_ids(obj)
% this function get the CPGs in 'results' and filer out the osc CPGs and
% the CPGs with missdefined periods.

% get CPG periods:
periods = horzcat(obj.results(:).periods);

% Filter CPG's where not both signals oscillating:
osc_ids = ~isnan(periods);
osc_ids = osc_ids(1,:) & osc_ids(2,:);

% Filter CPG's where the is a big difference between hip and ankle:
periods_ratios = (periods(1,:)./periods(2,:));
diff_ids = (periods_ratios >  0.85) & (periods_ratios <  1.15); 

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

ids = osc_ids & diff_ids;

obj.ids_des_period = ids &...
    ( (periods(1,:) > obj.MML.perLimOut(1,1)) &...
    (periods(1,:) < obj.MML.perLimOut(1,2)) );

obj.ids_not_des_period = ids & ...
    ( (periods(1,:) < obj.MML.perLimOut(1,1)) |...
    (periods(1,:) > obj.MML.perLimOut(1,2)) );

obj.ids_osc = ids;
obj.ids_des_period;

end

