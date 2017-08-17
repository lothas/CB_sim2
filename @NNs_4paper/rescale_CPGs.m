function obj = rescale_CPGs(obj)
% this function acts exactly as the rescaling function in MOOGA

seq = vertcat(obj.results(:).seq);


ids_nan = isnan(obj.periods);
ids_inRange = (obj.periods > obj.MML.perLim(1)) &...
    (obj.periods < obj.MML.perLim(2));

ids_to_change = ~ids_nan &  ~ids_inRange;

% Select new random period within desired range
des_period = obj.MML.perLim(1) +...
    rand(size(obj.periods,2),1)*(obj.MML.perLim(2)-obj.MML.perLim(1));

des_period_2change = (des_period(ids_to_change))';
period_2change = obj.periods(ids_to_change);

% Scale Tr, Ta to obtain desired period
ratio = des_period_2change./period_2change;

t_rescaled = seq(:,1);
t_rescaled(ids_to_change,1) = seq(ids_to_change,1).*...
    ratio';

ids_tau_out_of_range = (t_rescaled(:,1) < obj.MML.Gen.Range(1,1)) |...
        (t_rescaled(:,1) > obj.MML.Gen.Range(2,1));

num_of_outOfRange = sum(ids_tau_out_of_range);
maxRange = ones(1,num_of_outOfRange)*obj.MML.Gen.Range(1,1);
minRange = ones(1,num_of_outOfRange)*obj.MML.Gen.Range(2,1);
t_rescaled(ids_tau_out_of_range,1) = min(max(t_rescaled(ids_tau_out_of_range,1),maxRange'),...
    minRange');

obj.tau_rescaled = t_rescaled;

end

