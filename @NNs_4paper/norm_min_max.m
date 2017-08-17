function vec_norm = norm_min_max(obj,vec,paramName)
% this function normalize the CPG parametrs based on the min max values
%            p     - p_min
%   p_norm = -------------
%            p_max - p_min
%
% the CPG parameters range is taken from the genomeRanges in 'MML' class
%
% Inputs:
% 'vec' - row vector with the parameters values
% 'pramName' - can be every name from 'seqOrder'

paramID = strcmp(paramName,obj.seqOrder);

samplNum = size(vec,1);

minRange = ones(samplNum,1) * obj.MML.Gen.Range(1,paramID);
maxRange = ones(samplNum,1) * obj.MML.Gen.Range(2,paramID);

vec_norm = (vec - minRange) ./ (maxRange - minRange);

end

