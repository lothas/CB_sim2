function [normParam] = param_norm_minMax(obj,gen_num,whichCase)
% this function takes the CPG param and normalize them according to:
%            p     - p_min
%   p_norm = -------------
%            p_max - p_min
%
% the CPG parameters range is taken from the genomeRanges in 'MML' class
%
% Inputs:
% 'gen_num' - the generation number

% TODO: we don't need 'MML'. the 'GA' already contain the raanges!
%       change it.

paramRanges = obj.MML.Gen.Range;
% paramNames = MML.Gen.Keys;
allSeqs = obj.data{1,whichCase}.GA.Seqs;
Param = zeros(size(allSeqs,1),size(allSeqs,2));
normParam = zeros(size(allSeqs,1),size(allSeqs,2));

for i=1:23
    Param(:,i) = allSeqs(:,i,gen_num);
    minParam = paramRanges(1,i);
    maxParam = paramRanges(2,i);
    normParam(:,i) = ( Param(:,i) - minParam ) ./ ( maxParam - minParam );
end


end

