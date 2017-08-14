function [normParam] = param_norm_minMax(MML,data,gen_num)
% this function takes the CPG param and normalize them according to:
%            p     - p_min
%   p_norm = -------------
%            p_max - p_min
%
% the CPG parameters range is taken from the genomeRanges in 'MML' class
%
% Inputs:
% 'MML' - Matsuoka sim class
% 'data' - contain the data from the GA results
% 'gen_num' - the generation number


paramRanges = MML.Gen.Range;
% paramNames = MML.Gen.Keys;
allSeqs = data.GA.Seqs;
Param = zeros(size(allSeqs,1),size(allSeqs,2));
normParam = zeros(size(allSeqs,1),size(allSeqs,2));

for i=1:23
    Param(:,i) = allSeqs(:,i,gen_num);
    minParam = paramRanges(1,i);
    maxParam = paramRanges(2,i);
    normParam(:,i) = ( Param(:,i) - minParam ) ./ ( maxParam - minParam );
end


end

