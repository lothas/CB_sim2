function plot_best_seq_over_gen(obj,whichGA,generations)
%PLOT_BEST_SEQ_OVER_GEN plot the first seq for each generation
% 
% INPUTS:
% *) 'whichGA' - the number of the wanted GA file from class
% *) 'generations' - row vector with the desired generations to show

genes = obj.data{1,whichGA}.GA.Seqs(:,:,generations);

best_genes = squeeze(genes(1,:,:));

SPnum = ceil(sqrt(size(generations,2))); % number of subplots 

figure
for i=1:size(generations,2)
    subplot(SPnum,SPnum,i)
    Title = ['gen #',generations(i)];
    obj.plot_seq(best_genes(:,i),Title)
end

end

