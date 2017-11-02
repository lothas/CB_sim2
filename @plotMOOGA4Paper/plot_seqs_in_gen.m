function plot_seqs_in_gen(obj,whichGA,generation,IDS)
%PLOT_SEQS_IN_GEN plot the whanted seqs for a given generation
% 
% INPUTS:
% *) 'whichGA' - the number of the wanted GA file from class
% *) 'generation' - the number of the desired generation
% *) 'IDS' - index row vector containing the seqs to show

genes = squeeze(obj.data{1,whichGA}.GA.Seqs(:,:,generation));
wanted_genes = genes(IDS,:);
wanted_genes = wanted_genes';

fits = squeeze(obj.data{1,whichGA}.GA.Fit(:,:,generation));
wanted_fits = fits(IDS,:);

SPnum = ceil(sqrt(size(IDS,2))); % number of subplots 

figure
for i=1:size(IDS,2)
    subplot(SPnum,SPnum,i);
    Title = ['gene #',num2str(IDS(i)),'   Fits: ',num2str(wanted_fits(i,1:3)),' '];
    obj.plot_seq(wanted_genes(:,i),Title);
end

end
