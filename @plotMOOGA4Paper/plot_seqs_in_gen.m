function plot_seqs_in_gen(obj,whichGA,generation,num2show)
%PLOT_SEQS_IN_GEN plot the whanted seqs for a given generation
% 
% INPUTS:
% *) 'whichGA' - the number of the wanted GA file from class
% *) 'generation' - the number of the desired generation
% *) 'IDS' - index row vector containing the seqs to show

genes = squeeze(obj.data{1,whichGA}.GA.Seqs(:,:,generation));
wanted_genes = genes';

fits = squeeze(obj.data{1,whichGA}.GA.Fit(:,:,generation));
wanted_fits = fits';

IDs = zeros(1,num2show);
count = 1;
for i=1:num2show
    IDs(1,i) = find(wanted_fits(count,:) == max(wanted_fits(count,:)),1);
    wanted_fits(count,IDs(1,i)) = 0;
    if count == 3 
        count = 1;
    else
        count = count + 1;
    end
end

SPnum = ceil(sqrt(num2show)); % number of subplots 

chosen_genes = genes(IDs,:)';
chosen_fits = fits(IDs,:)';

figure
for i=1:num2show
    subplot(SPnum,SPnum,i);
    Title = ['gene #',num2str(IDs(i)),'   Fits: ',num2str(chosen_fits(1:3,i)'),' '];
    obj.plot_seq(chosen_genes(:,i),Title);
end

end
