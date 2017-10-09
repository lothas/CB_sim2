function [ax] = twoVec_hist_compare(vec1,vec2,vec_name,...
    binNum,legendNames,typeOfGraph)

% This function get two vector (don't need to be in the same length)
%    and compare the histograms

% Inputs:
% *) 'vec1' and 'vec2' - two vectors
% *) 'vec_name' - the name of the parametrs in the two vectors
% *) 'binNum' - the number of bins in the histogram
% *) 'typeOfGraph' - 'bars' - histogram in bars form
%                    'plot' - hist in line  form
% *) 'legendNames' - the name of each vector                   

% check the hypothesis that both vector came from the same distribution:
[pvalue,rejection] = ranksum(vec1,vec2);
str1 = sprintf('p_value of %s is %0.3f and the hypotesis got %d',...
    vec_name,pvalue,rejection);
disp([str1,' real pval is ',num2str(pvalue)]);

% check "kullbacl-Leiblar" divergence (dist of the two distrubutions)
prepare_KL_div_calc(vec1,vec2)

switch typeOfGraph
    case 'bars'
        % % % Plot traditional Histograms:
        % figure;
        histogram(vec1,binNum,'Normalization','pdf'); grid minor; hold on;
        histogram(vec2,binNum,'Normalization','pdf');
        xlabel(vec_name);
        ylabel('probability density');
        title(['histogram of ',vec_name]);
        legend(legendNames);
        set(gca,'FontSize',13);
        hold off;
    case 'plot'
        % Plot Histograms with a smooth line instead of bars:
        [Est_pdf_vec1, Edges_pdf_vec1] = ...
            histcounts(vec1,binNum, 'Normalization','pdf');
        bin_center_vec1 = (Edges_pdf_vec1(1:end-1)+ Edges_pdf_vec1(2:end))/2;

        [Est_pdf_vec2, Edges_pdf_vec2] = ...
            histcounts(vec2,binNum, 'Normalization','pdf');
        bin_center_vec2 = (Edges_pdf_vec2(1:end-1)+ Edges_pdf_vec2(2:end))/2;
        % figure;
        plot(bin_center_vec1,Est_pdf_vec1,'LineWidth',2); hold on;
        plot(bin_center_vec2,Est_pdf_vec2,'LineWidth',2); grid minor;
        xlabel(vec_name);
        ylabel('probability density');
        title(['histogram of ',vec_name]);
        legend(legendNames);
        axis tight
        set(gca,'FontSize',13);
        fileName = sprintf('figure_TBD_Histogram_of_%s',vec_name);
        savefig(fileName)
        hold off;
end

ax = gca;

end

