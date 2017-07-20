function dist = KL_div_4paper(vec1,vec2)
% this function get two vectors (don't need to be the same length).
% than, it uses hist counts with predefined bins in order to create two
% distributions with the same x-axis (defined as bins center)
% and then it uses the 'kldiv.m' function to get the "Kullback-leibler divergence

vec1_min = min(vec1);
vec1_max = max(vec1);
vec2_min = min(vec2);
vec2_max = max(vec2);

tot_min = min(vec1_min,vec2_min);
tot_max = min(vec1_max,vec2_max);

binsNum = 20;

edges = linspace(tot_min,tot_max,binsNum);

Prob_vec1 = histcounts(vec1,edges, 'Normalization', 'probability');
Prob_vec2 = histcounts(vec2,edges, 'Normalization', 'probability');

bin_center = (edges(1:end-1)+ edges(2:end))/2;

dist = kldiv(bin_center,Prob_vec1,Prob_vec2);

figure;
plot(bin_center,Prob_vec1,'LineWidth',2); hold on;
plot(bin_center,Prob_vec2,'LineWidth',2);
legend('n-osc CPGs','osc CPGs');
hold off;

end

