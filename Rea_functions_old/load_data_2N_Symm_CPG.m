function [results_temp,periods_temp]=load_data_2N_Symm_CPG(howMuchData)
% this function load the data to work space
% for the 2neuron CPG (symmetric)

% inputs- "howMuchData": how many datasamples we want to load (before
% checking which ones hace periods!


fileNames = {'MatsRandomRes_2Neurons_08_02_2017.mat';};

for i=1:ceil((howMuchData/10000))
    if i>size(fileNames,1)
        disp('not enough data');
        break;
    end
    load(fileNames{i,1},'results');
    periods = horzcat(results(:).periods);
    if i == 1
        results_temp = results; 
        periods_temp = periods;
    else
        results_temp = horzcat(results_temp,results);
        periods_temp = horzcat(periods_temp,periods);
    end
    clear results periods
end

end