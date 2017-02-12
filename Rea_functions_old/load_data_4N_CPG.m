function [results_temp,periods_temp]=load_data_4N_CPG(howMuchData)
% this function load the data to work space
% for the 4neuron CPG (general case)

% inputs- "howMuchData": how many datasamples we want to load (before
% checking which ones hace periods!

fileNames = {'MatsRandomRes_16_12_2016.mat';
    'MatsRandomRes_18_12_2016.mat';
    'MatsRandomRes_19_12_2016.mat';
    'MatsRandomRes_20_12_2016.mat';
    'MatsRandomRes_21_12_2016.mat';
    'MatsRandomRes_25_12_2016.mat';
    'MatsRandomRes_27_12_2016.mat';
    'MatsRandomRes_01_01_2016.mat';
    'MatsRandomRes_02_01_2016.mat';
    'MatsRandomRes_30_01_2016.mat';
    'MatsRandomRes_31_01_2016.mat';
    'MatsRandomRes_01_02_2016.mat';
    'MatsRandomRes_02_02_2016.mat';
    'MatsRandomRes_03_02_2016.mat';
    'MatsRandomRes_04_02_2016.mat';};

for i=1:ceil((howMuchData/200000))
    load(fileNames{i,1},'results','periods');
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

