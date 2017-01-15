
clc; clear all; close all;

% Load varibles: for 4 neuron data
load('MatsRandomRes_16_12_2016.mat','results','periods')
results1 = results; periods1 = periods;
clear results periods
load('MatsRandomRes_18_12_2016.mat','results','periods')
results2 = results; periods2 = periods;
clear results periods
load('MatsRandomRes_19_12_2016.mat','results','periods')
results3 = results; periods3 = periods;
clear results periods
load('MatsRandomRes_20_12_2016.mat','results','periods')
results4 = results; periods4 = periods;
clear results periods
load('MatsRandomRes_21_12_2016.mat','results','periods')
results5 = results; periods5 = periods;
clear results periods
load('MatsRandomRes_25_12_2016.mat','results','periods')
results6 = results; periods6 = periods;
clear results periods

% concetrate the data:
results = horzcat(results1,results2,results3,results4,results5,results6);
periods = horzcat(periods1,periods2,periods3,periods4,periods5,periods6);

clear results1 results2 results3 results4 results5 results6
clear periods1 periods2 periods3 periods4 periods5 periods6

clc; clear all; close all;