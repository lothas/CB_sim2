close all; clear all;
cur_date = now;
File1 = ['GA_',datestr(cur_date,'mm_dd_hh_MM'),'_St1.mat'];
File2 = ['GA_',datestr(cur_date,'mm_dd_hh_MM'),'_St2.mat'];

% Run MOOGA with large population
GA_ROBIO( 4, 6000, 'GA_11_06_08_33_St1.mat', File1 );
% Optimize with smaller population
GA_ROBIO( 10, 1000, File1, File2 );