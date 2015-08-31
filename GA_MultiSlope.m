close all; clear all;
gen = 12;
pop = 800;
file_in = 'GA_-150_08_29_18_41.mat';
for i = [15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5]
    GA = GA_Slope( gen, pop, -i, file_in );
    file_in = GA.FileOut;
end