
close all; clear all; clc

% cases = {'GA only','GA + NN','GA + rescale','GA + NN + rescale'};
% fileIn = {[],[],[],[]};
% 
% for i=1:4
%     GA_try_Matsuoka(cases{1,i},fileIn{1,i});
% end

cases = {'GA + NN + rescale'};
fileIn = {[]};

for i=1:2
    GA_try_Matsuoka(cases{1,i},fileIn{1,i});
end