function results = load_results(results_fileName)
%LOAD_RESULTS take the name of the file/files which stores the results
%   and outputs the result structure
% 
% Inputs:
% *) 'results_fileName' - cell array containing the names of the files

load(results_fileName{1,1},'results','header');
disp('data file information:');
disp(header);
% if I have more than 1 data file:
for i=2:numel(results_fileName)
    data = load(results_fileName{1,i},'results');
    results = [results, data.results]; %#ok<AGROW>
end

end

