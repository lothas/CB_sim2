
function save_data_CSV(obj,caseNum,fileName)

    % save inputs and targets matrices to CSV file to be read by 
    % TENSORFLOW or Eureqa.
    
    % the file Header is printed to a seperate text file with the same name
    
   
    %% get inputs and targets names
    [Inputs_names,Targets_names] = ...
        obj.check_NN_case(caseNum,'period');
   
    %% make NN inputs & outputs matrices:
    [inputs,targets] = ...
        obj.prepare_NN_train_data(Inputs_names,Targets_names);
    
    %% make them coloum matrices (samples in rows)
    inputs = inputs';
    targets = targets';
    %% make files names:
    dataFileName = ['fileDate_',fileName,datestr(now,'mm_dd_hh_MM'),'.csv'];
    headerFileName = ['fileDate_',fileName,datestr(now,'mm_dd_hh_MM'),'.txt'];
    %% save data file as CVS
    csvwrite(dataFileName,[inputs,targets]);
    

    %% make header file:
    source_of_data = ['data file is: ',obj.results_fileName,'.mat'];
    % NOTE: make sure the you change this order according to the relevant
    %       sequence!
    
    % get inputs order:
    inputsOrder = [sprintf('The inputs order is: \n ')];
    for i=1:length(Inputs_names)
        inputsOrder = [inputsOrder,...
            sprintf(' %d ) " %s " \n ',i,Inputs_names{1,i})];
    end

    % get targets order:
    targetsOrder = [sprintf('The targets order is: \n ')];
    for i=1:length(Targets_names)
        targetsOrder = [targetsOrder,...
            sprintf(' %d ) " %s " \n ',i,Targets_names{1,i})];
    end
    
    fid = fopen(headerFileName,'wt');
    fprintf(fid,'\n');
    fprintf(fid,[source_of_data,'\n']);
    fprintf(fid,[inputsOrder,'\n']);
    fprintf(fid,[targetsOrder,'\n']);
    fprintf(fid,['inputs first, targets second','\n']);
    fclose(fid);


end
