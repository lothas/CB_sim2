function [ inputs,targ ] = prepareData_2Neurons( results,periods,ids,NNinput,NNtarget )
% this function takes the results and periods and prepare data data for NN
% for Matsuoka 2Neurons symmetric CPG
% this function replaces the "prepareData" function.

% inputs: 
% 1) results - the results structure as saved from jonathan's code.
% 2) periods - the period vector as saved from jonathan's code.
% 3) ids - ids of the genes that we want to work with.
% 4) NNinput - the NN inputs variables
% 5) NNtarget - the NN target.
% 6) 'What_or_W_Flag' - logic, if we want the Matsuoka weights to be
%                       normalized or not

% outputs:
% inputs - the input to the NN
% targ - the NN targets
%


MatsuokaParam = vertcat(results(ids).seq)'; % extracting the parameters from the structure to a matrix
MatsuokaParam = MatsuokaParam(1:5,:); % get rid of the the unimportant 'gen' (the k's for the control).

tau = MatsuokaParam(1,:);
b = MatsuokaParam(2,:);
c_i = MatsuokaParam(3,:); % ignore  MatsuokaParam(4,:) (not relevant!)
W_ij = MatsuokaParam(5,:);

Period = periods(ids);
freq = 1./periods(ids);

%% making the NNinputs:
for i=1:length(NNinput)
   switch NNinput{1,i};
       case 'tau'
           inputs(i,:) = tau;
       case 'b'
           inputs(i,:) = b;
       case 's'
           inputs(i,:) = c_i;
       case 'a'
           inputs(i,:) = W_ij;
       case 'period'
           inputs(i,:) = Period;
       case 'freq'
           inputs(i,:) = freq;
       otherwise
           warning('no such string');
   end
end

%% making the NNtarget:
for i=1:length(NNtarget)
   switch NNtarget{1,i};
       case 'tau'
           targ(i,:) = tau;
       case 'b'
           targ(i,:) = b;
       case 's'
           targ(i,:) = c_i;
       case 'a'
           targ(i,:) = W_ij;
       case 'period'
           targ(i,:) = Period;
       case 'freq'
           targ(i,:) = freq;
       otherwise
           warning('no such string');
   end
end

end

