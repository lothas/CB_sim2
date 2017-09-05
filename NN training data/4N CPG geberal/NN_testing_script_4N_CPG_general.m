
%% NN testing Script:

clear all; close all; clc;

% the order of the parametrs in CPG Sequence:
seqOrder = {'tau','b','c_1','c_2','c_3','c_4',...
    'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
    'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}'};

% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;
MML.nNeurons = 4;

% % change tau_a/tau_r to 12 (instead of 5)
MML.Sim.Con.tau_ratio = 12;

% file name for uploading:L
results_fileName = {'MatsRandomRes_4Neurons_Large_b_Large_W_All_osc'};

%% Load data:
load(results_fileName{1,1},'results','header');
disp('data file information:');
disp(header);
% if I have more than 1 data file:
for i=2:numel(results_fileName)
    data = load(results_fileName{1,i},'results');
    results = [results, data.results]; %#ok<AGROW>
end

clear data i
%% get and filter periods:

% get CPG periods:
periods = horzcat(results(:).periods);

% Filter CPG's where not both signals oscillating:
osc_ids_temp = ~isnan(periods);
osc_ids_temp = osc_ids_temp(1,:) & osc_ids_temp(2,:);
disp(['Number of non-osc CPGs: ',num2str(sum(~osc_ids_temp))]);

% Filter CPG's where the is a big difference between hip and ankle:
periods_ratios = (periods(1,:)./periods(2,:));
diff_ids = (periods_ratios >  0.85) & (periods_ratios <  1.15); 
disp(['Number of CPGs with not matching periods: (from the osc ones)',...
    num2str(sum(osc_ids_temp & ~diff_ids))]);
periods = mean(periods,1);

% % plot the distribution of the missdefined CPG periods:
if false
    figure;
    h=histogram(periods_ratios,100); grid minor;
    h.BinLimits = [0,2.5];
    h.BinWidth = 0.1;
    h.Normalization = 'pdf';
    xlabel('P_{hip} / P_{ankle} : the ratio between the two CPG outputs');
    ylabel('probability density');
    title('histogram of the ratio between the two CPG outputs');
    set(gca,'FontSize',10);
    savefig('figure_TBD_Histogram_of_ratio_between_periods_hipAnkle')
end

% % check that all of the parameters are in the genome range:
seq = (vertcat(results(:).seq))';
ids_in_genome_range = true(1,size(seq,2));
for n=1:MML.Gen.Length
    ids_temp = (seq(n,:) > MML.Gen.Range(1,n)) &...
        (seq(n,:) < MML.Gen.Range(2,n));
    ids_in_genome_range = ids_in_genome_range & ids_temp;
end
disp(['Number of CPGs with parameters not in range: ',...
    num2str(sum(~ids_in_genome_range))]);

num_of_osc_ids_exluded_param_range = sum(osc_ids_temp & diff_ids);
osc_ids = osc_ids_temp & diff_ids & ids_in_genome_range;

osc_inRange_ids = osc_ids &...
    ( (periods(1,:) > MML.perLimOut(1,1)) &...
    (periods(1,:) < MML.perLimOut(1,2)) );

num_of_osc_ids = sum(osc_ids);
num_of_inRange_ids = sum(osc_inRange_ids);

%% keep the good seq and periods
seq = seq(:,osc_ids);
periods = periods(:,osc_ids);

%% Prepare NN inputs and outputs:
input_names = {'periods','tau',...
    'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
    'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}'};

output_names = {'b'};

[sampl,targ] = ...
    prepare_NN_data(input_names,output_names,...
    seqOrder,seq,periods);

%% Neural Network:
architecture = [20];

net = fitnet(architecture);
net.trainFcn = 'trainbr';
net.trainParam.showWindow = 1; 

[net, tr] = train(net, sampl, targ);