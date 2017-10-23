function generate_GenomeFile(whichCase)
%GENERATE_GENOMEFILE simple and easy function to generates the genome files.
% all in one place. so different sripts and function wouldn't overwrite the
% genomeFile.
% 
% Input:
% *) 'whichCase' - can be:
%       #) '2N_symm' - for 2neuron symmetric CPG
%       #) '2N_general' - for 2neuron general CPG
%       #) '4N_tagaLike' - for 4neuron taga-like CPG
%       #) '4N_general' - for 4neuron general CPG
% 
% NOTE: change the "keys" if you need to change the GA optimization
% parameters.

% define Mutation strength:
MutDelta0 = 0.04;   MutDelta1 = 0.02;

nAnkle = 1;%1; % Number of ankle torques
nHip = 1;   % Number of hip torques
N = nAnkle+nHip;

switch whichCase
    case '2N_symm'
        genome_file = 'MatsuokaGenome_2Neuron_Symm.mat';
        maxAnkle = 0;   % Max ankle torque
        maxHip = 10;    % Max hip torque
        Mamp = [10,10];
        mamp = 0*Mamp;
        
        Mw = [10,10];
        mw = 0*Mw;
        
        Keys = {'\tau_r', 'beta','amp_2n_same_inputs',    '2neuron_symm_weights', 'ks_\tau',     'ks_c_2n_symm', 'IC_matsuoka';
                      1 ,      1,                   2,                         1,        1 ,          1,            0 };
        Range = {  0.02 ,    0.2,                  mw,                         0,      -10 ,         -1; % Min
                   0.25 ,    2.5,                  Mw,                        10,       10 ,          1}; % Max

        % Note: because of some old mistakes. the tonic input gene ('c') 
        %   is encoded with two values. only one of them is used.          
    case '2N_general'
        genome_file = 'MatsuokaGenome_2Neuron_General.mat';
        maxAnkle = 0;   % Max ankle torque
        maxHip = 10;    % Max hip torque
        Mamp = [10,10];
        mamp = 0*Mamp;
        
        Mw = [10,10];
%         Mw = [5,5];
        mw = 0*Mw;
        
        Keys = {'\tau_r', 'beta', 'amp_2n_dif_inputs',    '2neuron_general_weights', 'ks_\tau', 'ks_c_2n_general', 'IC_matsuoka';
                      1 ,      1,                   2,                            2,        1 ,                 2,            0 };
        Range = {  0.02 ,    0.2,                mamp,                           mw,      -10 ,           [-1,-1]; % Min
                   0.25 ,    2.5,                Mamp,                           Mw,       10 ,             [1,1]}; % Max


    case '4N_tagaLike'
        genome_file = 'MatsuokaGenome_4Neuron_tagaLike.mat';
        maxAnkle = 20;   % Max ankle torque
        maxHip = 8;    % Max hip torque
        Mamp = [maxAnkle*ones(1,2*nAnkle), maxHip*ones(1,2*nHip)];
        mamp = 0*Mamp;
        
        mw = 0*ones(1,4);
        Mw = 10*ones(1,4);
        
        % CPG strucute: (ALSO Symm W_ij = W_ji)
        %   H_F   H_E           % 
        % 4 O-----O 3           %   
        %    \    |             %   w = [0  , W12, 0  , 0  ; 
        %     \   |             %        W21, 0  , w23, W24;
        %      \  |             %        0  , 0  , 0  , W34;
        %       \ |             %        0  , 0  , w43, 0  ;
        % 1 O-----O 2           % w12=w21 = w1  
        %  A_F    A_E           % w23 = w2
        %                       % w24 = w3
        %                       % w43=w34 = w4
        
        % Final genome with tau_r + beta (constant tau_u/tau_v ratio) 
        Keys = {'\tau_r', 'beta', 'amp_4n_symm',   '4neuron_taga_like', 'ks_\tau',     'ks_c_4n_symm', 'IC_matsuoka';
                      1 ,      1,             1,                     4,        1 ,                 1 ,            0 };
        Range = {  0.02 ,    0.2,             0,                    mw,      -10 ,                 -1; % Min
                   0.25 ,    2.5,        maxHip,                    Mw,       10 ,                  1}; % Max

        % Note: the current encoding is for symmetric tonic inputs. 
        %   you can change it. but don't forget to change 
        %   the adaptation coef ('k_c') as well       
        
%         % Final genome with tau_r + beta (constant tau_u/tau_v ratio) 
%         Keys = {'\tau_r', 'beta',      'amp',   '4neuron_taga_like', 'ks_\tau',             'ks_c', 'IC_matsuoka';
%                       1 ,      1,          4,                     4,        1 ,                  4,            0 };
%         Range = {  0.02 ,    0.2,       mamp,                    mw,      -10 ,          -0.1*Mamp; % Min
%                    0.25 ,    2.5,       Mamp,                    Mw,       10 ,           0.1*Mamp}; % Max


    case '4N_general'
        genome_file = 'MatsuokaGenome_4Neuron_general.mat';
        maxAnkle = 20;   % Max ankle torque
        maxHip = 8;    % Max hip torque
        Mamp = [maxAnkle*ones(1,2*nAnkle), maxHip*ones(1,2*nHip)];
        mamp = 0*Mamp;
        Mw = 10*ones(1,12);
        mw = 0*Mw;
        
          % Final genome with tau_r + beta (constant tau_u/tau_v ratio) 
        Keys = {'\tau_r', 'beta', 'amp',   'weights', 'ks_\tau',     'ks_c', 'IC_matsuoka';
                      1 ,      1,    4 ,          12,        1 ,         4 ,            0 };
        Range = {  0.02 ,    0.2,  mamp,          mw,      -10 ,  -0.1*Mamp; % Min
                   0.25 ,    2.5,  Mamp,          Mw,       10 ,   0.1*Mamp}; % Max

    otherwise
        error('invalid input');
end

disp('The defined Genome is:')
disp('Tau_retio = 12');
disp(['case: ',whichCase]);
for i=1:(length(Keys)-1)
    disp([Keys{1,i},' :']);
    disp(['    number of elements: ',num2str(Keys{2,i})]);
    disp(['    min: [',num2str(Range{1,i}),']']);
    disp(['    Max: [',num2str(Range{2,i}),']']);
end

save(genome_file, 'nAnkle', 'nHip', 'maxAnkle', 'maxHip', ...
    'Mamp', 'mamp', 'N', 'Mw', 'mw', ...
    'MutDelta0', 'MutDelta1', 'Keys', 'Range');

end

