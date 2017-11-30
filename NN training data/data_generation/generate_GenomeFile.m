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

nAnkle = 1; % Number of ankle torques
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
    case '4N_tagaLike_1symmInput'
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
        Range = {  0.02 ,    0.2,             0,                    mw,      -10 ,                 -0.1*maxAnkle; % Min
                   0.25 ,    2.5,      maxAnkle,                    Mw,       10 ,                  0.1*maxAnkle}; % Max

    case '4N_tagaLike_generalInput'
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
        Keys = {'\tau_r', 'beta',      'amp_same4each_joint',   '4neuron_taga_like', 'ks_\tau','ks_c_same4each_joint', 'IC_matsuoka';
                      1 ,      1,                          2,                     4,        1 ,                     2,            0 };
        Range = {  0.02 ,    0.2,        0*[maxAnkle,maxHip],                    mw,      -10 , -0.1*[maxAnkle,maxHip]; % Min
                   0.25 ,    2.5,          [maxAnkle,maxHip],                    Mw,       10 ,  0.1*[maxAnkle,maxHip]}; % Max


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
               
    case '6N_tagaLike_2Ank_torques'
        nAnkle = 2; % Number of ankle torques
        nHip = 1;   % Number of hip torques
        N = nAnkle+nHip;

        genome_file = 'MatsuokaGenome_4Neuron_tagaLike.mat';
        maxAnkle = 20;   % Max ankle torque
        maxHip = 8;    % Max hip torque
        Mamp = [maxAnkle*ones(1,2*nAnkle), maxHip*ones(1,2*nHip)];
        mamp = 0*Mamp;
        
        mw = 0*ones(1,6);
        Mw = 5*ones(1,6);
        
        % CPG strucute: (ALSO Symm W_ij = W_ji)
        % 1) every pair of Extensor and reflexor neurons are connected to
        %   each other with symmetric weights.
        % 2) both of the ankles' Extensor neurons are connected to both of the hip neurons
        %
        %  A2_F    A2_E
        %(3) O-----O (4)
        %       /  |
        %      /   |
        %     /    |
        %    /     |
        %   H_F   H_E
        %(5) O-----O (6)  
        %     \    |             %   w = [0  , W12, 0  , 0  , 0  , 0; 
        %      \   |             %        W21, 0  , 0  , 0  , W25, W26;
        %       \  |             %        0  , 0  , 0  , W34, 0  , 0;
        %        \ |             %        0  , 0  , W43, 0  , W45, W46;
        %(1) O-----O (2)         %        0  , 0  , 0  , 0  , 0  , W56;
        %   A1_F    A1_E         %        0  , 0  , 0  , 0  , W65, 0;
        %                        
        %                        w12 = w21 = w34 = w43 = w1  
        %                        w56 = 65  = w2
        %                        w25 = w3
        %                        w26 = w4
        %                        w45 = w5
        %                        w46 = w6
        
        % Final genome with tau_r + beta (constant tau_u/tau_v ratio) 
        Keys = {'\tau_r', 'beta', 'amp_6n_symm',   '6neuron_taga_like', 'ks_\tau',     'ks_c_6n_symm', 'IC_matsuoka';
                      1 ,      1,             1,                     6,        1 ,                 1 ,            0 };
        Range = {  0.02 ,    0.2,             0,                    mw,      -10 ,                 -0.1*maxAnkle; % Min
                   0.25 ,    2.5,      maxAnkle,                    Mw,       10 ,                  0.1*maxAnkle}; % Max
case '6N_tagaLike_2Ank_torques_symm'
        nAnkle = 2; % Number of ankle torques
        nHip = 1;   % Number of hip torques
        N = nAnkle+nHip;

        genome_file = 'MatsuokaGenome_4Neuron_tagaLike.mat';
        maxAnkle = 20;   % Max ankle torque
        maxHip = 8;    % Max hip torque
        Mamp = [maxAnkle*ones(1,2*nAnkle), maxHip*ones(1,2*nHip)];
        mamp = 0*Mamp;
        
        mw = 0*ones(1,4);
        Mw = 5*ones(1,4);
        
        % CPG strucute: (ALSO Symm W_ij = W_ji)
        % 1) every pair of Extensor and reflexor neurons are connected to
        %   each other with symmetric weights.
        % 2) both of the ankles' Extensor neurons are connected to both of the hip neurons
        %
        %  A2_F    A2_E
        %(3) O-----O (4)
        %       /  |
        %      /   |
        %     /    |
        %    /     |
        %   H_F   H_E
        %(5) O-----O (6)  
        %     \    |             %   w = [0  , W12, 0  , 0  , 0  , 0; 
        %      \   |             %        W21, 0  , 0  , 0  , W25, W26;
        %       \  |             %        0  , 0  , 0  , W34, 0  , 0;
        %        \ |             %        0  , 0  , W43, 0  , W45, W46;
        %(1) O-----O (2)         %        0  , 0  , 0  , 0  , 0  , W56;
        %   A1_F    A1_E         %        0  , 0  , 0  , 0  , W65, 0;
        %                        
        %                        w12 = w21 = w34 = w43 = w1  
        %                        w56 = w65  = w2
        %                        w25 = w46 = w3
        %                        w26 = w45 = w4
        
        % Final genome with tau_r + beta (constant tau_u/tau_v ratio) 
        Keys = {'\tau_r', 'beta', 'amp_6n_symm',  '6neuron_taga_like_symm', 'ks_\tau',     'ks_c_6n_symm', 'IC_matsuoka';
                      1 ,      1,             1,                         4,        1 ,                 1 ,            0 };
        Range = {  0.02 ,    0.2,             0,                        mw,      -10 ,                 -0.1*maxAnkle; % Min
                   0.25 ,    2.5,      maxAnkle,                        Mw,       10 ,                  0.1*maxAnkle}; % Max


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
