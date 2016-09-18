classdef MatsuokaML
    %MATSUOKAML Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nNeurons = 4;
        NN = [];
        SVM = [];
        normParams = [];
        
        doPlot = 0;
        
        tStep = 0.025;
        tEnd = 50;
        
        absTol = 1e-8;
        relTol = 1e-7;
        
        perLim = [];
        perLimOut = [];
        
        Sim = [];
        Gen = [];
        
        selected_genes = {'weights'};
    end
    
    methods
        function obj = MatsuokaML()
            genome_file = 'MatsuokaGenome.mat';
            load(genome_file);
            Keys(:,strcmp(Keys(1,:),'IC_matsuoka')) = []; %#ok<NODEF>
            obj.Gen = Genome(Keys, Range);
            
            % Prepare Sim of Matsuoka CPG so Gen can decode and use a
            % genetic sequence
            obj.Sim.Mod.SetKeys = {};
            obj.Sim.Env.SetKeys = {};
            
            % Initialize the controller
            obj.Sim.Con = Matsuoka;
            obj.Sim.Con.startup_t = 1.5; % Give some time for the neurons to converge
            % before applying a torque
            obj.Sim.Con.FBType = 0; % no slope feedback
            obj.Sim.Con.nPulses = N;
            obj.Sim.Con.stDim = 4*N;
            obj.Sim.Con = obj.Sim.Con.SetOutMatrix([nAnkle,nHip]);
            obj.Sim.Con.MinSat = [-maxAnkle,-maxHip];
            obj.Sim.Con.MaxSat = [ maxAnkle, maxHip];
            
%             obj.Sim = obj.GenDecode(obj.Sim, [0.1, 6.0, 2.2, 5.9, 4.7, 1.52, 1.94, 0.63, 1.22, 0.1, -0.17, 0, 1.88, 3.67, 0.89, 0, 1.57, 0, 0, 0, 0, 0]);
        end
        
    end
    
end

