classdef MatsuokaML
    %MATSUOKAML Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nNeurons = 4;
        NN = [];
        SVM = [];
        normParams = [];
        TrTaRatio = 5;
        
        doPlot = 0;
        
        tStep = 0.01;
        tEnd = 6;
        
        absTol = 1e-8;
        relTol = 1e-7;
        
        perLim = [];
        perLimOut = [];
    end
    
    methods
    end
    
end

