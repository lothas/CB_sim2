function SlopeSim()
    close all; %clear all;
    Data = load('SlopeGenomes.mat');
    TopN = Data.TopN;
    Seqs = Data.Seqs;
    Slopes = Data.Slopes;
    
    % range of slopes is -16 to +8 = -4 +/-12
    slope = -4*pi/180;
    sloped = slope*180/pi;
    alpha = 12*pi/180;
    
    % Build look-up table
    selid = randsample(1:TopN,1);
    LUT.Slopes = Slopes*pi/180;
    LUT.Freq = reshape(Seqs(selid,5,:),[],1);
    LUT.Amp = reshape(Seqs(1,[6,9,12],:),3,[])';
    LUT.Offset = reshape(Seqs(1,[7, 10, 13],:),3,[])';
    LUT.Duration = reshape(Seqs(1,[8, 11, 14],:),3,[])';
    
    Data = load('GA_-40_08_27_03_58.mat');
    wSim = deepcopy(Data.GA.Sim);
    wSim.Graphics = 1;
    slid = find(abs(Slopes-sloped)==min(abs(Slopes-sloped)), 1, 'first');
    wSim = Data.GA.Gen.Decode(wSim,Seqs(selid,:,slid));
    
    wSim.Con.FBType = 3; % 3 - look-up table
    wSim.Con.LUT = LUT;
    
    % Set up the terrain
    freq = 2*pi/500;
    amp = tan(alpha)/freq;
    wSim.Env = Terrain(1,amp,freq);
    wSim.Env.Set('start_slope',slope,'end_slope',slope);
    wSim.Mod.xS = -3;
    wSim.Mod.yS = wSim.Mod.xS*tan(slope);
    
    wSim = wSim.Init();
    wSim.Mod.LegShift = wSim.Mod.Clearance;
    wSim.Con = wSim.Con.Reset(wSim.IC(wSim.ConCo));
%     wSim.Con = wSim.Con.Adaptation(slope);
    
    wSim.Run();
end