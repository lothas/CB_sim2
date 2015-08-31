FileName = 'SlopeGenomes.mat';
if exist(FileName,'file') == 2
    Data = load(FileName);
    TopN = Data.TopN;
    TopS = Data.TopS;
    Slopes = Data.Slopes;
    Fits = Data.Fits;
    Seqs = Data.Seqs;
    SFits = Data.SFits;
    SSeqs = Data.SSeqs;
    Signals = Data.Signals;
    SSignals = Data.SSignals;
else
    Filenames = {'GA_-170_08_29_22_10.mat';
    %     'GA_-165_08_29_21_51.mat';
        'GA_-160_08_29_21_31.mat';
    %     'GA_-155_08_29_21_03.mat';
        'GA_-150_08_29_18_41.mat';
    %     'GA_-145_08_29_18_06.mat';
        'GA_-140_08_29_17_12.mat';
    %     'GA_-135_08_29_16_17.mat';
        'GA_-130_08_29_15_10.mat';
    %     'GA_-125_08_29_14_43.mat';
        'GA_-120_08_29_13_48.mat';
    %     'GA_-115_08_29_12_31.mat';
        'GA_-110_08_29_11_07.mat';
    %     'GA_-105_08_29_10_18.mat';
        'GA_-100_08_29_00_42.mat';
    %     'GA_-95_08_28_23_20.mat';
        'GA_-90_08_28_21_53.mat';
    %     'GA_-85_08_28_20_10.mat';
        'GA_-80_08_28_17_24.mat';
        'GA_-70_08_28_15_26.mat';
        'GA_-60_08_28_13_08.mat';
        'GA_-50_08_28_10_42.mat';
        'GA_-40_08_27_03_58.mat';
        'GA_-30_08_27_01_00.mat';
        'GA_-20_08_26_21_55.mat';
        'GA_-10_08_26_19_02.mat';
        'GA_0_08_22_17_51.mat';
        'GA_10_08_25_11_22.mat';
        'GA_20_08_25_13_16.mat';
        'GA_30_08_25_17_31.mat';
        'GA_40_08_25_18_39.mat';
        'GA_50_08_26_10_19.mat';
        'GA_60_08_26_11_52.mat';
        'GA_70_08_26_13_37.mat';
        'GA_80_08_26_15_18.mat';
        'GA_83_08_26_18_00.mat'};

    Slopes = [-17:1:-8,-7:1:8,8.3];
    NSlopes = length(Slopes);

    TopN = 20;
    TopS = 3;
    Fits = zeros(TopN,3,NSlopes);
    Seqs = zeros(TopN,14,NSlopes);
    SFits = zeros(TopS,3,NSlopes);
    SSeqs = zeros(TopS,14,NSlopes);

    for i = 1:length(Slopes)
        Data = load(Filenames{i});
        ThisTop = Data.GA.GetTopPop(TopN);
        Fits(:,:,i) = Data.GA.Fit(ThisTop,:,Data.GA.Progress);
        Seqs(:,:,i) = Data.GA.Seqs(ThisTop,:,Data.GA.Progress);

        % Select a few with specific values
        IDs = [];
        MeanPcnt = 1;
        while length(IDs)<TopS
            IDs = find(Fits(:,1,i)>=MeanPcnt*mean(Fits(:,1,i)) & ...
                        Fits(:,2,i)>=MeanPcnt*mean(Fits(:,2,i)),3);
            MeanPcnt = MeanPcnt*0.9;
        end
        SFits(:,:,i) = Fits(IDs,:,i);
        SSeqs(:,:,i) = Seqs(IDs,:,i);
    end

    % Build signal
    Signals = cell(NSlopes,1);
    for i = 1:length(Slopes)
        Signal = [];
        Time = [];
        SigLength = 0;
        for g = 1:TopN
            wSim = deepcopy(Data.GA.Sim);
            wSim = Data.GA.Gen.Decode(wSim,Seqs(g,:,i));
            [ThisTime, ThisSig] = wSim.Con.GetTorqueSig();
            ThisSigLength = length(ThisSig);
            if ThisSigLength>SigLength
                Signal = [Signal;zeros(ThisSigLength-SigLength,2)];
                Signal = Signal+ThisSig;
                SigLength = ThisSigLength;
                Time = ThisTime;
            else
                Signal(1:ThisSigLength,:) = Signal(1:ThisSigLength,:)+ThisSig;
            end
        end
        Signals{i} = Signal;
    end
    
    % Build signal
    SSignals = cell(NSlopes,1);
    for i = 1:length(Slopes)
        Signal = [];
        Time = [];
        SigLength = 0;
        for g = 1:TopS
            wSim = deepcopy(Data.GA.Sim);
            wSim = Data.GA.Gen.Decode(wSim,SSeqs(g,:,i));
            [ThisTime, ThisSig] = wSim.Con.GetTorqueSig();
            ThisSigLength = length(ThisSig);
            if ThisSigLength>SigLength
                Signal = [Signal;zeros(ThisSigLength-SigLength,2)];
                Signal = Signal+ThisSig;
                SigLength = ThisSigLength;
                Time = ThisTime;
            else
                Signal(1:ThisSigLength,:) = Signal(1:ThisSigLength,:)+ThisSig;
            end
        end
        SSignals{i} = Signal;
    end
    
    save(FileName,'TopN','TopS','Slopes','Fits','Seqs','SFits','SSeqs',...
        'Signals','SSignals');
end

figure
hold on
for i = 1:length(Slopes)
    scatter3(Slopes(i)*ones(1,TopN),Fits(:,1,i),Fits(:,2,i));
end
xlabel('Slope');
ylabel('Velocity fitness');
zlabel('Energy fitness');

figure
for i = 1:length(Slopes)
    subplot(1,2,1)
    hold on
    scatter(Slopes(i)*ones(1,TopN),Fits(:,1,i));
    subplot(1,2,2)
    hold on
    scatter(Slopes(i)*ones(1,TopN),Fits(:,2,i));
end
xlabel('Slope');
ylabel('Energy fitness');
subplot(1,2,1)
xlabel('Slope');
ylabel('Velocity fitness');

figure
hold on
for i = 1:length(Slopes)
    plot3(Slopes(i)*ones(1,ThisSigLength),Time,Signals{i}(:,1),'b');
end

figure
hold on
for i = 1:length(Slopes)
    plot3(Slopes(i)*ones(1,ThisSigLength),Time,Signals{i}(:,2),'r');
end

figure
hold on
for i = 1:length(Slopes)
    plot3(Slopes(i)*ones(1,ThisSigLength),Time,SSignals{i}(:,1),'b');
end

figure
hold on
for i = 1:length(Slopes)
    plot3(Slopes(i)*ones(1,ThisSigLength),Time,SSignals{i}(:,2),'r');
end