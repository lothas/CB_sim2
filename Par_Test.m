MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.01; % 0.01
MML.tEnd = 15; % 15

% Get a periodic sequence
while 1
    seq = MML.Gen.RandSeq(1);
    [out, ~, ~] = MML.runSim(seq);
    if all(~isnan(out.periods)) && out.perOK == 1
        break
    end
end
MML.tEnd = min(10*max(out.periods),15);

% Test period relation to tau and b
points = 20;
points2 = 3;
b = linspace(seq(2)/3,seq(2)*3,points);
b2 = linspace(seq(2)/3,seq(2)*3,points2);
tau = linspace(seq(1)/3,seq(1)*3,points);
tau2 = linspace(seq(1)/3,seq(1)*3,points2);
MML.Gen.Range(:,1:2) = [min(tau), min(b); max(tau), max(b)];

out_tau = repmat(out, points, points2); out_b = out_tau;
per_tau = zeros(points,points2); per_b = per_tau;
perOK_tau = per_tau; perOK_b = per_tau;

parfor i = 1:points
    for j = 1:points2
        seq_tau = seq;
        seq_tau(1) = tau(i);
        seq_tau(2) = b2(j);
%         [out_tau(i,j), sim, signal] = MML.runSim(seq_tau);
        [out_tau(i,j), ~, ~] = MML.runSim(seq_tau);
        if any(isnan(out_tau(i,j).periods))
            per_tau(i,j) = NaN;
        else
            per_tau(i,j) = max(out_tau(i,j).periods);
        end
        perOK_tau(i,j) = out_tau(i,j).perOK;
        
        seq_b = seq;
        seq_b(1) = tau2(j);
        seq_b(2) = b(i);
%         [out_tau(i,j), sim, signal] = MML.runSim(seq_b);
        [out_b(i,j), ~, ~] = MML.runSim(seq_b);
        if any(isnan(out_b(i,j).periods))
            per_b(i,j) = NaN;
        else
            per_b(i,j) = max(out_b(i,j).periods);
        end
        perOK_b(i,j) = out_b(i,j).perOK;
    end
end

% tau_leg = cellstr(num2str(tau2(:)))';
% b_leg = cellstr(num2str(b2(:)))';
tau_leg = cellfun(@(x)sprintf('tau = %.3f',x), ...
    num2cell(tau2)', 'UniformOutput', 0);
b_leg = cellfun(@(x)sprintf('b = %.3f',x), ...
    num2cell(b2)', 'UniformOutput', 0);

figure
hold on
h = zeros(1,points2);
for i = 1:points2
    good_tau2 = perOK_tau(:,i) == 1;
    h(i) = plot(tau, per_tau(:,i));
    scatter(tau(good_tau2), per_tau(good_tau2,i), '*')
end
legend(h, b_leg)
xlabel('\tau')
ylabel('Period')

figure
hold on
h = zeros(1,points2);
for i = 1:points2
    good_b2 = perOK_b(:,i) == 1;
    h(i) = plot(b, per_b(:,i));
    scatter(b(good_b2), per_b(good_b2,i), '*')
end
legend(h, tau_leg)
xlabel('b')
ylabel('Period')

% Test period relation to all weights
w_ids = 7:18;
pointsW = 50;
k = logspace(-1, 0.5,pointsW);
weights = seq(w_ids);
maxw = max(weights)*max(k);
minw = min(min(weights)*min(k), min(weights)*max(k));
MML.Gen.Range(:,w_ids) = repmat([minw; maxw], 1, length(w_ids));

out_W = repmat(out, pointsW, 1);
per_W = zeros(pointsW,1); perOK_W = per_W;
parfor i = 1:pointsW
    seqW = seq;
    seqW(w_ids) = k(i)*seqW(w_ids);
    [out_W(i), ~, ~] = MML.runSim(seqW);
%     [out_W(i), sim, signal] = MML.runSim(seqW);
    
    if any(isnan(out_W(i).periods))
        per_W(i) = NaN;
    else
        per_W(i) = max(out_W(i).periods);
    end
    perOK_W(i) = out_W(i).perOK;
end

good_W2 = perOK_W == 1;

figure;
semilogx(k,per_W);
hold on
semilogx(k(good_W2),per_W(good_W2),'*');
xlabel('k'); ylabel('period');
title('Effect of increasing all weights times k')

save('ParamTest', 'seq', 'points', 'points2', 'tau', 'tau2', ...
    'b', 'b2', 'out_tau', 'out_b', 'per_tau', 'per_b', ...
    'perOK_tau', 'perOK_b', 'pointsW', 'k', 'out_W', 'per_W', ...
    'perOK_W');

% Test period relation to weights going into first neuron
w_ids = 7:18;
wk_ids = 7:9;
points = 50;
k = logspace(-1, 0.5,points);
weights = seq(w_ids);
maxw = max([max(weights)*max(k), MML.Gen.Range(2,w_ids)]);
minw = min([min(weights)*min(k), min(weights)*max(k), ...
    MML.Gen.Range(1,w_ids)]);
MML.Gen.Range(:,w_ids) = repmat([minw; maxw], 1, length(w_ids));

out = cell(points,1); sim = out; signal = out;
periods = zeros(points,1);
parfor i = 1:points
    seq2 = seq;
    seq2(wk_ids) = k(i)*seq2(wk_ids);
    [out{i}, sim{i}, signal{i}] = MML.runSim(seq2);
    
    if any(isnan(out{i}.periods))
        periods(i) = NaN;
    else
        periods(i) = max(out{i}.periods);
    end
end

out = cell2mat(out);
good_ids = find([out.perOK]);
gPers = max([out(good_ids).periods]);
gk = k(good_ids);

figure;
semilogx(gk,gPers,'*-');
hold on
semilogx(k,periods,'o');
xlabel('k'); ylabel('period');

% Test period relation to weights coming from first neuron
w_ids = 7:18;
wk_ids = [10,13,16];
points = 50;
k = logspace(-1, 0.5,points);
weights = seq(w_ids);
maxw = max([max(weights)*max(k), MML.Gen.Range(2,w_ids)]);
minw = min([min(weights)*min(k), min(weights)*max(k), ...
    MML.Gen.Range(1,w_ids)]);
MML.Gen.Range(:,w_ids) = repmat([minw; maxw], 1, length(w_ids));

out = cell(points,1); sim = out; signal = out;
periods = zeros(points,1);
parfor i = 1:points
    seq2 = seq;
    seq2(wk_ids) = k(i)*seq2(wk_ids);
    [out{i}, sim{i}, signal{i}] = MML.runSim(seq2);
    
    if any(isnan(out{i}.periods))
        periods(i) = NaN;
    else
        periods(i) = max(out{i}.periods);
    end
end

out = cell2mat(out);
good_ids = find([out.perOK]);
gPers = max([out(good_ids).periods]);
gk = k(good_ids);

figure;
semilogx(gk,gPers,'*-');
hold on
semilogx(k,periods,'o');
xlabel('k'); ylabel('period');