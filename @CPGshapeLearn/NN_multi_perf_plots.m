function NN_multi_perf_plots(obj)
% plot the performance of each NN when we train a lot

% TODO- fix bugs here!!

NNs = obj.multi_NN.NNs;
NNs_tr = obj.multi_NN.net_perf;
outNames = obj.multi_NN.outNames;

num_of_NNs = length(NNs);

for j=1:num_of_NNs
    subplot(3,3,j)
    plotperform(TR)
end

for j=1:num_of_NNs
    subplot(3,3,j)
    NN = NNs{1,j};
    tr = NNs_tr{1,j};
    test_ind = tr.testMask{1,1};
    test_ind = ~isnan(test_ind);
    NN_out = NN(obj.sampl);
    NN_test_out = NN_out(test_ind);
    NN_test_targ = obj.targ(test_ind);
    plotregression(NN_test_targ, NN_test_out, 'Testing')
end
end

