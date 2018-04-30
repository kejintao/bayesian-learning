
load SampleTree.mat;
load PubMed2000_binary_test.mat

test_num = 400;
Xid_test = xId_test(1:test_num);
Xcnt_test = xCnt_test(1:test_num);

%%
% for i = 1:length(Xid_test)
%     Xid_test{i} = Xid_test{i} + 1;
% end
% 
% D = length(Xid_test);
% perctest = .1;
% Xid = Xid_test;
% Xcnt = Xcnt_test;
% Xcnt_test = cell(1,D);
% Xid_test = cell(1,D);
% for d = 1:D
%     numW = sum(Xcnt{d});
%     numTest = floor(perctest*numW);
%     [a,b] = sort(rand(1,numW));
%     wordVec = [];
%     for i = 1:length(Xid{d})
%         wordVec = [wordVec Xid{d}(i)*ones(1,Xcnt{d}(i))];
%     end
%     wordTestVec = wordVec(b(1:numTest));
%     wordTrainVec = wordVec(b(numTest+1:end));
%     Xid{d} = unique(wordTrainVec);
%     Xcnt{d} = histc(wordTrainVec,Xid{d});
%     Xid_test{d} = unique(wordTestVec);
%     Xcnt_test{d} = histc(wordTestVec,Xid_test{d});
% end
% num_test = 0;
% for i = 1:length(Xcnt_test)
%     num_test = num_test + sum(Xcnt_test{i});
% end
    D = size(Xid_test,2);
    [llikhood,C_d] = nHDP_test(Xid_test,Xcnt_test,Tree,.1);
    Binary_mean = llikhood/D;
    disp([num2str(Binary_mean)]);
