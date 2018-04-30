% initialize / use a large subset of documents (e.g., 10,000) contained in Xid and Xcnt to initialize
num_topics = [20 3 3];
scale = 1000;
% load Nips1000_binary_train.mat

load PubMed2000_binary_train.mat
train_num = 1500;
xId = xId_train(1:train_num);
xCnt = xCnt_train(1:train_num);

%%
t1=cputime;
%  for i = 1:length(xWords)
%      Xid{i} = xWords{i} + 1; I have deal with it when producing them
%  end
%Xid = xId;
%Xcnt = xCnt;
Tree = nHDP_init(xId,xCnt,num_topics,scale);
for i = 1:length(Tree)
    if Tree(i).cnt == 0
        Tree(i).beta_cnt(:) = 0;
    end
    vec = gamrnd(ones(1,length(Tree(i).beta_cnt)),1);
    Tree(i).beta_cnt = .95*Tree(i).beta_cnt + .05*scale*vec/sum(vec);
end

% main loop / to modify this, at each iteration send in a new subset of docs
% contained in Xid_batch and Xcnt_batch
beta0 = .1; % this parameter is the Dirichlet base distribution and can be played with
for i = 1:1000
    iter_string = sprintf('__________iteration: %d _____________',i);
    disp(iter_string)
%     disp(['__________iteration_____________']);
    [a,b] = sort(rand(1,length(xId)));
    rho = (1+i)^-.75; % step size can also be played with
    Xid_batch = xId(b(1:500));
    Xcnt_batch = xCnt(b(1:500));
    Tree = nHDP_step(Xid_batch,Xcnt_batch,Tree,scale,rho,beta0);
end
TimeNips1k=cputime-t1;
disp(['Finished and totaltime is : ' num2str(TimeNips1k/60)]);
save SampleTree Tree
