function [llikhood,C_d] = nHDP_test(Xid_test,Xcnt_test,Tree,beta0)
% Written by John Paisley, jpaisley@berkeley.edu

Voc = length(Tree(1).beta_cnt);
tot_tops = length(Tree);
D = length(Xid_test);

% collects statistics for updating the tree
B_up = zeros(tot_tops,Voc);
weight_up = zeros(tot_tops,1);

gamma1 = 5; % top-level DP concentration
gamma2 = 1; % second-level DP concentration
gamma3 = (1/3); % beta switches
gamma4 = (2/3); %

% put info from Tree struct into matrix and vector form

M = zeros(tot_tops,Voc);
cnt_top = zeros(tot_tops,1);
for i = 1:size(M,1)
    if sum(isnan(Tree(i).beta_cnt))
        disp(['NaN Here']);
    Tree(i).beta_cnt = zeros(1,length(Tree(i).beta_cnt));
    disp(['reset']);
    end
    M(i,:) = Tree(i).beta_cnt + beta0;
    M(i,:) = M(i,:)/sum(M(i,:));
    cnt_top(i) = Tree(i).cnt;
end

[ElnB,ElnPtop,id_parent,id_me] = func_process_tree(Tree,beta0,gamma1);

% Family tree indicator matrix. Tree_mat(j,i) = 1 indicates that node j is 
% along the path from node i to the root node
Tree_mat = zeros(tot_tops);
for i = 1:tot_tops
    bool = 1;
    idx = i;
    while bool
        idx = find(id_me==id_parent(idx));
        if ~isempty(idx) %id_me(idx) ~= log(2)
            Tree_mat(idx,i) = 1;
        else
            bool = 0;
        end        
    end
end
level_penalty = psi(gamma3) - psi(gamma3+gamma4) + sum(Tree_mat,1)'*(psi(gamma4) - psi(gamma3+gamma4));
llikhood = 0;
temp_vec = zeros(1,D);
N = 0;
all_weights = zeros(tot_tops,1);
% main loop
for d = 1:D
    N = N + sum(Xcnt_test{d});
    ElnB_d = ElnB(:,Xid_test{d});                                            % pick out words in document for penalty
    ElnV = psi(1) - psi(1+gamma2);
    Eln1_V = psi(gamma2) - psi(1+gamma2);
    ElnP_d = zeros(tot_tops,1); - inf;                                   % -inf removes non-activated topics by giving them zero probability
    ElnP_d(id_parent==log(2)) = ElnV+psi(gamma3)-psi(gamma3+gamma4);    % activate first level of topics

%   % select subtree
%     bool = 1;
%     idx_pick = [];
%     Lbound = [];
%     vec_DPweights = zeros(tot_tops,1);                                  % ElnP_d minus the level penalty
%     while bool
%         idx_active = find(ElnP_d > -inf);                                             % index of active (selected and potential) nodes
%         penalty = ElnB_d(idx_active,:) + repmat(ElnP_d(idx_active),1,length(Xid_test{d}));
%         C_act = penalty;
%         penalty = penalty.*repmat(Xcnt_test{d},size(penalty,1),1);
%         ElnPtop_act = ElnPtop(idx_active);
%         if isempty(idx_pick)
%             score = sum(penalty,2) + ElnPtop_act;
%             [temp,idx_this] = max(score);
%             idx_pick = idx_active(idx_this);                                          % index of selected nodes
%             Lbound(end+1) = temp - ElnPtop_act(idx_this);
%         else
%             temp = zeros(tot_tops,1); 
%             temp(idx_active) = (1:length(idx_active))';
%             idx_clps = temp(idx_pick);                                                % index of selected nodes within active nodes
%             num_act = length(idx_active);
%             vec = max(penalty(idx_clps,:),[],1);
%             C_act = C_act - repmat(vec,num_act,1);
%             C_act = exp(C_act);
%             numerator = C_act.*penalty;
%             numerator = numerator + repmat(sum(numerator(idx_clps,:),1),num_act,1);
%             denominator = C_act + repmat(sum(C_act(idx_clps,:),1),num_act,1);
%             vec = sum(C_act(idx_clps,:).*log(eps+C_act(idx_clps,:)),1);
%             H = log(denominator) - (C_act.*log(C_act+eps) + repmat(vec,num_act,1))./denominator;
%             score = sum(numerator./denominator,2) + ElnPtop_act + H*Xcnt_test{d}';
%             score(idx_clps) = -inf;          
%             [temp,idx_this] = max(score);
%             idx_pick(end+1) = idx_active(idx_this);
%             Lbound(end+1) = temp - ElnPtop_act(idx_this);
%         end
%         idx_this = find(id_parent == id_parent(idx_pick(end)));
%         [t1,t2] = intersect(idx_this,idx_pick);
%         idx_this(t2) = [];
%         vec_DPweights(idx_this) = vec_DPweights(idx_this) + Eln1_V;
%         ElnP_d(idx_this) = ElnP_d(idx_this) + Eln1_V;
%         idx_add = find(id_parent == id_me(idx_pick(end)));
%         vec_DPweights(idx_add) = ElnV;
%         ElnP_d(idx_add) = ElnV + level_penalty(idx_add);
%         bool2 = 1;
%         idx = idx_pick(end);
%         while bool2
%             if ~isempty(idx) %id_me(idx) ~= log(2)
%                 ElnP_d(idx_add) = ElnP_d(idx_add) + vec_DPweights(idx);
%                 idx = find(id_me == id_parent(idx));
%             else
%                 bool2 = 0;
%             end
%         end
%         if length(Lbound) > 5
%             if abs(Lbound(end)-Lbound(end-1))/abs(Lbound(end-1)) < 10^-3 || length(Lbound) == 25
%                 bool = 0;
%             end
%         end
%     end

idx_pick = 1:tot_tops;
    
  % learn document parameters for subtree
    T = length(idx_pick);
    ElnB_d = ElnB(idx_pick,Xid_test{d});
    ElnP_d = 0*ElnP_d(idx_pick) - 1;
    cnt_old = zeros(length(idx_pick),1);
    bool_this = 1;
    num = 0;
    while bool_this
        num = num+1;
        C_d = ElnB_d + repmat(ElnP_d,1,length(Xid_test{d}));
        C_d = C_d - repmat(max(C_d,[],1),T,1);
        C_d = exp(C_d);
        C_d = C_d./repmat(sum(C_d,1),T,1);
        cnt = C_d*Xcnt_test{d}';
%         ElnP_d = func_doc_weight_up(cnt,id_parent(idx_pick),gamma2,gamma3,gamma4,Tree_mat(idx_pick,idx_pick));

        T = length(cnt);
        ElnP_d = zeros(T,1);

        bin_cnt1 = cnt;
        bin_cnt0 =full(sparse(Tree_mat)*cnt);
        Elnbin1 = psi(bin_cnt1+gamma3) - psi(bin_cnt1+bin_cnt0+gamma3+gamma4);
        Elnbin0 = psi(bin_cnt0+gamma4) - psi(bin_cnt1+bin_cnt0+gamma3+gamma4);

        stick_cnt = bin_cnt1+bin_cnt0;
        partition = unique(id_parent);
        for i = 1:length(partition)
            idx = find(id_parent==partition(i));
            t1 = stick_cnt(idx);
            
            ElnP_d(idx) =  psi(t1+cnt_top(idx)/sum(cnt_top(idx))) - psi(1+sum(t1));
        end
        this = ElnP_d + Elnbin1 + full(sparse(Tree_mat')*(Elnbin0 + ElnP_d));
        ElnP_d = this;

        if num > 10
            if sum(abs(cnt-cnt_old))/sum(cnt) < .5*10^-2 || num == 25
                bool_this = 0;
            end
        end
        cnt_old = cnt;
    end

    
    Tree_this = Tree_mat(idx_pick,idx_pick);
    id_par_this = id_parent(idx_pick);
    
    bin_cnt1 = cnt;
    bin_cnt0 = full(sparse(Tree_this)*cnt);
%     idx = find(bin_cnt0<.01);
    Ebin1 = (bin_cnt1+gamma3)./(bin_cnt1+bin_cnt0+gamma3+gamma4);
    Ebin0 = (bin_cnt0+gamma4)./(bin_cnt1+bin_cnt0+gamma3+gamma4);
%     Ebin1(idx) = 1-eps;
%     Ebin0(idx) = eps;
    Elnbin1 = log(Ebin1);
    Elnbin0 = log(Ebin0);
    
    stick_cnt = bin_cnt1+bin_cnt0;
    partition = unique(id_par_this);
    for i = 1:length(partition)
        idx = find(id_par_this==partition(i));
        t1 = stick_cnt(idx);
        weights = log(t1+cnt_top(idx)/sum(cnt_top(idx))) - log(sum(t1)+1);
        ElnP_d(idx) = weights;
    end
    this = ElnP_d + Elnbin1 + full(sparse(Tree_this')*(Elnbin0 + ElnP_d));
    this = this - max(this);
    P_d = exp(this);
    P_d = P_d/sum(P_d);

    vec = P_d'*M(idx_pick,Xid_test{d});
    lliknow = log(vec)*Xcnt_test{d}';
    llikhood = llikhood + lliknow;

end
