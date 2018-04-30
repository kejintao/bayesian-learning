% display vocabulary results. the vocabulary is in a cell called vocab.

[ElnB,ElnPtop,id_parent,id_me] = func_process_tree(Tree,beta0,5); % this only needs to be done once
fid = fopen('./16MarTopics/News5kTopics3.nhdp.txt', 'wt');
Allcases = size(Tree,2);
for idx=1:Allcases   % pick a topic index to show results for
idx_p = find(id_me==id_parent(idx));
idx_c = find(id_parent==id_me(idx));
disp('*** This node ***');
fprintf(fid,'*** This node ***\n');
fprintf(fid,'ThisNode :Index %d Count %f\n',idx,Tree(idx).cnt);
disp(['Count ' num2str(Tree(idx).cnt)]);
%fprintf(fid,'Count\t%f\n',Tree(idx).cnt);

[a,b] = sort(Tree(idx).beta_cnt,'descend');
for w = 1:10
    disp(['   ' vocab{b(w)}]);
    fprintf(fid,'\t%s\n',vocab{b(w)});
end
disp('*** Parent node ***');
fprintf(fid,'*** Parent node ***\n');
if isempty(idx_p)
    disp('No parent');
    fprintf(fid, 'No parent\n');
else
    fprintf(fid,'Parent :Index %d Count %f\n',idx_p,Tree(idx_p).cnt);
    [a,b] = sort(Tree(idx_p).beta_cnt,'descend');
    for w = 1:10
        disp(['   ' vocab{b(w)}]);
        fprintf(fid,'\t%s\n',vocab{b(w)});
    end
end
disp('*** Child nodes ***');
fprintf(fid,'*** Child nodes ***\n');
if isempty(idx_c)
    disp('No children');
    fprintf(fid, 'No children\n');

else
    for i = 1:length(idx_c)
        disp(['Child ' num2str(i) ' : Count ' num2str(Tree(idx_c(i)).cnt) ' : Index ' num2str(idx_c(i))]);
        fprintf(fid, 'Child %d :Index %d Count %f\n',i, idx_c(i),Tree(idx_c(i)).cnt);      
        [a,b] = sort(Tree(idx_c(i)).beta_cnt,'descend');
        for w = 1:10
            disp(['   ' vocab{b(w)}]);
            fprintf(fid,'\t%s\n',vocab{b(w)});
        end
    end
end
end
fclose(fid);
