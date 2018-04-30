
for i=1:length(xCnt)
    one_doc = xCnt{1, i};
    for j=1:length(one_doc)
        if one_doc(1, j) == 0
            fprintf('one zero\n')
        end
    end
end