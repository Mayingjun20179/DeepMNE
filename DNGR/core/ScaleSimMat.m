%% Scale Similar Matrix by Row %%
%相似性矩阵归一化，每一列的和为1
function W = ScaleSimMat(W)

%scale 
W = W - diag(diag(W));  %diagonal elements must be 0
D = sum(W);    %degree matrix
D(find(D>eps)) = D(find(D>eps)).^(-1);
D1 = diag(D);
W = D1*W;

end