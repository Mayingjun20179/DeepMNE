%% Scale Similar Matrix by Row %%
%�����Ծ����һ����ÿһ�еĺ�Ϊ1
function W = ScaleSimMat(W)

%scale 
W = W - diag(diag(W));  %diagonal elements must be 0
D = sum(W);    %degree matrix
D(find(D>eps)) = D(find(D>eps)).^(-1);
D1 = diag(D);
W = D1*W;

end