function rep = MDA_train(xx)
%%%xx表示所有网络构成的结构体矩阵
%d表示融合后网络的维度
% d = floor(0.1 * size(xx{1},2));
d = 20;
%%%%第一层，对所有的输入网络进行编码
Nn = length(xx);  %相似网络的个数
x = [];   %编码层的输出
for i=1:Nn
    dim = size(xx{i},2);
    [sae0,opts0] = BuildNets0_opt(dim);
    sae0.ae{i} = nntrain(sae0.ae{1}, xx{i}, xx{i}, opts0);
    t = nnff(sae0.ae{i}, xx{i}, xx{i});
    xx0 = t.a{2};
    x = [x,xx0(:,2:end)];
end


% 第二，输入一个自编码结构，input x，和opt参数
dim = size(x,2);
[sae,opts,nnsize] = BuildNets_opt(dim,d);
sae = saetrain(sae, x, opts);
rep = GenRep(x, sae, nnsize);  
end
