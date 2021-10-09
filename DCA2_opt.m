function S_S = DCA2_opt(S)
%%输入
%S为结构体里面包含多个相似矩阵
%reset_prob 为重启概率
%d为需要得到的特征向量的维数
%%输出
%S_F  为特征矩阵
%Step 1. Randomly Surf to Generate K steps Transition Matrix
Kstep = 3;
alpha = 0.98;

for i=1:length(S)
    Mk{i} = RandSurf(S{i}, Kstep, alpha);
end
% Mk = S;
%Step 2. Get PPMI Matrix
for i=1:length(Mk)
    PPMI{i} = GetPPMIMatrix(Mk{i});   
%     alpha = 1/size(PPMI{i},1);
%     PPMI{i} = log(PPMI{i}+alpha)-log(alpha);
end

%Step 3. Compress the Dimension using SDAE
rep = MDA_train(PPMI);
S_S = KSNS_opt(rep);
end

