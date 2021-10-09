function scores = MNLMF_opt(Y,Y0,sim_U,sim_V,test_data,option)
cfix = option(4);      %增强系数
train_set = cfix*Y;                                % c*Y
train_set1 = (cfix-1)*Y + ones(size(Y));   % (c-1)*Y+1

%%%%%%%    构建拉普拉斯邻居正则化矩阵
[LM,LD,dsMat,tsMat] =  construct_neighborhood(sim_U, sim_V);   


%%%%%%%    利用AGD计算U和V
[U,V] = AGD_optimization(train_set,train_set1,LM,LD,option);    
U = gather(U);  V = gather(V); 
%%%%%补全
[U1,V1] = complete_opt(Y0,dsMat,tsMat,U,V);
pscores = exp(U1*V1')./(1+exp(U1*V1'));
scores = arrayfun(@(x,y)tx_opt(pscores,x,y),test_data(:,1),test_data(:,2));

%%%%%     计算得分
% tic
% scores = predict_scores(Y0,dsMat,tsMat,U,V,test_data);
% toc

end

function [LM,LD,dsMat,tsMat] = construct_neighborhood(sim_U, sim_V)
%%%%计算拉普拉斯矩阵的时候构建的邻居矩阵
%diseaseMat这个是药物相似矩阵，metMat这个是靶标相似矩阵
K1 = 5;           %计算邻域拉普拉斯矩阵时邻居的个数
dsMat = sim_U - diag(diag(sim_U));      %%去对角线的相似矩阵
tsMat = sim_V - diag(diag(sim_V));
if K1 > 0    %K1为邻居点的个数
    S1 = get_nearest_neighbors(dsMat, K1);   %根据公式（8）构建邻居矩阵
    LM = laplacian_matrix(S1);
    S2 = get_nearest_neighbors(tsMat, K1);    %根据公式（9）构建邻居矩阵
    LD = laplacian_matrix(S2);
else  
    LM = laplacian_matrix(dsMat);
    LD = laplacian_matrix(tsMat);
end
end

function X = get_nearest_neighbors(S,K1)
%%%%根据公式（8）和（9）构建邻居矩阵
[m, n] = size(S);
X = gpuArray.zeros(m, n);
for i = 1:m
    [~,b] = sort(S(i,:),'descend');
    ii = b(1:min(K1,n));     %只取排名前5的邻居
    X(i,ii) = S(i, ii);
end
end


function L = laplacian_matrix(S)
%%%根据公式（10）下面的那个公式构建拉普拉斯矩阵 
%%% S为取邻居后的相似矩阵
x = sum(S);
y = sum(S');
L = 0.5*(diag(x+y) - (S+S')); % neighborhood regularization matrix
end




function [U,V] = AGD_optimization(train_set,train_set1,LM,LD,option)
train_set = gpuArray(train_set);
train_set1 = gpuArray(train_set1);
LM = gpuArray(LM);
LD = gpuArray(LD);

%%%参数
theta = 1;
num_factors = option(1);
[nm,nd] = size(train_set);
max_iter=100;   %AGD算法的迭代次数
%%%%%%生成U和V
%生成初始值
seed = 1;
randn('state',seed)
U = sqrt(1/num_factors)*gpuArray(randn(nm,num_factors));
randn('state',seed)
V = sqrt(1/num_factors)*gpuArray(randn(nd,num_factors));

dg_sum = gpuArray.zeros(size(U));
tg_sum = gpuArray.zeros(size(V));
last_log = log_likelihood(U,V,LM,LD,option,train_set,train_set1);   %%按照公式（12）计算对数似然估计值的相反数
for t =  1:max_iter
    %%%%更新U
    dg =  deriv_opt(train_set,train_set1,U,V,'disease',LM,LD,option); %公式（13）中关于U的导数的相反数
    dg_sum = dg_sum + dg.^2;   %步骤8中的φ
    vec_step_size = theta*gpuArray.ones(size(dg_sum))./ sqrt(dg_sum);   %迭代步长γ
    U = U + vec_step_size .* dg;    %步骤9
    %%%%更新V
    tg = deriv_opt(train_set,train_set1,U,V,'met',LM,LD,option);  
    tg_sum = tg_sum + tg.^2;
    vec_step_size = theta*gpuArray.ones(size(tg_sum)) ./ sqrt(tg_sum);
    V = V + vec_step_size .* tg;
    %%%%
    curr_log = log_likelihood(U,V,LM,LD,option,train_set,train_set1);   %%根据当前U，V计算公式（12）的相反数
    delta_log = (curr_log-last_log)/abs(last_log);   %改变的相对误差
    if abs(delta_log) < 1e-5
        break;
    end
    last_log = curr_log; 
end
end

function vec_deriv = deriv_opt(train_set,train_set1,U,V,name,LM,LD,option)

lata = option(2);
ar = option(3);

%%%按照公式（13）计算U和V的导数的相反数
if strcmp(name,'disease')==1
    vec_deriv = train_set*V;   %c*Y*V
else
    vec_deriv = train_set'*U;  %c*Y'*V
end
A = exp(U*V');
A = A./(A + ones(size(train_set)));   %%计算P
A = train_set1.* A;     %%(（c-1）Y+1).*P 
if strcmp(name,'disease') == 1
    vec_deriv = vec_deriv - A * V;
    vec_deriv = vec_deriv - (ar*U+lata*LM*U);   %%%
else
    vec_deriv = vec_deriv - A'*U;
    vec_deriv =  vec_deriv - (ar*V+lata*LD*V);
end
end

function loglik = log_likelihood(U,V,LM,LD,option,train_set,train_set1)


%%%%%此公式为公式（12）的相反数
%%%%intmat 就是c*Y
lata = option(2);
ar = option(3);
Z = U*V';
loglik = sum(sum(train_set1.*log(1+exp(Z))-train_set.*Z))+...
    1/2*trace(U'*(ar*eye(size(LM))+lata*LM)*U)+...
    1/2*trace(V'*(ar*eye(size(LD))+lata*LD)*V);
    

end


%%%%%%%%%%%%%%特征补全
function [U1,V1] = complete_opt(train_set,sim_U,sim_V,U,V)
%K2 表示没有与其他相连的药物或靶标计算连接时使用的邻居个数
%test_data  为两列的矩阵，第一列表示药物的索引，第二列表示靶标的索引
%N表示若当前药物或靶标没有与任何靶标和药物相连的时候，选取N个最邻近邻居特征的加权平均最为它的特征。
%%%%train_diseases表示不是全0的那些行的下标   
%%%%train_mets表示不是全0的那些列的下标 
U1 = U;
V1 = V;
K = 10;           %表示没有与其他相连的药物或靶标计算连接时使用的邻居个数

flagu = sum(train_set,2);    inu0 = find(flagu==0);    inuu = find(flagu>0);   
flagv = sum(train_set,1);   inv0 = find(flagv==0);    inuv = find(flagv>0);     

ar = 0.9;  ar = ar.^(0:K-1);  %加权
%%%%%对U进行补全
for i=1:length(inu0)
    [~,sortu] = sort(sim_U(inu0(i),inuu),'descend');
    inuk = inuu(sortu(1:K));   %选取的邻域的特征索引
%     U1(inu0(i),:) = (sim_U(inu0(i),inuk))*U(inuk,:)/sum(sim_U(inu0(i),inuk));
    U1(inu0(i),:) = (ar.*sim_U(inu0(i),inuk))*U(inuk,:)/sum(ar.*sim_U(inu0(i),inuk));
end

%%%%%对U进行补全
for i=1:length(inv0)   
    [~,sortv] = sort(sim_V(inv0(i),inuv),'descend');
    invk = inuv(sortv(1:K));   %选取的邻域的特征索引
%     V1(inv0(i),:) = (sim_V(inv0(i),invk))*V(invk,:)/sum(sim_V(inv0(i),invk)); 
    V1(inv0(i),:) = (ar.*sim_V(inv0(i),invk))*V(invk,:)/sum(ar.*sim_V(inv0(i),invk));     
end

end
% 
% 
% 
% 
% 
% 
% 
function c = tx_opt(A,b,c)
c = A(b,c);
end