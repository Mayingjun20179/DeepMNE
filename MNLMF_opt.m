function scores = MNLMF_opt(Y,Y0,sim_U,sim_V,test_data,option)
cfix = option(4);      %��ǿϵ��
train_set = cfix*Y;                                % c*Y
train_set1 = (cfix-1)*Y + ones(size(Y));   % (c-1)*Y+1

%%%%%%%    ����������˹�ھ����򻯾���
[LM,LD,dsMat,tsMat] =  construct_neighborhood(sim_U, sim_V);   


%%%%%%%    ����AGD����U��V
[U,V] = AGD_optimization(train_set,train_set1,LM,LD,option);    
U = gather(U);  V = gather(V); 
%%%%%��ȫ
[U1,V1] = complete_opt(Y0,dsMat,tsMat,U,V);
pscores = exp(U1*V1')./(1+exp(U1*V1'));
scores = arrayfun(@(x,y)tx_opt(pscores,x,y),test_data(:,1),test_data(:,2));

%%%%%     ����÷�
% tic
% scores = predict_scores(Y0,dsMat,tsMat,U,V,test_data);
% toc

end

function [LM,LD,dsMat,tsMat] = construct_neighborhood(sim_U, sim_V)
%%%%����������˹�����ʱ�򹹽����ھӾ���
%diseaseMat�����ҩ�����ƾ���metMat����ǰб����ƾ���
K1 = 5;           %��������������˹����ʱ�ھӵĸ���
dsMat = sim_U - diag(diag(sim_U));      %%ȥ�Խ��ߵ����ƾ���
tsMat = sim_V - diag(diag(sim_V));
if K1 > 0    %K1Ϊ�ھӵ�ĸ���
    S1 = get_nearest_neighbors(dsMat, K1);   %���ݹ�ʽ��8�������ھӾ���
    LM = laplacian_matrix(S1);
    S2 = get_nearest_neighbors(tsMat, K1);    %���ݹ�ʽ��9�������ھӾ���
    LD = laplacian_matrix(S2);
else  
    LM = laplacian_matrix(dsMat);
    LD = laplacian_matrix(tsMat);
end
end

function X = get_nearest_neighbors(S,K1)
%%%%���ݹ�ʽ��8���ͣ�9�������ھӾ���
[m, n] = size(S);
X = gpuArray.zeros(m, n);
for i = 1:m
    [~,b] = sort(S(i,:),'descend');
    ii = b(1:min(K1,n));     %ֻȡ����ǰ5���ھ�
    X(i,ii) = S(i, ii);
end
end


function L = laplacian_matrix(S)
%%%���ݹ�ʽ��10��������Ǹ���ʽ����������˹���� 
%%% SΪȡ�ھӺ�����ƾ���
x = sum(S);
y = sum(S');
L = 0.5*(diag(x+y) - (S+S')); % neighborhood regularization matrix
end




function [U,V] = AGD_optimization(train_set,train_set1,LM,LD,option)
train_set = gpuArray(train_set);
train_set1 = gpuArray(train_set1);
LM = gpuArray(LM);
LD = gpuArray(LD);

%%%����
theta = 1;
num_factors = option(1);
[nm,nd] = size(train_set);
max_iter=100;   %AGD�㷨�ĵ�������
%%%%%%����U��V
%���ɳ�ʼֵ
seed = 1;
randn('state',seed)
U = sqrt(1/num_factors)*gpuArray(randn(nm,num_factors));
randn('state',seed)
V = sqrt(1/num_factors)*gpuArray(randn(nd,num_factors));

dg_sum = gpuArray.zeros(size(U));
tg_sum = gpuArray.zeros(size(V));
last_log = log_likelihood(U,V,LM,LD,option,train_set,train_set1);   %%���չ�ʽ��12�����������Ȼ����ֵ���෴��
for t =  1:max_iter
    %%%%����U
    dg =  deriv_opt(train_set,train_set1,U,V,'disease',LM,LD,option); %��ʽ��13���й���U�ĵ������෴��
    dg_sum = dg_sum + dg.^2;   %����8�еĦ�
    vec_step_size = theta*gpuArray.ones(size(dg_sum))./ sqrt(dg_sum);   %����������
    U = U + vec_step_size .* dg;    %����9
    %%%%����V
    tg = deriv_opt(train_set,train_set1,U,V,'met',LM,LD,option);  
    tg_sum = tg_sum + tg.^2;
    vec_step_size = theta*gpuArray.ones(size(tg_sum)) ./ sqrt(tg_sum);
    V = V + vec_step_size .* tg;
    %%%%
    curr_log = log_likelihood(U,V,LM,LD,option,train_set,train_set1);   %%���ݵ�ǰU��V���㹫ʽ��12�����෴��
    delta_log = (curr_log-last_log)/abs(last_log);   %�ı��������
    if abs(delta_log) < 1e-5
        break;
    end
    last_log = curr_log; 
end
end

function vec_deriv = deriv_opt(train_set,train_set1,U,V,name,LM,LD,option)

lata = option(2);
ar = option(3);

%%%���չ�ʽ��13������U��V�ĵ������෴��
if strcmp(name,'disease')==1
    vec_deriv = train_set*V;   %c*Y*V
else
    vec_deriv = train_set'*U;  %c*Y'*V
end
A = exp(U*V');
A = A./(A + ones(size(train_set)));   %%����P
A = train_set1.* A;     %%(��c-1��Y+1).*P 
if strcmp(name,'disease') == 1
    vec_deriv = vec_deriv - A * V;
    vec_deriv = vec_deriv - (ar*U+lata*LM*U);   %%%
else
    vec_deriv = vec_deriv - A'*U;
    vec_deriv =  vec_deriv - (ar*V+lata*LD*V);
end
end

function loglik = log_likelihood(U,V,LM,LD,option,train_set,train_set1)


%%%%%�˹�ʽΪ��ʽ��12�����෴��
%%%%intmat ����c*Y
lata = option(2);
ar = option(3);
Z = U*V';
loglik = sum(sum(train_set1.*log(1+exp(Z))-train_set.*Z))+...
    1/2*trace(U'*(ar*eye(size(LM))+lata*LM)*U)+...
    1/2*trace(V'*(ar*eye(size(LD))+lata*LD)*V);
    

end


%%%%%%%%%%%%%%������ȫ
function [U1,V1] = complete_opt(train_set,sim_U,sim_V,U,V)
%K2 ��ʾû��������������ҩ���б��������ʱʹ�õ��ھӸ���
%test_data  Ϊ���еľ��󣬵�һ�б�ʾҩ����������ڶ��б�ʾ�б������
%N��ʾ����ǰҩ���б�û�����κΰб��ҩ��������ʱ��ѡȡN�����ڽ��ھ������ļ�Ȩƽ����Ϊ����������
%%%%train_diseases��ʾ����ȫ0����Щ�е��±�   
%%%%train_mets��ʾ����ȫ0����Щ�е��±� 
U1 = U;
V1 = V;
K = 10;           %��ʾû��������������ҩ���б��������ʱʹ�õ��ھӸ���

flagu = sum(train_set,2);    inu0 = find(flagu==0);    inuu = find(flagu>0);   
flagv = sum(train_set,1);   inv0 = find(flagv==0);    inuv = find(flagv>0);     

ar = 0.9;  ar = ar.^(0:K-1);  %��Ȩ
%%%%%��U���в�ȫ
for i=1:length(inu0)
    [~,sortu] = sort(sim_U(inu0(i),inuu),'descend');
    inuk = inuu(sortu(1:K));   %ѡȡ���������������
%     U1(inu0(i),:) = (sim_U(inu0(i),inuk))*U(inuk,:)/sum(sim_U(inu0(i),inuk));
    U1(inu0(i),:) = (ar.*sim_U(inu0(i),inuk))*U(inuk,:)/sum(ar.*sim_U(inu0(i),inuk));
end

%%%%%��U���в�ȫ
for i=1:length(inv0)   
    [~,sortv] = sort(sim_V(inv0(i),inuv),'descend');
    invk = inuv(sortv(1:K));   %ѡȡ���������������
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