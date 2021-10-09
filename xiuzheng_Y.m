%%%%关于连接矩阵处理(对于整个矩阵进行处理)
function Y = xiuzheng_Y(Y0,sim_U,sim_V,cv_flag)  %ar为衰减系数
K = 10;
ar = 0.9; 
ind_u0 = find(sum(Y0,2)==0);   %单独的行
ind_v0 = find(sum(Y0,1)==0);    %单独的列
[Nu,Nv] = size(Y0);

%%
sim_U(1:(Nu+1):end) = 0;      %对角线变为0
sim_V(1:(Nv+1):end) = 0; 
sim_U(ind_u0,ind_u0) = 0;       %单独疾病和代谢的相似度为0
sim_V(ind_v0,ind_v0) = 0;

[~,indu] = sort(sim_U,2,'descend');

Yu = zeros(Nu,Nv);
ar = ar.^(0:K-1);
for i=1:Nu    
    near_inu = indu(i,1:K);
    Yu(i,:) = (ar.*sim_U(i,near_inu))*Y0(near_inu,:)/sum(sim_U(i,near_inu));
end

[~,indv] = sort(sim_V,1,'descend');
ar = ar.^(0:K-1);
Yv = zeros(Nu,Nv);
for j=1:Nv    
    near_inv = indv(1:K,j);
    Yv(:,j) = Y0(:,near_inv)*(ar'.*sim_V(near_inv,j))/sum(sim_V(near_inv,j));
end 
Y = (Yu+Yv)/2;
Y = max(Y,Y0);

end