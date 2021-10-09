function main_cv(cv_flag,DATA)

%生成参数
%option(1)  子空间维度，
%option(2) 拉格朗日乘子正则化参数
%option(3)  这个是潜在空间特征正则化参数
%option(4)  重要性水平参数

% seed = 1:100:1000;
% factor = [50,100];
% lata = 10.^[-2:2];
% ar = 10.^[-2:2];
% cfix = 1:2:9;
% option = [];
% for i1 = 1:length(factor)
%     for i2 = 1:length(lata)
%         for i3 = 1:length(ar)
%             for i4 = 1:length(cfix)
%                 for i5 = 1:length(seed)
%                     option = [option;[factor(i1),lata(i2),ar(i3),cfix(i4),seed(i5)]];
%                 end
%             end
%         end
%     end
% end
seed = 1:100:1000;
option = [100.0000   10.0000   10.0000    3.0000];
option = [repmat(option,10,1),seed'];


%%生成交叉实验的测试集
cv = 5;  %cv -折交叉实验
%NRLMF算法计算结果
[N,NCW] = size(option);
result = zeros(N,NCW+11);
parfor i = 1:N
    i
    seeds = option(i,5);
    jieguo = cv_seed(DATA,seeds,cv_flag,cv,option(i,1:4));
    jieguo = [option(i,:),jieguo];
    result(i,:) = jieguo;
end

if cv_flag == 1
    CVa_result10 = result; 
    save  CVa_result10 CVa_result10;
elseif cv_flag == 2
    CVL_result10 = result; 
    save  CVL_result10 CVL_result10;
elseif cv_flag == 3
    CVd_result10 = result; 
    save  CVd_result10 CVd_result10; 
end
   

end



%%%%计算在不同的种子下的预测结果
function result = cv_seed(DATA,seeds,cv_flag,nF,option)
Ns = length(seeds);

result = 0;
for i=1:Ns
    cv_data = cross_validation(DATA,seeds(i),cv_flag,nF);
    result0 = five_cross(DATA,cv_data,option,cv_flag);
    result = result + result0;
end
result = result/Ns;
[option,result]
end


%%%%%%%%%%%%%%%%%%%%%%生成交叉实验的数据集
%%%%生成不同种子时的判断矩阵（测试集的索引为0，其他为1）
%%%%测试集的索引，测试集的标签
%intMat 交互矩阵
%seed   表示生成随机数的种子
%cv_flag     表示选择哪一种交叉实验
%cv_flag =1  表示对所有交互行交叉实验
%cv_flag = 2 表示对行交叉实验
%cv_flag = 3 表示对列交叉实验
%cv    表示交叉实验的次数（5折，10折）
function cv_data = cross_validation(DATA,seeds,cv_flag,cv)

intMat = DATA.interaction;
Ns = length(seeds);
cv_data = cell(Ns,cv);  %cv_data{i,j}表示种子为seeds(i)时，第j个交叉实验对应的判别矩阵W
%测试点的位置 test_data，测试点的标签 test_label
if cv_flag<=3
    for i= 1:Ns
        [Nu,Nv] = size(intMat);      [Nx1,Ny1] = find(intMat==1);  Nuv = length(Nx1);
        rand('state',seeds(i))
        if cv_flag == 1
            index = randperm(Nuv)';
        elseif cv_flag == 2
            index = randperm(Nu)';
        elseif cv_flag == 3
            index = randperm(Nv)';
        end
        step = floor(length(index)/cv);   %
        for j = 1:cv   %计算每次交叉实验对应的测试集
            if j < cv
                ii = index((j-1)*step+1:j*step);
            else
                ii = index((j-1)*step+1:end);
            end
            if cv_flag == 1    %对所有交互交叉
                test_data = [Nx1(ii),Ny1(ii)];   %这是1/5的1的位置
                [indx0,indy0] = find(intMat==0); %这是全部的0的位置
                %%%%%修改为1和0的个数一样多的
                N1 = size(test_data,1);
                N0 = length(indx0);  indx00 = randperm(N0); 
                test_data = [test_data;[indx0(indx00(1:N1)),indy0(indx00(1:N1))]];   %这是1/5的1的位置
%                 test_data = [test_data;[indx0,indy0]];   %这是1/5的1的位置
            elseif cv_flag == 2   %对所有行交叉
                test_data=[];
                for k=1:length(ii)
                    test_data = [test_data;[ii(k)*ones(Nv,1),[1:Nv]']];
                end
            elseif cv_flag==3     %对所有列交叉
                test_data=[];
                for k=1:length(ii)
                    test_data = [test_data;[[1:Nu]',ii(k)*ones(Nu,1)]];
                end
            end
            test_label = [];
            W = ones(size(intMat));
            for s=1:size(test_data,1)
                test_label = [test_label;intMat(test_data(s,1),test_data(s,2))];
                W(test_data(s,1),test_data(s,2)) = 0;
            end
            cv_data{i,j} = {W, test_data, test_label};
        end
    end
end

end


