 %%%%%%%%%%%%%整体进行删除
clc
clear

% % %%%%%%%%%%%% 
% load('lnc_dis_inf2019.mat')
% DATA_train = chuli_opt(lnc_dis_inf2019);
% save DATA_train DATA_train;

load DATA_train.mat
%%%%%%%%%%%%%%%%所有交互的推断
cv_flag = 1;
main_cv(cv_flag,DATA_train);

%%%%%计算不同参数下的均值
factor = [50,100];
lata = 10.^[-2:2];
ar = 10.^[-2:2];
cfix = 1:2:9;
option = [];
for i1 = 1:length(factor)
    for i2 = 1:length(lata)
        for i3 = 1:length(ar)
            for i4 = 1:length(cfix)
                option = [option;[factor(i1),lata(i2),ar(i3),cfix(i4)]];
            end
        end
    end
end

load('CVa_result10.mat')
CVA_mean = [];
CVA_std = [];
for i=1:size(option,1)
   flag = ismember(CVa_result10(:,1:4),option(i,:),'rows');  ind = find(flag==1);
   result_mean = [option(i,:),mean(CVa_result10(ind,6:end))];
   result_std = [option(i,:),std(CVa_result10(ind,6:end))];
   CVA_mean = [CVA_mean;result_mean];
   CVA_std = [CVA_std;result_std];   
end

[~,ind] = max(CVA_mean(:,5));
CVA_mean(ind,:)
CVA_std(ind,:)


% load('CVa_result.mat')
% [~,ind] = max(CVa_result(:,5));
% 0.9505    0.9462    0.8787
% 0.0029    0.0026    0.0039

%%%%lncRNA的推断
cv_flag = 2;
load DATA_train.mat
main_cv(cv_flag,DATA_train);
%%%%%%%%%新lncRNA的推断
% top5       top10     top15     top20     top25     top30   top35    top40
% 0.5000    0.4700    0.4667    0.4710    0.4696    0.4733  0.4817    0.4800
%0.1054    0.0598 0.0372    0.0296    0.0230    0.0245    0.0164    0.0168





%%%%%%%%%新疾病的推断
cv_flag = 3;
main_cv(cv_flag,DATA_train);
% top5       top10     top15     top20     top25     top30   top35    top40
% 0.5160    0.4960    0.4640    0.4420 0.4248    0.4000    0.3943    0.3830
% 0.0548    0.0310    0.0435    0.0343 0.0359    0.0263    0.0214    0.0197
