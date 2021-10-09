function result = five_cross(DATA,cv_data,option,cv_flag)
interaction = DATA.interaction;
cv = length(cv_data);
result = 0;

for k=1:cv
    %%%%%%%%%%%%%%%%%%%%%%%%%
    train_set = interaction.*cv_data{k}{1};    %��ʵ�����Ӿ������в��Լ�λ�ñ�Ϊ0
    Y = xiuzheng_Y(train_set,DATA.lnc_sim,DATA.dis_sim);
    %%%%%%����lncRNA���ƾ���
    KSNS_SL = KSNS_opt(Y);  
    lnc_sim = KSNS_SL;
    %%%���㼲�������ƾ���(�����ø�˹���ںϣ�
    KSNS_Sd = KSNS_opt(Y');
    dis_sim = KSNS_Sd;
    %%%%MNLMF
    scores = MNLMF_opt(Y,train_set,lnc_sim,dis_sim,cv_data{k}{2},option);
    
    %Result Evaluation
    result0 = evaluate_opt(scores,cv_data{k}{3});
    result = result + result0;    
end
result = result/cv
end


function result = evaluate_opt(score,label)
%%%%����TN,TP,FN,FP
sort_predict_score=unique(sort(score));
score_num = length(sort_predict_score);
Nstep = min(score_num,2000);
threshold=sort_predict_score(ceil(score_num*(1:Nstep)/(Nstep+1)));

threshold=threshold';
threshold = threshold(end:-1:1);
threshold_num=length(threshold);
TN=zeros(threshold_num,1);
TP=zeros(threshold_num,1);
FN=zeros(threshold_num,1);
FP=zeros(threshold_num,1);

for i=1:threshold_num
    tp_index=logical(score>=threshold(i) & label==1);
    TP(i,1)=sum(tp_index);
    
    tn_index=logical(score<threshold(i) & label==0);
    TN(i,1)=sum(tn_index);
    
    fp_index=logical(score>=threshold(i) & label==0);
    FP(i,1)=sum(fp_index);
    
    fn_index=logical(score<threshold(i) & label==1);
    FN(i,1)=sum(fn_index);
end


%%%%%����AUPR
SEN=TP./(TP+FN);
PRECISION=TP./(TP+FP);
recall=SEN;
x=recall;
y=PRECISION;
[x,index]=sort(x);
y=y(index,:);
x = [0;x];  y = [1;y];
x(end+1,1)=1;  y(end+1,1)=0;
AUPR=0.5*x(1)*(1+y(1));
for i=1:threshold_num
    AUPR=AUPR+(y(i)+y(i+1))*(x(i+1)-x(i))/2;
end
AUPR_xy = [x,y];

%%%%%����AUC
AUC_x = FP./(TN+FP);      %FPR
AUC_y = TP./(TP+FN);      %tpR
[AUC_x,ind] = sort(AUC_x);
AUC_y = AUC_y(ind);
AUC_x = [0;AUC_x];
AUC_y = [0;AUC_y];
AUC_x = [AUC_x;1];
AUC_y = [AUC_y;1];

AUC = 0.5*AUC_x(1)*AUC_y(1)+sum((AUC_x(2:end)-AUC_x(1:end-1)).*(AUC_y(2:end)+AUC_y(1:end-1))/2);

AUCxy = [AUC_x(:),AUC_y(:)];

%%%%%��������ָ��
temp_accuracy=(TN+TP)/length(label);   %%׼ȷ��
temp_sen=TP./(TP+FN);    %%��ʵ�������ӣ�Ԥ����ȷ����ȷ��
recall = temp_sen;
temp_spec=TN./(TN+FP);   %%��ʵ�����ӣ�Ԥ����ȷ��
temp_precision=TP./(TP+FP); %%Ԥ�������ӵ���ȷ��
temp_f1=2*recall.*temp_precision./(recall+temp_precision);
[~,index]=max(temp_f1);
%%%%����F1��ߵ�����ֵ��
f1=temp_f1(index);
%%%%����÷�ǰ10��ǰ15��ǰ20���ٻ���
[~,ind] = sort(score,'descend');

precision_top5 = sum(label(ind(1:5)))/5;
precision_top10 = sum(label(ind(1:10)))/10;
precision_top15 = sum(label(ind(1:15)))/15;
precision_top20 = sum(label(ind(1:20)))/20;
precision_top25 = sum(label(ind(1:25)))/25;
precision_top30 = sum(label(ind(1:30)))/30;
precision_top35 = sum(label(ind(1:35)))/35;
precision_top40 = sum(label(ind(1:40)))/40;

result=[AUPR,AUC,f1,precision_top5,precision_top10,precision_top15,precision_top20,...
    precision_top25,precision_top30,precision_top35,precision_top40];

end

