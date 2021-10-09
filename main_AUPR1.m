clear
clc
load('CVa_result10.mat')

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
for i=1:size(option,1)
   flag = ismember(CVa_result10(:,1:4),option(i,:),'rows');  ind = find(flag==1);
   result_mean = [option(i,:),mean(CVa_result10(ind,6:end))];
   CVA_mean = [CVA_mean;result_mean];
end




idex = 4;
name = 'AUPR';


result_MNLMF = CVA_mean;

[~,ind] = max(result_MNLMF(:,idex+1));
best_option = result_MNLMF(ind,:)


%%%%计算范围

%%%%想把r-c进行联合
flag2 = ismember(result_MNLMF(:,[2,3]),best_option([2,3]),'rows');
ind2 = find(flag2==1);
r_c_result =  result_MNLMF(ind2,:);
[~,ind] = unique(r_c_result(:,1:4),'rows')
r_c_result = r_c_result(ind,:);
[min(r_c_result(:,idex+1)),max(r_c_result(:,idex+1))]



%%%% lata_ar联合
flag2 = ismember(result_MNLMF(:,[1,4]),best_option([1,4]),'rows');
ind2 = find(flag2==1);
lata_ar_result =  result_MNLMF(ind2,:);
[~,ind] = unique(lata_ar_result(:,1:4),'rows')
lata_ar_result = lata_ar_result(ind,:);
[min(lata_ar_result(:,idex+1)),max(lata_ar_result(:,idex+1))]


%%%%%%%%%%%%%%%%%%%%%%%%%%画图
set(gcf,'unit','centimeters','Position',[3,3,40,22]);          %%设置整张图片在电脑窗口的位置和大小
pos = [0.08,0.17,0.33,0.79;
    0.51,0.17,0.41,0.79];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% r-c
subplot('Position',pos(1,:))
hold on
ind50 = find(r_c_result(:,1)==50);
ind100 = find(r_c_result(:,1)==100);
h(1) = plot(r_c_result(ind50,4),r_c_result(ind50,idex+1),'-r ','LineWidth',2);
h(2) = plot(r_c_result(ind100,4),r_c_result(ind100,idex+1),'-b','LineWidth',2);
h(3) = plot(r_c_result(ind50,4),r_c_result(ind50,idex+1),'r . ','MarkerSize',18);
h(4) = plot(r_c_result(ind100,4),r_c_result(ind100,idex+1),'b o','MarkerSize',8);
bestx = [0,best_option(4),best_option(4)];
besty = [best_option(5),best_option(5),0];
% h(5) = plot(bestx,besty,'k--','LineWidth',0.5);
legend(h(1:2),{['\itr', '\rm = 50'],['\itr', '\rm = 100']},'Fontname','Times New Roman','FontSize',18)

xlim([0.5,9.5]);
ylim([0.94,0.952])

ykedu  = [0.94:0.002:0.952];
set(gca,'YTick',ykedu)
yy = cell(1,length(ykedu));
for i=1:length(ykedu)
    yy{i} = num2str(ykedu(i));
end
set(gca,'YTickLabel',yy,'Fontname','Times New Roman','FontSize',16)

set(gca,'XTick',1:2:9)
set(gca,'XTickLabel',{'1','3','5','7','9'},'Fontname','Times New Roman','FontSize',16)
% title(['Sensitivity analysis of parameters ','\itr',' \rmand ','\itc'],...
%     'Fontname','Times New Roman','FontSize',14)

xlabel('\itc','Fontname','Times New Roman','FontSize',20);
ylabel('AUPR','Fontname','Times New Roman','FontSize',18);
box on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% lata_ar
subplot('Position',pos(2,:))
hold on
aupr = reshape(lata_ar_result(:,5),5,5);
aupr = aupr - 0.8;
h = bar3(aupr,0.5);
for n=1:numel(h)
    cdata=get(h(n),'zdata');
    set(h(n),'cdata',cdata,'facecolor','interp')
end

set(gca,'XTick',1:5)
set(gca,'XTickLabel',{'10^{-2}','10^{-1}','10^0','10^1','10^2'},'Fontname','Times New Roman','FontSize',16)
set(gca,'YTick',1:5)
set(gca,'YTickLabel',{'10^{-2}','10^{-1}','10^0','10^1','10^2'},'Fontname','Times New Roman','FontSize',16)
xlim([0.2,5.7]);
ylim([0.2,5.7]);
zlim([0,0.2]);
zlabel_value = 0:0.05:0.2;
true_zlabel = cell(1,length(zlabel_value));
for i=1:length(zlabel_value)
    true_zlabel{i} = num2str(zlabel_value(i)+0.8);
end
set(gca,'ZTick',zlabel_value)
set(gca,'ZTickLabel',true_zlabel,'Fontname','Times New Roman','FontSize',16)

xlabel('\it\lambda','Fontname','Times New Roman','FontSize',20);
ylabel('\it\alpha','Fontname','Times New Roman','FontSize',20);
zlabel('AUPR','Fontname','Times New Roman','FontSize',18);
set(gca,'clipping','on');
view([-30.4 24.4]);
hco = colorbar('position',[0.94 0.2 0.02 0.7]);   %maxaupr  = 0.2231 minaupr = 0.05
ylabel_value = 0:0.05:0.2;
ylabel_name = cell(1,length(ylabel_value));
for i=1:length(ylabel_name)
    ylabel_name{i} = num2str(ylabel_value(i)+0.8);
end
set(hco,'YTick',ylabel_value);
set(hco,'YTickLabel',ylabel_name,'Fontname','Times New Roman','FontSize',16);

% title(['Sensitivity analysis of parameters ','\it\lambda',' \rmand ','\it\alpha'],...
%     'Fontname','Times New Roman','FontSize',14)
box off
grid on