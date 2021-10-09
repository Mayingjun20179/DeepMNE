function main_cv(cv_flag,DATA)

%���ɲ���
%option(1)  �ӿռ�ά�ȣ�
%option(2) �������ճ������򻯲���
%option(3)  �����Ǳ�ڿռ��������򻯲���
%option(4)  ��Ҫ��ˮƽ����

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


%%���ɽ���ʵ��Ĳ��Լ�
cv = 5;  %cv -�۽���ʵ��
%NRLMF�㷨������
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



%%%%�����ڲ�ͬ�������µ�Ԥ����
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


%%%%%%%%%%%%%%%%%%%%%%���ɽ���ʵ������ݼ�
%%%%���ɲ�ͬ����ʱ���жϾ��󣨲��Լ�������Ϊ0������Ϊ1��
%%%%���Լ������������Լ��ı�ǩ
%intMat ��������
%seed   ��ʾ���������������
%cv_flag     ��ʾѡ����һ�ֽ���ʵ��
%cv_flag =1  ��ʾ�����н����н���ʵ��
%cv_flag = 2 ��ʾ���н���ʵ��
%cv_flag = 3 ��ʾ���н���ʵ��
%cv    ��ʾ����ʵ��Ĵ�����5�ۣ�10�ۣ�
function cv_data = cross_validation(DATA,seeds,cv_flag,cv)

intMat = DATA.interaction;
Ns = length(seeds);
cv_data = cell(Ns,cv);  %cv_data{i,j}��ʾ����Ϊseeds(i)ʱ����j������ʵ���Ӧ���б����W
%���Ե��λ�� test_data�����Ե�ı�ǩ test_label
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
        for j = 1:cv   %����ÿ�ν���ʵ���Ӧ�Ĳ��Լ�
            if j < cv
                ii = index((j-1)*step+1:j*step);
            else
                ii = index((j-1)*step+1:end);
            end
            if cv_flag == 1    %�����н�������
                test_data = [Nx1(ii),Ny1(ii)];   %����1/5��1��λ��
                [indx0,indy0] = find(intMat==0); %����ȫ����0��λ��
                %%%%%�޸�Ϊ1��0�ĸ���һ�����
                N1 = size(test_data,1);
                N0 = length(indx0);  indx00 = randperm(N0); 
                test_data = [test_data;[indx0(indx00(1:N1)),indy0(indx00(1:N1))]];   %����1/5��1��λ��
%                 test_data = [test_data;[indx0,indy0]];   %����1/5��1��λ��
            elseif cv_flag == 2   %�������н���
                test_data=[];
                for k=1:length(ii)
                    test_data = [test_data;[ii(k)*ones(Nv,1),[1:Nv]']];
                end
            elseif cv_flag==3     %�������н���
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


