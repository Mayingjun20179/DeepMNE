function rep = MDA_train(xx)
%%%xx��ʾ�������繹�ɵĽṹ�����
%d��ʾ�ںϺ������ά��
% d = floor(0.1 * size(xx{1},2));
d = 20;
%%%%��һ�㣬�����е�����������б���
Nn = length(xx);  %��������ĸ���
x = [];   %���������
for i=1:Nn
    dim = size(xx{i},2);
    [sae0,opts0] = BuildNets0_opt(dim);
    sae0.ae{i} = nntrain(sae0.ae{1}, xx{i}, xx{i}, opts0);
    t = nnff(sae0.ae{i}, xx{i}, xx{i});
    xx0 = t.a{2};
    x = [x,xx0(:,2:end)];
end


% �ڶ�������һ���Ա���ṹ��input x����opt����
dim = size(x,2);
[sae,opts,nnsize] = BuildNets_opt(dim,d);
sae = saetrain(sae, x, opts);
rep = GenRep(x, sae, nnsize);  
end
