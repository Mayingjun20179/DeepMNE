function S_S = DCA2_opt(S)
%%����
%SΪ�ṹ���������������ƾ���
%reset_prob Ϊ��������
%dΪ��Ҫ�õ�������������ά��
%%���
%S_F  Ϊ��������
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

