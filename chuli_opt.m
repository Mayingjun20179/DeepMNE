%%%%�洢������miRNA�����ƾ�������Ӿ���
function L_D = chuli_opt(lnc_dis_inf2019)
%L_D���ж�Ӧ����lncRNA���ж�Ӧ����disease
%%%������
dis_sim =  DCA2_opt({lnc_dis_inf2019.dis_inf.tree_sim,lnc_dis_inf2019.dis_inf.GOMF_sim,lnc_dis_inf2019.dis_inf.gen_sim});
%%%����lncRNA
lnc_sim_d2 = lnc_dis_inf2019.lnc_inf.d2sim;
lnc_fea_exp = lnc_dis_inf2019.lnc_inf.exp_fea;   lnc_sim_exp = KSNS_opt(lnc_fea_exp);
lnc_fea_pc = lnc_dis_inf2019.lnc_inf.pc_fea;    lnc_sim_pc = KSNS_opt(lnc_fea_pc);
lnc_sim =  DCA2_opt({lnc_sim_exp,lnc_sim_pc});
%%%�������
L_D.interaction = lnc_dis_inf2019.interaction;
L_D.lnc_sim = lnc_sim;
L_D.dis_sim = dis_sim;

end



