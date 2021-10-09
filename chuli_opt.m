%%%%存储疾病和miRNA的相似矩阵和连接矩阵
function L_D = chuli_opt(lnc_dis_inf2019)
%L_D中行对应的是lncRNA，列对应的是disease
%%%处理疾病
dis_sim =  DCA2_opt({lnc_dis_inf2019.dis_inf.tree_sim,lnc_dis_inf2019.dis_inf.GOMF_sim,lnc_dis_inf2019.dis_inf.gen_sim});
%%%处理lncRNA
lnc_sim_d2 = lnc_dis_inf2019.lnc_inf.d2sim;
lnc_fea_exp = lnc_dis_inf2019.lnc_inf.exp_fea;   lnc_sim_exp = KSNS_opt(lnc_fea_exp);
lnc_fea_pc = lnc_dis_inf2019.lnc_inf.pc_fea;    lnc_sim_pc = KSNS_opt(lnc_fea_pc);
lnc_sim =  DCA2_opt({lnc_sim_exp,lnc_sim_pc});
%%%处理关联
L_D.interaction = lnc_dis_inf2019.interaction;
L_D.lnc_sim = lnc_sim;
L_D.dis_sim = dis_sim;

end



