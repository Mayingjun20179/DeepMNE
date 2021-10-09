# DeepMNE
The Code and data of "Deep Multi-network Embedding for lncRNA-Disease Association prediction"

Data:

lnc_dis_inf2019. mat: Total data file, Matlab mat format file

Disease:

lnc_dis_inf2019.dis_inf:  Disease  information

lnc_dis_inf2019.dis_inf.nameID: The name and ID of disease

lnc_dis_inf2019.dis_inf.gen_sim: Gene function similarity of disease

lnc_dis_inf2019.dis_inf.GOMF_sim: GO functional similarity of diseases

lnc_dis_inf2019.dis_inf.tree_sim: Semantic similarity of diseases

LncRNA:

lnc_dis_inf2019.lnc_inf: LncRNA information

lnc_dis_inf2019.lnc_inf.exp_fea: The expression profile feature of lncRNA

lnc_dis_inf2019.lnc_inf.pc_fea: The PC-PseDNC features of the lncRNA 

Code：
xiuzheng_Y:Calculation of weighted K nearest Neighbor Profiles (WKNNP)

KSNS_opt: Compute kernel neighborhood similarity.

DNGR： Reference from “S. Cao, W. Lu, and Q. Xu, “Deep neural networks for learning graph representations,” in Proceedings of the Thirtieth AAAI Conference on Artificial Intelligence (AAAI-16), 2016, pp. 1145-1152.”

mainZ：The main function of DeepMNE

Import all the codes into the matlab path, run mainZ, the results of DeepMNE in Table 1 and Table 2 can be automatically calculated

