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

Import all the codes into the matlab path, run mainZ, the results of DeepMNE in Table 1 and Table 2 can be automatically calculated.

GRGMF.zip: Reference from "Z. C. Zhang, X. F. Zhang, M. Wu, L. Ou-Yang, X. M. Zhao, and X. L. Li, “A graph regularized generalized matrix factorization model for predicting links in biomedical bipartite networks,” Bioinformatics, vol. 36, no. 11, pp. 3474-3481, Jun 1, 2020."

LDGRNMF.zip: Reference from "G. Li, J. Luo, C. Liang, Q. Xiao, P. Ding, and Y. Zhang, “Prediction of lncRNA-disease associations based on network consistency projection,” IEEE Access, vol. 7, pp. 58849-58856, 2019."

LDA-LNSUBRW.zip: Reference from "G. Xie, J. Jiang, and Y. Sun, “LDA-LNSUBRW: lncRNA-disease association prediction based on linear neighborhood similarity and unbalanced bi-random walk,” IEEE/ACM Trans Comput Biol Bioinform, vol. PP, Sep 1, 2020."

NCPLDA.zip: Reference from "G. Li, J. Luo, C. Liang, Q. Xiao, P. Ding, and Y. Zhang, “Prediction of lncRNA-disease associations based on network consistency projection,” IEEE Access, vol. 7, pp. 58849-58856, 2019."

