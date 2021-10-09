%% this is the configuration file of stacked denoising autoencoder %%
function [sae0,opts0,nnsize0] = BuildNets0_opt(dim)
%1层自编码网络
%第一层节点数为样本的维度，第二层有200个节点
nnsize0 = [dim,floor(0.3*dim),floor(0.1*dim)];
len = length(nnsize0);

rand('state',0)
sae0 = saesetup(nnsize0);

for i = 1: len - 1   
    sae0.ae{i}.activation_function       = 'sigm';
    sae0.ae{i}.dropoutFraction           = 0;          %  Dropout fraction, only used for fine tuning
    sae0.ae{i}.momentum                  = 0;          %  Momentum
    sae0.ae{1}.scaling_learningRate      = 0.95;          %  Scaling factor for the learning rate (each epoch)
    sae0.weightPenaltyL2                  = 0.5;
    sae0.ae{i}.nonSparsityPenalty        = 0;          %  0 indicates Non sparsity penalty
    sae0.ae{i}.sparsityTarget            = 0.01;       %  Sparsity target
    sae0.ae{i}.inputZeroMaskedFraction   = 0;        %  Used for Denoising AutoEncoders
end

sae0.ae{1}.learningRate              = 0.25;
sae0.ae{2}.learningRate              = 0.25;
opts0.numepochs = 100;   %训练次数
opts0.batchsize = 50;    %
end