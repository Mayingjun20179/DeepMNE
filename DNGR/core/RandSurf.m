%% Randomly Surf %%与随机游走很相似
%for more details, pls see our paper

function M = RandSurf(A, max_step, alpha)
num_nodes = length(A);
A = ScaleSimMat(A);   %A的对角线为0，每一列的和为1

P0 = eye(num_nodes, num_nodes);
P = P0;
M = zeros(num_nodes, num_nodes);

for i = 1: max_step
    P = alpha*P*A + (1-alpha)*P0;
    M = M + P;
end

end