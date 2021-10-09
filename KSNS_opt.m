%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   KSNS
function S = KSNS_opt(X,sim)
if nargin==1
    distance = pdist(X); 
    distance = squareform(distance);  %º∆À„≈∑ Ωæ‡¿Î
    sim = max(distance(:))-distance;
end
neighbor_num = floor(0.3*size(sim,1));
% X = gpuArray(X);
% sim = gpuArray(sim);
feature_matrix = X;
nearst_neighbor_matrix = KSNS_neighbors(sim,neighbor_num);
S = jisuanW_opt(feature_matrix,nearst_neighbor_matrix);
% class(S)
S =  gather(S);
S = double(S);
end

function nearst_neighbor_matrix=KSNS_neighbors(sim,neighbor_num)
%%%%nearst_neighbor_matrix£∫  represent Neighbor matrix
  N = size(sim,1);
  D = sim-diag(inf*ones(1,N));
  [~, si]=sort(D,2,'descend');
  nearst_neighbor_matrix=gpuArray.zeros(N,N);
  index=si(:,1:neighbor_num);
  for i=1:N
      nearst_neighbor_matrix(i,index(i,:))=1;     
  end
end

%The iterative process of this algorithm
function [W,objv] = jisuanW_opt(feature_matrix,nearst_neighbor_matrix)
lata1 = 4;  lata2 = 1;
X=feature_matrix';  % each column of X represents a sample, and each behavior is an indicator
[~,N] = size(X);    % N represents the number of samples
C = nearst_neighbor_matrix';
rand('state',1);
W = rand(N,N);
W = W- diag(diag(W));
W = W./repmat(sum(W),N,1);
G  = jisuan_Kel(X);
G(isnan(G))=0;
G = G/max(G(:));

W = gpuArray(W);
G = gpuArray(G);


WC1 = W'*G*W-2*W*G+G;
WC = sum(diag(WC1))/2;
wucha = WC + norm(W.*(1-C),'fro')^2*lata1/2 +  norm(W,'fro')^2*lata2/2;
objv = wucha;
jingdu = 0.0001;
error = jingdu*(1+lata1+lata2);   %Iteration error threshold
we = 1;      %Initial value of error
gen=1;
while  gen<100 && we>error
    %gen
    FZ = G+lata1*C.*W;
    FM = G*W+lata1*W+lata2*W;        
    FM(FM==0) = eps;  W = FZ./FM.*W;  
%     W = FZ./(FM+eps).*W;  
    WC1 = W'*G*W-2*W*G+G;
    WC = sum(diag(WC1))/2;
    objv1 = WC + norm(W.*(1-C),'fro')^2*lata1/2 +  norm(W,'fro')^2*lata2/2;    
    we = abs(objv1-objv(end));
    objv = [objv,objv1];
    gen = gen+1;
end
W=W'; 
W = matrix_normalize(W);
end

function W = matrix_normalize(W)
K = 10;
W(isnan(W))=0;
W(1:(size(W,1)+1):end)=0;
for round=1:K
    SW = sum(W,2);
    ind = find(SW>0);
    SW(ind) = 1./sqrt(SW(ind));
    D1 = diag(SW);
    W=D1*W*D1;
end
W = gather(W);
end

function K  =jisuan_Kel(X)
%X Columns represent samples, and rows represent features
X = X';
sA = (sum(X.^2, 2));
sB = (sum(X.^2, 2));
K = exp(bsxfun(@minus,bsxfun(@minus,2*X*X', sA), sB')/mean(sA));
end
