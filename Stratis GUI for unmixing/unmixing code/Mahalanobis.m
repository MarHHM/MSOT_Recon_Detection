function [d] = Mahalanobis(mixed)

%format of mixed is lambda * x * y;   
% switch size(size(mixed))
%     case 3

S = size(mixed); 
M = S(1); N = S(2);
miu = zeros(M,1);
T = zeros(M,N);
G = zeros(M);

%remmean
for i = 1:M
    miu(i,1) = 1/N *sum(mixed(i,:));
    T(i,:) = mixed(i,:)-miu(i,1); 
end 

%corrected
for j = 1:N
    G = T(:,j) * T(:,j)' + G;
end 

G = G / (N-1);% figure,imagesc(G);
G_inv = inv(G);

d = zeros(N,1);
for j = 1:N
    d(j) = (mixed(:,j)- miu(:,1))'* G_inv * (mixed(:,j)- miu(:,1));
end

end 