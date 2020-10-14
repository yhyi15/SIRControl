clear;
clc;

n = 3;
delta = 0.25;

A=zeros(n,n);
G1 = graph(A);
G2 = graph(A);
for i = 1:n-1
    for j = i+1:n
        edgeCoin1 = binornd(1,0.5);
        edgeCoin2 = binornd(1,0.5);
        if edgeCoin1 ==1
            G1 = addedge(G1, i,j, 1);
        end
        if edgeCoin2 ==1
            G2 = addedge(G2, i,j, 1);
        end
    end
end

A1 = adjacency(G1,'weighted');
A2 = adjacency(G2,'weighted');

dmax1 = max(sum(A1));
dmax2 = max(sum(A2));

dmax = max(dmax1,dmax2);

beta = delta/(dmax+1);

betaList = beta* ones(n,1);
deltaList = delta * ones(n,1);
B = diag(betaList);
D = diag(deltaList);

x0 = zeros(n, 1);
r0 = zeros(n, 1);

%set initial conditions
s = 1;
%choose s seeds
S = randsample(n,s);
% initiate x and r
for i = 1: s
    x0(i) = 1;
end
for i = 1:n
    r0(i) = 0;
end

X0 = diag(x0);
R0 = diag(r0);
I = eye(n);
M1 = I - D + (I-X0-R0)*B*A1;
M2 = I - D + (I-X0-R0)*B*A2;

sigma1 = ones(1,n)* (M1+D-I) * ((I-M1)\x0);
sigma2 = ones(1,n)* (M2+D-I) * ((I-M2)\x0);
sigma3 = ones(1,n)* ((M1+M2)/2+D-I) * ((I-(M1+M2)/2)\x0);

if (sigma1+sigma2)/2<sigma3-10^-6
    disp("nonconvex");
end
disp("end");