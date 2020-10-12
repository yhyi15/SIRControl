% optimizing deterministic spread model
% input: a contact graph
% algorithm: minimize upper bound 
%\hat{sigma}= 1^T *(M+D-I)(I-M)^(-1) * x(0)
% output: sigma =  |(x(t)+ r(t))-(x(0)+ r(0))|_1 and \hat{sigma}
clear;
clc;

% read contact graph from file
% G = readGraph(filename);
n = 200;

%paramaters
delta = 0.25;
beta = 0.1/log(n);

%create an arbitrary graph
%
G = graph;
for i = 1:n-1
    for j = i+1:n
        edgeCoin = binornd(1,2*log(n)/n);
        if edgeCoin ==1
            G = addedge(G, i,j, 0.2 + 0.8*rand());
        end
    end
end

% take the (presumably) largest connected component
[bins, binsizes] = conncomp(G);
gccSize = max(binsizes);
idx = binsizes(bins) == gccSize;
SG = subgraph(G, idx);
%reorder nodes in the gcc
order = 1:gccSize;
GCC = reordernodes(SG, order);


%set initial conditions
s = 30;
%choose s seeds
S = randsample(gccSize,s);
% initiate x and r
x0 = zeros(gccSize, 1);
r0 = zeros(gccSize, 1);
for i = 1: s
    x0(i) = rand()/10;
end
for i = 1:gccSize
    r0(i) = rand()/100;
end

%input to the algorithm: an edge set, an integer k
edgelist = table2array(GCC.Edges);
m = numedges(GCC);
nn = numnodes(GCC);
%randomly choose 1/5 edges as candidate set
q = floor(m/5);
qidx = randsample(m, q);
Q = edgelist(qidx, :);
k = min(q, 20);
%matrices
A = adjacency(G,'weighted');
betaList = beta* ones(nn,1);
deltaList = delta * ones(nn,1);
B = diag(betaList);
D = diag(deltaList);
X0 = diag(x0);
R0 = diag(r0);
I = eye(nn);
M = I - D + (I-X0-R0)*B*A;
% choose deleting set using the surrogate
P=zeros(k,3);
for i = 1:k
    disp(i)
    %iterate over all remaining candidates
    curr = ones(1, nn)* (M+D-I)*((I-M)\x0);
    for e = 1:q
        Ainc = zeros(nn,nn);
        Ainc(Q(e,1), Q(e,2))= -A(Q(e,1), Q(e,2));
        Ainc(Q(e,2), Q(e,1))= -A(Q(e,2), Q(e,1));
        Mt = M + (I-X0-R0)*B*Ainc;
        sigmahat = ones(1, nn)* (Mt+D-I)*((I-Mt)\x0);
        if sigmahat<curr
            curr = sigmahat;
            choice = e;
            Mc = Mt;
        end
    end
    disp(sigmahat);
    P(i,:) = Q(choice,:);
    Q(choice,:)=[];
    q=q-1;
    M = Mc;
end
sigmahat = ones(1,nn)* (M-D-I) * ((I-M)\x0);


%run dynamics to calculate sigma(P)
%rounds
rounds = 10000;
%before deleting edges
x = x0;
r = r0;
for i = 1:rounds
    x = x+diag(ones(nn,1)-x-r)*B*A*x - D*x;
    r = r+ D*x;
end
sigma = ones(1,nn)*(x+r-x0-r0);
disp(sigma);

for e= 1:k
    A(P(e,1),P(e,2))=0;
    A(P(e,2),P(e,1))=0;
end

%run dynamics to calculate sigma(P)
%rounds
%rounds = 10000;
%after deleting edges
x = x0;
r = r0;
for i = 1:rounds
    x = x+diag(ones(nn,1)-x-r)*B*A*x - D*x;
    r = r+ D*x;
end
sigma = ones(1,nn)*(x+r-x0-r0);
disp(sigma);

