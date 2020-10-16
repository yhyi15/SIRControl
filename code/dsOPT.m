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
            G = addedge(G, i,j, 1);
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
    x0(i) = rand()/2;
end
for i = 1:gccSize
    r0(i) = rand()/20;
end

%input to the algorithm: an edge set, an integer k
edgelist = table2array(GCC.Edges);
m = numedges(GCC);
nn = numnodes(GCC);
%randomly choose 1/5 edges as candidate set
constraint = 2;
qsize=floor(m/2);
q=qsize;
qidx = randsample(m, q);
Q = edgelist(qidx, :);
k = floor(qsize/3);
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

disp("%%%%%%%upper bound%%%%%%")
disp("original:")
sigmahat = ones(1,nn)* (M+D-I) * ((I-M)\x0);
disp(sigmahat);

% choose deleting set using the surrogate
P=zeros(k,3);
%disp("%%%%%%%%%%%%%%%%%%%%%%%%%")
disp("greedy:")
for i = 1:k
    %disp(i)
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
    %disp(sigmahat);
    P(i,:) = Q(choice,:);
    %disp(choice)
    Q(choice,:)=[];
    q=q-1;
    M = Mc;
end
sigmahat = ones(1,nn)* (M+D-I) * ((I-M)\x0);
disp(sigmahat);

%to compare with random choices
%disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
disp("random choices:");
q2=qsize;
M = I - D + (I-X0-R0)*B*A;
P2=zeros(k,3);
Q2=edgelist(qidx, :);
for i = 1:k
    %disp(i)
    %iterate over all remaining candidates
    %curr = ones(1, nn)* (M+D-I)*((I-M)\x0);
    e =randi(q2,1);
    %disp(e)
    
    Ainc = zeros(nn,nn);
    Ainc(Q2(e,1), Q2(e,2))= -A(Q2(e,1), Q2(e,2));
    Ainc(Q2(e,2), Q2(e,1))= -A(Q2(e,2), Q2(e,1));
    M = M + (I-X0-R0)*B*Ainc;
    %sigmahat = ones(1, nn)* (Mt+D-I)*((I-Mt)\x0);
    
    %disp(sigmahat);
    P2(i,:) = Q2(e,:);
    Q2(e,:)=[];
    q2=q2-1;
end
sigmahat = ones(1,nn)* (M+D-I) * ((I-M)\x0);
disp(sigmahat);



%%%% to reduce maximum degree by removing arbitrary edges incident to the
%%%% maximum degree node
disp("removing an edge incident to maximum degree nodes")
q3=qsize;
M = I - D + (I-X0-R0)*B*A;
P3=zeros(k,3);
Q3=edgelist(qidx, :);
GCCCopy=GCC;
for i = 1:k
    %disp(i)
    %iterate over all remaining candidates
    %curr = ones(1, nn)* (M+D-I)*((I-M)\x0);
    dmax = 0;
    for e = 1:q
        dupdate = max(degree(GCCCopy, Q3(e, 1)), degree(GCCCopy, Q3(e, 2)));
        %dupdate = max(degree(GCCCopy, Q3(e, 1)), degree(GCCCopy, Q3(e, 2)));
        if dupdate>dmax
            choice = e;
            dmax = dupdate;
        end
    end
    
    Ainc = zeros(nn,nn);
    Ainc(Q3(e,1), Q3(e,2))= -A(Q3(e,1), Q3(e,2));
    Ainc(Q3(e,2), Q3(e,1))= -A(Q3(e,2), Q3(e,1));
    M = M + (I-X0-R0)*B*Ainc;
    
    %disp(sigmahat);
    P3(i,:) = Q3(choice,:);
    %disp(choice)
    Q3(choice,:)=[];
    GCCCopy = rmedge(GCCCopy, Q3(choice,1), Q3(choice,2));
    q3=q3-1;
end
sigmahat = ones(1,nn)* (M+D-I) * ((I-M)\x0);
disp(sigmahat);


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
disp("%%%%%%%DS real%%%%%%")
disp("original:")
disp(sigma);

A2=A;
A3=A;

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
disp("greedy actual:")
disp(sigma);

for e= 1:k
    A2(P2(e,1),P2(e,2))=0;
    A2(P2(e,2),P2(e,1))=0;
end

%run dynamics to calculate sigma(P)
%rounds
%rounds = 10000;
%after deleting edges
x = x0;
r = r0;
for i = 1:rounds
    x = x+diag(ones(nn,1)-x-r)*B*A2*x - D*x;
    r = r+ D*x;
end
sigma = ones(1,nn)*(x+r-x0-r0);
disp("random actual:")
disp(sigma);

for e= 1:k
    A3(P3(e,1),P3(e,2))=0;
    A3(P3(e,2),P3(e,1))=0;
end
x = x0;
r = r0;
for i = 1:rounds
    x = x+diag(ones(nn,1)-x-r)*B*A3*x - D*x;
    r = r+ D*x;
end
sigma = ones(1,nn)*(x+r-x0-r0);
disp("degree reducing actural:")
disp(sigma);
