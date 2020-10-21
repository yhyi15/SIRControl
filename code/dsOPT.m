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

%parameters
delta = 0.25;
beta = 0.1/log(n);

%create an arbitrary graph
%
G = graph;
for i = 1:n-1
    for j = i+1:n
        edgeCoin = binornd(1,2*log(n)/n);
        if edgeCoin ==1
            G = addedge(G, i,j, 0.2);%/(2.5*log(n)));
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
s = 20;
%choose s seeds
S = randsample(gccSize,s);
% initiate x and r
x0 = zeros(gccSize, 1);
r0 = zeros(gccSize, 1);
for i = 1: s
    x0(S(i)) = 0.8 +0.2*rand();
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
A = adjacency(GCC,'weighted');
betaList = 0.5*beta* ones(nn,1) + 1*beta*rand([nn,1]);
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
originalRst = sigmahat;
disp(sigmahat);

% choose deleting set using the surrogate
P=zeros(k,3);
greedyRst = zeros(k+1,1);
greedyRst(1) = originalRst;
%disp("%%%%%%%%%%%%%%%%%%%%%%%%%")
disp("greedy:")
initInv = inv(I-M);
curInv = initInv;
for i = 1:k
    %disp(i)
    %iterate over all remaining candidates
    curr = ones(1, nn)* (M+D-I)*(curInv*x0);
    minval=curr;
    % = 0;
    curLft = ones(1,nn) * (M+D-I);
    curRt = curInv * x0;
    result = zeros(q, 1);
    for e = 1:q
        %Ainc = zeros(nn,nn);
        %Ainc(Q(e,1), Q(e,2))= -A(Q(e,1), Q(e,2));
        %Ainc(Q(e,2), Q(e,1))= -A(Q(e,2), Q(e,1));
        %Mt = M + (I-X0-R0)*B*Ainc;
        %sigmahat = ones(1, nn)* (Mt+D-I)*((I-Mt)\x0);
        %disp(sigmahat);
        ei = Q(e,1); ej = Q(e,2);
        Ri = curInv(ei,:);
        Rj = curInv(ej,:);
        Ci = curInv(:,ei);
        Cj = curInv(:,ej);
        Aij = B(ei,ei)*A(ei,ej);
        Aji = B(ej,ej)*A(ej,ei);
        lft = curLft;
        lft(ej)= lft(ej) - (1-x0(ei)-r0(ei))*Aij;
        lft(ei)= lft(ei) - (1-x0(ej)-r0(ej))*Aji;
        infectNum = numUpdate(nn, Ri,Rj,Ci,Cj, x0, r0, lft, curRt, ei, ej, Aij, Aji);
        result(e)=infectNum;
        %disp(infectNum);
        %if infectNum<minval
        %    minval = infectNum;
        %    choice = e;
        %    tempLft = lft;
        %    flg =1;
        %end
    end
    for e = 1:q
        if result(e)<minval
            minval = result(e);
            choice = e;
        end
    end
    %disp(choice);
    %disp(sigmahat);
    greedyRst(i+1) = minval;
    P(i,:) = Q(choice,:);
    M(Q(choice,1),Q(choice,2)) = 0;
    M(Q(choice,2),Q(choice,1)) = 0;
    Aij = B(Q(choice,1),Q(choice,1))*A(Q(choice,1),Q(choice,2));
    Aji = B(Q(choice,2),Q(choice,2))*A(Q(choice,2),Q(choice,1));
    curLft(Q(choice,2))= curLft(Q(choice,2)) - (1-x0(Q(choice,1))-r0(Q(choice,1)))*Aij;
    curLft(Q(choice,1))= curLft(Q(choice,1)) - (1-x0(Q(choice,2))-r0(Q(choice,2)))*Aji;
    curInv = invUpdate(curInv, x0, r0, Q(choice,1),Q(choice,2),Aij, Aji);
    %disp(choice)
    Q(choice,:)=[];
    q=q-1;
end
sigmahat = ones(1,nn)* (M+D-I) * ((I-M)\x0);
disp(sigmahat);

%to compare with random choices
%disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
disp("random choices:");
randRst = zeros(k+1,1);
randRst(1) = originalRst;
q2=qsize;
M = I - D + (I-X0-R0)*B*A;
P2=zeros(k,3);
Q2=edgelist(qidx, :);
curInv = initInv;
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
    %%sigmahat = ones(1, nn)* (Mt+D-I)*((I-Mt)\x0);
    choice = e;
    P2(i,:) = Q2(choice,:);
    M(Q2(choice,1),Q2(choice,2)) = 0;
    M(Q2(choice,2),Q2(choice,1)) = 0;
    Aij = B(Q2(choice,1),Q2(choice,1))*A(Q2(choice,1),Q2(choice,2));
    Aji = B(Q2(choice,2),Q2(choice,2))*A(Q2(choice,2),Q2(choice,1));
    curLft(Q2(choice,2))= curLft(Q2(choice,2)) - (1-x0(Q2(choice,1))-r0(Q2(choice,1)))*Aij;
    curLft(Q2(choice,1))= curLft(Q2(choice,1)) - (1-x0(Q2(choice,2))-r0(Q2(choice,2)))*Aji;
    curInv = invUpdate(curInv, x0, r0, Q2(choice,1),Q2(choice,2),Aij, Aji);
    
    sigmahat = (ones(1,nn)* (M+D-I)) * (curInv * x0);
    randRst(i+1)=sigmahat;
    
    %disp(sigmahat);
    P2(i,:) = Q2(e,:);
    Q2(e,:)=[];
    q2=q2-1;
end
sigmahat = (ones(1,nn)* (M+D-I)) * (curInv * x0);
disp(sigmahat);



%%%% to reduce maximum degree by removing arbitrary edges incident to the
%%%% maximum degree node
disp("removing an edge incident to maximum degree nodes")
maxDRst = zeros(k+1,1);
maxDRst(1) = originalRst;
q3=qsize;
M = I - D + (I-X0-R0)*B*A;
P3=zeros(k,3);
Q3=edgelist(qidx, :);
GCCCopy=GCC;
curInv = initInv;
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
    Ainc(Q3(choice,1), Q3(choice,2))= -A(Q3(choice,1), Q3(choice,2));
    Ainc(Q3(choice,2), Q3(choice,1))= -A(Q3(choice,2), Q3(choice,1));
    M = M + (I-X0-R0)*B*Ainc;
    
    choice = e;
    P3(i,:) = Q3(choice,:);
    M(Q3(choice,1),Q3(choice,2)) = 0;
    M(Q3(choice,2),Q3(choice,1)) = 0;
    Aij = B(Q3(choice,1),Q3(choice,1))*A(Q3(choice,1),Q3(choice,2));
    Aji = B(Q3(choice,2),Q3(choice,2))*A(Q3(choice,2),Q3(choice,1));
    curLft(Q3(choice,2))= curLft(Q3(choice,2)) - (1-x0(Q3(choice,1))-r0(Q3(choice,1)))*Aij;
    curLft(Q3(choice,1))= curLft(Q3(choice,1)) - (1-x0(Q3(choice,2))-r0(Q3(choice,2)))*Aji;
    curInv = invUpdate(curInv, x0, r0, Q3(choice,1),Q3(choice,2),Aij, Aji);
    
    sigmahat = (ones(1,nn)* (M+D-I)) * (curInv * x0);
    maxDRst(i+1)=sigmahat;
    
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
rounds = 1000;
%before deleting edges
x = x0;
r = r0;
for i = 1:rounds
    xt = x+diag(ones(nn,1)-x-r)*B*A*x - D*x;
    r = r+ D*x;
    x=xt;
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
    xt = x+diag(ones(nn,1)-x-r)*B*A*x - D*x;
    r = r+ D*x;
    x=xt;
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
    xt = x+diag(ones(nn,1)-x-r)*B*A2*x - D*x;
    r = r+ D*x;
    x=xt;
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
    xt = x+diag(ones(nn,1)-x-r)*B*A3*x - D*x;
    r = r+ D*x;
    x=xt;
end
sigma = ones(1,nn)*(x+r-x0-r0);
disp("degree reducing actural:")
disp(sigma);
