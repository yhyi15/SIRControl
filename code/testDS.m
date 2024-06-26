clear;
clc;

% read contact graph from file
% G = readGraph(filename);
%n = 200;

%parameters
%delta = 0.35;
%beta = 0.03;

rdEdges = readmatrix('data/randGraphTestModels.txt');

%create an arbitrary graph
%
n = max(max(rdEdges));
m = size(rdEdges, 1);
G = graph(rdEdges(:,1),rdEdges(:,2),ones(m,1),n);

GCC = G;
gccSize = numnodes(GCC);
maxd = max(degree(G));
disp(maxd);
beta = 0.25/(maxd+1);
disp(beta);
delta = 0.25;
%set initial conditions
s = 5;
%choose s seeds
%S = randsample(gccSize,s);
S = [16, 42, 6, 28, 3];
% initiate x and r
x0 = zeros(gccSize, 1);
r0 = zeros(gccSize, 1);
for i = 1: s
    x0(S(i)) = 1;
end
for i = 1:gccSize
    r0(i) = 0;
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
    if(max(x)<=10^-4)
        break;
    end
end
sigma = ones(1,nn)*(x+r-x0-r0);
originalRstReal = sigma;
disp("%%%%%%%DS real%%%%%%")
disp("original:")
disp(sigma);

A2=A;
A3=A;

%for e= 1:k
%    A(P(e,1),P(e,2))=0;
%    A(P(e,2),P(e,1))=0;
%end

%run dynamics to calculate sigma(P)
%rounds
%rounds = 10000;
%after deleting edges
greedyRstReal = zeros(k+1,1);
greedyRstReal(1) = originalRstReal;
for e = 1:k
    A(P(e,1),P(e,2))=0;
    A(P(e,2),P(e,1))=0;
    x = x0;
    r = r0;
    for i = 1:rounds
        xt = x+diag(ones(nn,1)-x-r)*B*A*x - D*x;
        r = r+ D*x;
        x=xt;
        if(max(x)<=10^-4)
            break;
        end
    end
    sigma = ones(1,nn)*(x+r-x0-r0);
    greedyRstReal(e+1)=sigma;
end
disp("greedy actual:")
disp(sigma);


%%%%%%%%%%%%%%%%%%%%%% actual %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epochs = 15000;
rounds =2000;
x=x0;
r=r0;
prob = B*A2;
%contagion = zeros(nn,nn);

Gcount = 0;
for i = 1:epochs
%    disp(i);
    x = x0;
    r = r0;
    for j = 1:rounds
        randMatrix = rand(nn);
        contagion = zeros(nn,nn);
        for ni = 1:nn
            for nj = 1:nn
                if ni~=nj && randMatrix(ni,nj)<=prob(ni,nj)
                    contagion(ni,nj)=1;
                end
            end
        end
        xt = contagion*x + x + r;
        for kit =1:nn
            if xt(kit)>1
                xt(kit) = 1;
            end
        end
        recRand = rand([nn,1]);
        recover = zeros(nn,1);
        for kit = 1:nn
            if recRand(kit)<= delta
                recover(kit) = 1;
            end
        end
        xt = xt - recover.*x -r;
        r = r + recover.*x;
        x=xt;
        if x == zeros(nn,1)
            break;
        end
    end
    Gcount = Gcount + ones(1,nn)*r-s;
end
MarkovOriginal = Gcount/epochs;

greedyRstMarkov = zeros(k+1,1);
greedyRstMarkov(1) = MarkovOriginal;

for e = 1:k
    A2(P(e,1),P(e,2))=0;
    A2(P(e,2),P(e,1))=0;
    prob = B*A2;
    Gcount = 0;
    for i = 1:epochs
        %    disp(i);
        x = x0;
        r = r0;
        for j = 1:rounds
            randMatrix = rand(nn);
            contagion = zeros(nn,nn);
            for ni = 1:nn
                for nj = 1:nn
                    if ni~=nj && randMatrix(ni,nj)<=prob(ni,nj)
                        contagion(ni,nj)=1;
                    end
                end
            end
            xt = contagion*x + x + r;
            for kit =1:nn
                if xt(kit)>1
                    xt(kit) = 1;
                end
            end
            recRand = rand([nn,1]);
            recover = zeros(nn,1);
            for kit = 1:nn
                if recRand(kit)<= delta
                    recover(kit) = 1;
                end
            end
            xt = xt - recover.*x -r;
            r = r + recover.*x;
            x=xt;
            if x == zeros(nn,1)
                break;
            end
        end
        Gcount = Gcount + ones(1,nn)*r-s;
    end
    MarkovRst = Gcount/epochs;
    greedyRstMarkov(e+1)=MarkovRst;
end

outFile0 = fopen('testModel/ds-info.txt','w');
outFile1 = fopen('testModel/ds-upper.txt','w');
outFile2 = fopen('testModel/ds-real.txt','w');
outFile3 = fopen('testModel/ds-markov.txt','w');

numRmv = (1:k+1)'-ones(k+1,1);
fprintf(outFile0, '%d\n%d\n%d\n%d\n%d\n%8f\n%8f\n',nn, m, qsize, k, s, beta, delta);
for i = 1:k+1
    fprintf(outFile1,'%5f\n', greedyRst(i));
    fprintf(outFile2,'%5f\n', greedyRstReal(i));
    fprintf(outFile3,'%5f\n', greedyRstMarkov(i));
end

fclose(outFile0);
fclose(outFile1);
fclose(outFile2);
fclose(outFile3);
