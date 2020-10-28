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
A2 = A;
beta = beta;
betaList = beta* ones(nn,1);
deltaList = delta * ones(nn,1);
B = diag(betaList);
D = diag(deltaList);
X0 = diag(x0);
R0 = diag(r0);
I = eye(nn);
M = I - D + (I-X0-R0)*B*A;


%%%%%%%%%%%%%%%%%%%%%% ic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps = 0.2;
rounds = 5*nn*log(nn)/(eps^2);

count = 0;

disp(rounds);


P = zeros(k,3);
Ggreedy = graph(A);

A = 1/(1-(1-delta)*(1-beta))*beta*A;


disp("%%%%%%%%%%%% Original %%%%%%%%%%%%%%%%%")
cur_count = 0;
for round = 1:10*rounds
    terminal = randsample(nn,1);
    
    stack = zeros(nn, 1);
    stack(1) = terminal;
    top = 1;
    activated = zeros(nn,1);
    
    while top>0
        %pop
        cur = stack(top);
        top = top -1;
        if activated(cur)==0
            nb = neighbors(Ggreedy, cur);
            for iter = 1:size(nb,1)
                %push
                t = nb(iter);
                contCoin = rand();
                if (contCoin<=A(cur,t))
                    %push
                    top = top+1;
                    stack(top) = t;
                end
            end
            activated(cur)=1;
        end
    end
    for seed_iter= 1:s
        if activated(S(seed_iter))==1
            cur_count=cur_count+1;
            break;
        end
    end
end
infect = cur_count*nn/10/rounds;
disp("before removal:")
%disp(cur_count);
disp(infect-s);

originalRst = infect -s;

greedyRst = zeros(k+1,1);
greedyRst(1) = originalRst;

disp("%%%%%%%%%%%%%%% Greedy %%%%%%%%%%%%%%%%%%%%%");

disp(strcat("k: ", num2str(k)));
for i = 1:k
    disp(strcat("deleting the ",num2str(i),"-th edge"));
    cur_count = 0;
    count = zeros(q,1);
    for round = 1:rounds
        succ = 0;
        contEdge = zeros(nn,2);
        contNumEdges = 0;
        %if(mod(round,100)==0)
        %    disp(round/rounds);
        %end
        terminal = randsample(nn,1);
        
        stack = zeros(nn, 1);
        stack(1) = terminal;
        top = 1;
        activated = zeros(nn,1);
        
        while top>0
            %pop
            cur = stack(top);
            top = top -1;
            if activated(cur)==0
                nb = neighbors(Ggreedy, cur);
                for iter = 1:size(nb,1)
                    %push
                    t = nb(iter);
                    contCoin = rand();
                    if (contCoin<=A(cur,t))
                        %push
                        top = top+1;
                        stack(top) = t;
                        contNumEdges = contNumEdges + 1;
                        contEdge(contNumEdges,:) = [cur, t];
                    end
                end
                activated(cur)=1;
            end
        end
        for iter_c = 1:s
            if activated(S(iter_c))==1
                succ = 1;
                cur_count=cur_count+1;
                break;
            end
        end
        
        if succ ==0
            continue;
        end

        for e =1:q
            sn = Q(e,1);
            tn = Q(e,2);
            in_h = findedgeCont(contEdge, contNumEdges,sn,tn);
            if in_h ==0
                count(e) = count(e)+succ;
                continue;
            end
            stack_h = zeros(contNumEdges, 1);
            stack_h(1)= terminal;
            top_h = 1;
            activated_h = zeros(nn,1);
            while top_h>0
                cur_h = stack_h(top_h);
                top_h = top_h -1;
                if activated_h(cur_h)==0
                    nb_h = neighborsCont(contEdge, contNumEdges, cur_h);
                    %disp(size(nb_h));
                    for iter_h = 1:size(nb_h,1)
                        %push
                        t_h = nb_h(iter_h);
                        if (~((cur_h==sn && t_h == tn) || (cur_h==tn && t_h == sn)))
                            %push
                            top_h = top_h+1;
                            stack_h(top_h) = t_h;
                        end
                    end
                    activated_h(cur_h)=1;
                end
            end
            for iter_h_ss = 1:s
                if activated_h(S(iter_h_ss))==1
                    count(e)=count(e)+1;
                    break;
                end
            end
        end
    end
    %disp("%%%%%%%%%%%%%%%")
    %disp(cur_count);
    [minv, idx] = min(count);
    %disp(minv);
    P(i,:) = Q(idx,:);
    Ggreedy = rmedge(Ggreedy,Q(idx,1),Q(idx,2));
    Q(idx,:)=[];
    q=q-1;
    
    greedyRst(i+1) = icEval(Ggreedy, nn, rounds,S ,s, A)-s;
end



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

outFile0 = fopen('testModel/ic-info.txt','w');
outFile1 = fopen('testModel/ic-upper.txt','w');
outFile2 = fopen('testModel/ic-markov.txt','w');

numRmv = (1:k+1)'-ones(k+1,1);
fprintf(outFile0, '%d\n%d\n%d\n%d\n%d\n%8f\n%8f\n',nn, m, qsize, k, s, beta, delta);
for i = 1:k+1
    fprintf(outFile1,'%5f\n', greedyRst(i));
    fprintf(outFile2,'%5f\n', greedyRstMarkov(i));
end

fclose(outFile0);
fclose(outFile1);
fclose(outFile2);
