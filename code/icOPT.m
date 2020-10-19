% optimizing deterministic spread model
% input: a contact graph
% algorithm: minimize upper bound 
%\hat{sigma}= 1^T *(M+D-I)(I-M)^(-1) * x(0)
% output: sigma =  |(x(t)+ r(t))-(x(0)+ r(0))|_1 and \hat{sigma}
clear;
clc;

% read contact graph from file
% G = readGraph(filename);
n = 100;

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
            G = addedge(G, i,j, (0.1 + 0.9*rand())/(2.5*log(n)));
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
betaList = 0.2*beta* ones(nn,1) + 0.8*beta*rand([nn,1]);
deltaList = delta * ones(nn,1);
B = diag(betaList);
D = diag(deltaList);
X0 = diag(x0);
R0 = diag(r0);
I = eye(nn);
M = I - D + (I-X0-R0)*B*A;


%calculating original expected infections
%randomly pick a node from the network
eps = 0.1;
rounds = nn*log(nn)/(eps^2);

count = 0;

disp(rounds);


P = zeros(k,3);
Ggreedy = GCC;

cur_count = 0;
for round = 1:rounds
    Ghat = graph;
    if(mod(round,5000)==0)
        disp(round/rounds);
    end
    terminal = randsample(nn,1);
    
    stack = zeros(nn, 1);
    stack(1) = terminal;
    top = 1;
    visited = zeros(nn,1);
    while top>0
        %pop
        cur = stack(top);
        top = top -1;
        flag = 0;
        for iter_s = 1:s
            if cur == S(iter_s)
                flag = 1;
                break;
            end
        end
        if flag ==1
            cur_count = cur_count+1;
            break;
        end
        nb = neighbors(Ggreedy, cur);
        for iter = 1:size(nb,1)
            %push
            t = nb(iter);
            contCoin = rand();
            if (visited(t)==0 && contCoin<=A(cur,t))
                %push
                top = top+1;
                stack(top) = t;
                Ghat = addedge(Ghat, cur,t, 1);
                visited(t) = 1;
            end
        end
    end
end
infect = cur_count*nn/rounds;
disp("%%%%%%%%%%%%%% before removal %%%%%%%%%%%%%")
disp(infect);


for i = 1:k
    cur_count = 0;
    count = zeros(q,1);
    for round = 1:rounds
        Ghat = graph;
        if(mod(round,5000)==0)
            disp(round/rounds);
        end
        terminal = randsample(nn,1);
        
        stack = zeros(nn, 1);
        stack(1) = terminal;
        top = 1;
        visited = zeros(nn,1);
        while top>0
            %pop
            cur = stack(top);
            top = top -1;
            flag = 0;
            for iter_s = 1:s
                if cur == S(iter_s)
                    flag = 1;
                    break;
                end
            end
            if flag ==1
                cur_count = cur_count+1;
                break;
            end
            nb = neighbors(Ggreedy, cur);
            for iter = 1:size(nb,1)
                %push
                t = nb(iter);
                contCoin = rand();
                if (visited(t)==0 && contCoin<=A(cur,t))
                    %push
                    top = top+1;
                    stack(top) = t;
                    Ghat = addedge(Ghat, cur,t, 1);
                    visited(t) = 1;
                end
            end
        end
        for e =1:q
            s = Q(e,1);
            t = Q(e,2);
            exist = findedge(Ghat,s,t);
            if exist>0
                Ghat = rmedge(Ghat,s,t);
                components=conncomp(Ghat);
                for iter_ss=1:s
                    if components(terminal)==components(S(iter_ss))
                        count(e)= count(e)+1;
                    end
                end
                Ghat = addedge(Ghat,s,t,1);
            end
        end
    end
    [minv, idx] = min(count);
    P(k,:) = Q(idx,:);
    Q(idx,:)=[];
    q=q-1;
    Ggreedy = rmedge(Ggreedy,Q(idx,1),Q(idx,2));
end

cur_count = 0;
for round = 1:rounds
    Ghat = graph;
    if(mod(round,5000)==0)
        disp(round/rounds);
    end
    terminal = randsample(nn,1);
    
    stack = zeros(nn, 1);
    stack(1) = terminal;
    top = 1;
    visited = zeros(nn,1);
    while top>0
        %pop
        cur = stack(top);
        top = top -1;
        flag = 0;
        for iter_s = 1:s
            if cur == S(iter_s)
                flag = 1;
                break;
            end
        end
        if flag ==1
            cur_count = cur_count+1;
            break;
        end
        nb = neighbors(Ggreedy, cur);
        for iter = 1:size(nb,1)
            %push
            t = nb(iter);
            contCoin = rand();
            if (visited(t)==0 && contCoin<=A(cur,t))
                %push
                top = top+1;
                stack(top) = t;
                Ghat = addedge(Ghat, cur,t, 1);
                visited(t) = 1;
            end
        end
    end
end

infect = cur_count*nn/rounds;
disp("%%%%%%%%%%%% after removal%%%%%%%%")
disp(infect);


