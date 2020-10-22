% optimizing deterministic spread model
% input: a contact graph
% algorithm: minimize upper bound 
%\hat{sigma}= 1^T *(M+D-I)(I-M)^(-1) * x(0)
% output: sigma =  |(x(t)+ r(t))-(x(0)+ r(0))|_1 and \hat{sigma}
clear;
clc;

rdEdges = readmatrix('data/BBC_contact_network.csv');

%create an arbitrary graph
%
n = max(max(rdEdges));
m = size(rdEdges, 1);
G = graph(rdEdges(:,1),rdEdges(:,2),0.15*ones(m,1),n);
% for i = 1:n-1
%     for j = i+1:n
%         edgeCoin = binornd(1,2*log(n)/n);
%         if edgeCoin ==1
%             G = addedge(G, i,j, 1);%/(2.5*log(n)));
%         end
%     end
% end

preEdges = rdEdges;

% while max(degree(G))>8
%     dmax=0;
%     for e=1:m
%         dupdate = max(degree(G, preEdges(e, 1)), degree(G, preEdges(e, 2)));
%         if dupdate>dmax
%             dmax=dupdate;
%             idx=e;
%         end
%     end
%     G = rmedge(G,preEdges(idx,1),preEdges(idx,2));
%     preEdges(idx,:)=[];
%     m=m-1;
% end

% % take the (presumably) largest connected component
% [bins, binsizes] = conncomp(G);
% gccSize = max(binsizes);
% idx = binsizes(bins) == gccSize;
% SG = subgraph(G, idx);
% %reorder nodes in the gcc
% order = 1:gccSize;
% GCC = reordernodes(SG, order);
GCC = G;
gccSize = numnodes(GCC);
maxd = max(degree(G));
disp(maxd);
beta = 1/(maxd+1);
disp(beta);
delta = 1;


%set initial conditions
s = 5;
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
k = floor(qsize/2);

%matrices
A = adjacency(GCC,'weighted');
betaList = beta* ones(nn,1) + 0.8*beta*rand([nn,1]);
deltaList = delta * ones(nn,1);
B = diag(betaList);
D = diag(deltaList);
X0 = diag(x0);
R0 = diag(r0);
I = eye(nn);
M = I - D + (I-X0-R0)*B*A;


%calculating original expected infections
%randomly pick a node from the network
eps = 0.2;
rounds = 0.01*nn*log(nn)/(eps^2);

count = 0;

disp(rounds);


P = zeros(k,3);
Ggreedy = graph(A);


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
                        if 1>0 %(~((cur_h==sn && t_h == tn) || (cur_h==tn && t_h == sn)))
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
    P(k,:) = Q(idx,:);
    Ggreedy = rmedge(Ggreedy,Q(idx,1),Q(idx,2));
    Q(idx,:)=[];
    q=q-1;
    
    greedyRst(i+1) = icEval(Ggreedy, nn, rounds,S ,s, A)-s;
end


cur_count = 0;
for round = 1:10*rounds
    %Ghat = graph;
    %if(mod(round,5000)==0)
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
disp("greedy after removal:")
disp(infect-s);


randRst = zeros(k+1,1);
randRst(1) = originalRst;

%to compare with random choices
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
disp("random choices:");
q2=qsize;
P2=zeros(k,3);
Q2=edgelist(qidx, :);
Grand = graph(A);
for i = 1:k
    disp(strcat("deleting the ",num2str(i),"-th edge"));
    %disp(i)
    %iterate over all remaining candidates
    %curr = ones(1, nn)* (M+D-I)*((I-M)\x0);
    e =randi(q2,1);
    %disp(e)
    
    Grand = rmedge(Grand, Q2(e,1), Q2(e,2));
    P2(i,:) = Q2(e,:);
    Q2(e,:)=[];
    q2=q2-1;
    randRst(i+1) = icEval(Grand, nn, rounds,S ,s, A)-s;
end


cur_count = 0;
for round = 1:10*rounds
    %Ghat = graph;
    %if(mod(round,5000)==0)
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
        if(activated(cur)==0)
            nb = neighbors(Grand, cur);
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
disp("rand after removal:")
disp(infect-s);


maxDRst = zeros(k+1,1);
maxDRst(1) = originalRst;

%to compare with max degree
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
disp("max dgree:");
q3=qsize;
P3=zeros(k,3);
Q3=edgelist(qidx, :);
GmaxD = graph(A);
for i = 1:k
    disp(strcat("deleting the ",num2str(i),"-th edge"));
    %%%%%%%
    dmax = 0;
    for e = 1:q
        dupdate = max(degree(GmaxD, Q3(e, 1)), degree(GmaxD, Q3(e, 2)));
        %dupdate = max(degree(GCCCopy, Q3(e, 1)), degree(GCCCopy, Q3(e, 2)));
        if dupdate>dmax
            choice = e;
            dmax = dupdate;
        end
    end
    
    GmaxD = rmedge(GmaxD, Q3(choice,1), Q3(choice,2));
    P3(i,:) = Q3(choice,:);
    Q3(choice,:)=[];
    q3=q3-1;
    maxDRst(i+1) = icEval(GmaxD, nn, rounds,S ,s, A)-s;
end


cur_count = 0;
for round = 1:10*rounds
    %Ghat = graph;
    %if(mod(round,5000)==0)
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
        if(activated(cur)==0)
            nb = neighbors(GmaxD, cur);
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
disp("maxD after removal:")
disp(infect-s);

disp("printing results")
outFile0 = fopen('results/basicInfo_BBC_0.txt','w');
outFile1 = fopen('results/icGreedyResult_BBC_0.txt','w');
outFile2 = fopen('results/icRandResult_BBC_0.txt','w');
outFile3 = fopen('results/icMaxDResult_BBC_0.txt','w');

numRmv = (1:k+1)'-ones(k+1,1);
fprintf(outFile0, '%d\n%d\n%d\n%d\n%d\n',nn, m, q, k, s);
for i = 1:k+1
    fprintf(outFile1,'%5f\n', greedyRst(i));
    fprintf(outFile2,'%5f\n', randRst(i));
    fprintf(outFile3,'%5f\n', maxDRst(i));
end
fclose(outFile1);
fclose(outFile2);
fclose(outFile3);
disp("finished writing")    

