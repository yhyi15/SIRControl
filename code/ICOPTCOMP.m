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
G = graph(rdEdges(:,1),rdEdges(:,2),ones(m,1),n);
% for i = 1:n-1
%     for j = i+1:n
%         edgeCoin = binornd(1,2*log(n)/n);
%         if edgeCoin ==1
%             G = addedge(G, i,j, 1);%/(2.5*log(n)));
%         end
%     end
% end

preEdges = rdEdges;

while max(degree(G))>8
    dmax=0;
    for e=1:m
        dupdate = max(degree(G, preEdges(e, 1)), degree(G, preEdges(e, 2)));
        if dupdate>dmax
            dmax=dupdate;
            idx=e;
        end
    end
    G = rmedge(G,preEdges(idx,1),preEdges(idx,2));
    preEdges(idx,:)=[];
    m=m-1;
end

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
%S = randsample(gccSize,s);
S = [442,39,129,307,370]';
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
%qidx = randsample(m, q);
qidx = readmatrix('data/edgeSet.txt');
Q = edgelist(qidx, :);



%k = floor(qsize/2);
k=220;
%k =150;
%k=220;

%matrices
A = zeros(nn+1,nn+1);
A(1:nn,1:nn) = 1.1/log(nn)*adjacency(GCC,'weighted');
for sidx = 1:s
    A(nn+1, S(sidx))=1;
    A(S(sidx), nn+1)=1;
end
betaList = beta* ones(nn,1) + 0.8*beta*rand([nn,1]);
deltaList = delta * ones(nn,1);
B = diag(betaList);
D = diag(deltaList);
X0 = diag(x0);
R0 = diag(r0);
I = eye(nn);
%M = I - D + (I-X0-R0)*B*A;


%calculating original expected infections
%randomly pick a node from the network
eps = 0.2;
rounds = 0.1*nn*log(nn)/(eps^2);

count = 0;

disp(rounds);


P = zeros(k,3);
Ggreedy = graph(A);
GgreedyRho = graph(A);
for sidx =1:s
    Ggreedy = addedge(Ggreedy, nn+1, S(sidx),1);
end


disp("%%%%%%%%%%%% Original %%%%%%%%%%%%%%%%%")
% cur_count = 0;
% for round = 1:10*rounds
%     terminal = randsample(nn,1);
%     
%     stack = zeros(nn, 1);
%     stack(1) = terminal;
%     top = 1;
%     visited = zeros(nn+1,1);
%     visited(terminal)=1;
%     while top>0
%         %pop
%         cur = stack(top);
%         if cur == nn+1
%             cur_count = cur_count+1;
%             break;
%         end
%         top = top -1;
%         %if activated(cur)==0
%             nb = neighbors(Ggreedy, cur);
%             for iter = 1:size(nb,1)
%                 %push
%                 t = nb(iter);
%                 contCoin = rand();
%                 if (visited(t)==0 && contCoin<=A(cur,t))
%                     %push
%                     top = top+1;
%                     stack(top) = t;
%                     visited(t)=1;
%                 end
%             end
%             %activated(cur)=1;
%         %end
%     end
% %     for seed_iter= 1:s
% %         if activated(S(seed_iter))==1
% %             cur_count=cur_count+1;
% %             break;
% %         end
% %     end
% end
% infect = cur_count*nn/10/rounds;
infect = icEval2(Ggreedy, nn, rounds,S ,s, A);
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
        Ghat = graph;
        %contEdge = zeros(nn,2);
        %contNumEdges = 0;
        %if(mod(round,100)==0)
        %    disp(round/rounds);
        %end
        terminal = randsample(nn,1);
        
        stack = zeros(nn, 1);
        stack(1) = terminal;
        top = 1;
        visited = zeros(nn+1,1);
        visited(terminal)=1;
        while top>0
            %pop
            cur = stack(top);
            if cur == nn+1
                succ=1;
                cur_count=cur_count+1;
                break;
            end
            top = top -1;
            %if activated(cur)==0
                nb = neighbors(Ggreedy, cur);
                for iter = 1:size(nb,1)
                    %push
                    t = nb(iter);
                    contCoin = rand();
                    if (contCoin<=A(cur,t))
                        Ghat = addedge(Ghat, cur,t);
                        %push
                        if visited(t)==0
                            top = top+1;
                            stack(top) = t;
                        end
                        %contNumEdges = contNumEdges + 1;
                        %contEdge(contNumEdges,:) = [cur, t];
                    end
                end
                %activated(cur)=1;
            %end
        end
%         for iter_c = 1:s
%             if activated(S(iter_c))==1
%                 succ = 1;
%                 cur_count=cur_count+1;
%                 break;
%             end
%         end
        
        if succ ==0
            continue;
        end
        %slist = Q(:,1);
        %tlist = Q(:,2);
        %contEdge = table2array(Ghat.Edges);
        %contNumEdges = numedges(Ghat);
        %for sidx = 1:s
        %    Ghat = addedge(Ghat, S(sidx), nn+1);
        %end
        for e =1:q
            %Sset = S;
            sn = Q(e,1);
            tn = Q(e,2);
            %in_h = findedgeCont(contEdge, contNumEdges,sn,tn);
            in_h = findedge(Ghat,sn,tn);
            if in_h ==0
                count(e) = count(e)+succ;
                continue;
            end
            stack_h = zeros(nn+1, 1);
            stack_h(1)= terminal;
            top_h = 1;
            visited_h = zeros(nn+1,1);
            visited_h(terminal)=1;
            reached = 0;
            while top_h>0
                cur_h = stack_h(top_h);
                if cur_h == nn+1
                    reached = 1;
                    break;
                end
                top_h = top_h -1;
                %if activated_h(cur_h)==0
                    %nb_h = neighborsCont(contEdge, contNumEdges, cur_h);
                    nb_h = neighbors(Ghat, cur_h);
                    %disp(size(nb_h));
                    for iter_h = 1:size(nb_h,1)
                        %push
                        t_h = nb_h(iter_h);
                        if (visited_h(t_h)==0 && ~((cur_h==sn && t_h == tn) || (cur_h==tn && t_h == sn)))
                            %push
                            top_h = top_h+1;
                            stack_h(top_h) = t_h;
                            visited_h(t_h)=1;
                        end
                    end
                %end
            end
            %for iter_h_ss = 1:s
            %    if activated_h(S(iter_h_ss))==1
            %        count(e)=count(e)+1;
            %        break;
            %    end
            %end
            count(e) = count(e)+reached;
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
    
    greedyRst(i+1) = icEval2(Ggreedy, nn, rounds,S ,s, A)-s;
    disp(greedyRst(i+1));
end


cur_count = 0;
for round = 1:10*rounds
    terminal = randsample(nn,1);
    
    stack = zeros(nn, 1);
    stack(1) = terminal;
    top = 1;
    visited = zeros(nn+1,1);
    visited(terminal)=1;
    while top>0
        %pop
        cur = stack(top);
        if cur == nn+1
            cur_count = cur_count+1;
            break;
        end
        top = top -1;
        %if activated(cur)==0
            nb = neighbors(Ggreedy, cur);
            for iter = 1:size(nb,1)
                %push
                t = nb(iter);
                contCoin = rand();
                if (visited(t)==0 && contCoin<=A(cur,t))
                    %push
                    top = top+1;
                    stack(top) = t;
                    visited(t)=1;
                end
            end
            %activated(cur)=1;
        %end
    end
%     for seed_iter= 1:s
%         if activated(S(seed_iter))==1
%             cur_count=cur_count+1;
%             break;
%         end
%     end
end
infect = cur_count*nn/10/rounds;
disp("greedy after removal:")
disp(infect-s);


% greedyRstRho = zeros(k+1,1);
% greedyRstRho(1) = originalRst;
% 
% Q4=edgelist(qidx, :);
% q4=qsize;
% 
% disp("%%%%%%%%%%%%%%% Greedy RHO%%%%%%%%%%%%%%%%%%%%%");
% disp(strcat("k: ", num2str(k)));
% for i = 1:k
%     disp(strcat("deleting the ",num2str(i),"-th edge"));
%     cur_count = 0;
%     count = zeros(q,1);
%     round = 1;
%     %numcyc = 0;
%     edgelistTemp = table2array(GgreedyRho.Edges);
%     m = numedges(GgreedyRho);
%     while round<rounds
%         succ = 0;
%         %Ghat = graph;
%         %contEdge = zeros(nn,2);
%         %contNumEdges = 0;
%         %if(mod(round,100)==0)
%         %    disp(round/rounds);
%         %end
%         
%         terminal = randsample(nn,1);
%         Ghat = graph;
%         Ghat = addnode(Ghat, nn);
%         for ed =1: m
%             sst = edgelistTemp(ed,1);
%             ttt = edgelistTemp(ed,2);
%             contCoin2 = rand();
%             if (contCoin2<=A(sst,ttt) && sst~=nn+1 && ttt~=nn+1 )
%                 Ghat = addedge(Ghat, sst, ttt);
%             end
%         end
%         [bins, binsizes] = conncomp(Ghat);
%         if max(binsizes)>=9*log(nn)*25
%             continue;
%         end
%         collision = 0;
%         for si = 1:s-1
%             for sj = si+1:s
%                 if bins(S(si))==bins(S(sj))
%                     collision = 1;
%                 end
%             end
%         end
%         if collision ==1
%             continue;
%         end
%         round =round+1;
%         GhatIdx = find(bins==bins(terminal));
%         %Ghat = subgraph(Ghat, GhatIdx);
%         nhat = numnodes(Ghat);
%         mhat = numedges(Ghat);
%         edgelisthat = table2array(Ghat.Edges);
%         Ghhat=graph;
%         for editer = 1:mhat
%             if bins(edgelisthat(editer,1))== bins(terminal)
%                 Ghhat = addedge(Ghhat, edgelisthat(editer,1), edgelisthat(editer,2));
%             end
%         end
%         Ghat = Ghhat;
%         %disp("cycles:")
%         %disp(mhat - nhat+1);
%         
%         if nhat<=1 || mhat >= nhat
%            succ=0;
%            continue;
%         end
%         for ssidx =1:s
%             if bins(S(ssidx))==bins(terminal)
%                 succ =1;
%             end
%             Ghat = addedge(Ghat, nn+1, S(ssidx));
%         end
%         if succ ==0
%             continue;
%         end
% 
%         cur_count = cur_count+1;
%         [path , dis]= shortestpath(Ghat, terminal, nn+1, 'Method','unweighted');
%         for e =1:q4
%             %Sset = S;
%             sn = Q4(e,1);
%             tn = Q4(e,2);
%             reached = 1;
%             for einp = 1:dis
%                 if ((path(einp) == sn && path(einp+1)== tn)...
%                         ||(path(einp) == tn && path(einp+1)== sn))
%                     reached = 0;
%                 end
%             end
%             count(e)=count(e)+reached;
%         end
%     end
%     %disp("%%%%%%%%%%%%%%%")
%     %disp(cur_count);
%     [minv, idx] = min(count);
%     %disp(size(count));
%     %disp(numcyc);
%     %disp(minv);
%     P(i,:) = Q4(idx,:);
%     GgreedyRho = rmedge(GgreedyRho,Q4(idx,1),Q4(idx,2));
%     Q4(idx,:)=[];
%     q4=q4-1;
%     
%     greedyRstRho(i+1) = icEval2(GgreedyRho, nn, rounds,S ,s, A)-s;
%     disp(greedyRstRho(i+1));
% end
% 

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
    randRst(i+1) = icEval2(Grand, nn, rounds,S ,s, A)-s;
end


cur_count = 0;
for round = 1:10*rounds
    terminal = randsample(nn,1);
    
    stack = zeros(nn, 1);
    stack(1) = terminal;
    top = 1;
    visited = zeros(nn+1,1);
    visited(terminal)=1;
    while top>0
        %pop
        cur = stack(top);
        if cur == nn+1
            cur_count = cur_count+1;
            break;
        end
        top = top -1;
        %if activated(cur)==0
            nb = neighbors(Grand, cur);
            for iter = 1:size(nb,1)
                %push
                t = nb(iter);
                contCoin = rand();
                if (visited(t)==0 && contCoin<=A(cur,t))
                    %push
                    top = top+1;
                    stack(top) = t;
                    visited(t)=1;
                end
            end
            %activated(cur)=1;
        %end
    end
%     for seed_iter= 1:s
%         if activated(S(seed_iter))==1
%             cur_count=cur_count+1;
%             break;
%         end
%     end
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
    for e = 1:q3
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
    maxDRst(i+1) = icEval2(GmaxD, nn, rounds,S ,s, A)-s;
end


cur_count = 0;
for round = 1:10*rounds
    terminal = randsample(nn,1);
    
    stack = zeros(nn, 1);
    stack(1) = terminal;
    top = 1;
    visited = zeros(nn+1,1);
    visited(terminal)=1;
    while top>0
        %pop
        cur = stack(top);
        if cur == nn+1
            cur_count = cur_count+1;
            break;
        end
        top = top -1;
        %if activated(cur)==0
            nb = neighbors(GmaxD, cur);
            for iter = 1:size(nb,1)
                %push
                t = nb(iter);
                contCoin = rand();
                if (visited(t)==0 && contCoin<=A(cur,t))
                    %push
                    top = top+1;
                    stack(top) = t;
                    visited(t)=1;
                end
            end
            %activated(cur)=1;
        %end
    end
%     for seed_iter= 1:s
%         if activated(S(seed_iter))==1
%             cur_count=cur_count+1;
%             break;
%         end
%     end
end
infect = cur_count*nn/10/rounds;
disp("maxD after removal:")
disp(infect-s);

disp("printing results")
outFile0 = fopen('results/basicInfo_BBC_21.txt','w');
outFile1 = fopen('results/icGreedyResult_BBC_23.txt','w');
outFile2 = fopen('results/icRandResult_BBC_23.txt','w');
outFile3 = fopen('results/icMaxDResult_BBC_23.txt','w');
%outFile4 = fopen('results/icGreedyResultRho_BBC_23.txt','w');

numRmv = (1:k+1)'-ones(k+1,1);
fprintf(outFile0, '%d\n%d\n%d\n%d\n%d\n',nn, m, qsize, k, s);
for i = 1:k+1
    fprintf(outFile1,'%5f\n', greedyRst(i));
    fprintf(outFile2,'%5f\n', randRst(i));
    fprintf(outFile3,'%5f\n', maxDRst(i));
    %fprintf(outFile4,'%5f\n', greedyRstRho(i));
end
fclose(outFile1);
fclose(outFile2);
fclose(outFile3);
disp("finished writing");    

