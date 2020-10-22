clear;
clc;

G = graph;
blockSize = 100;
numBlock = 5;

inProb = 0.5;
outProb = 0.4/(numBlock -1);
n = 100;
for  u =1:5
    for v = 1:5
        if u == v
            prob = 2*inProb*log(n)/blockSize;
        else
            prob = 2*(0.8+0.2*rand())*outProb*log(n)/blockSize;
        end
        for i = 1:n-1
            for j = i+1:n
                    edgeCoin = binornd(1, prob);
                if edgeCoin ==1
                    G = addedge(G, (u-1)*blockSize+i,(v-1)*blockSize+j, 1);%/(2.5*log(n)));
                end
            end
        end
    end
end

edges=table2array(G.Edges);

m=numedges(G);
outFile = fopen('SBMGraph.txt','w');

for i =1:m
    fprintf(outFile,'%d\t%d\n', edges(i,1),edges(i,2));
end
fclose(outFile);