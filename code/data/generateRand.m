clear;
clc;
G = graph;

n = 500;
%G = graph;
for i = 1:n-1
    for j = i+1:n
        %edgeCoin = binornd(1,1.75*log(n/5)/n);
        edgeCoin = binornd(1,0.8*log(n)/n);
        if edgeCoin ==1
            G = addedge(G, i,j, 1);%/(2.5*log(n)));
        end
    end
end

edges=table2array(G.Edges);

m=numedges(G);
outFileR = fopen('randGraph.txt','w');

for i =1:m
    fprintf(outFileR,'%d\t%d\n', edges(i,1),edges(i,2));
end
fclose(outFileR);