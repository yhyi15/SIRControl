clear;
clc;
N=500;
m=3;
m0=4;

A=BAGraph(N,m,m0);
G = graph(A);
edges=table2array(G.Edges);

m=numedges(G);
outFile = fopen('baGraph.txt','w');

for i =1:m
    fprintf(outFile,'%d\t%d\n', edges(i,1),edges(i,2));
end
fclose(outFile);