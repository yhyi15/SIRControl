fileID = fopen('results/icGreedyResult_BBC_23.txt','r');
formatSpec = '%f';
greedy = fscanf(fileID,formatSpec);

fileID = fopen('results/icRandResult_BBC_23.txt','r');
formatSpec = '%f';
rand_up = fscanf(fileID,formatSpec);


fileID = fopen('results/icMaxDResult_BBC_23.txt','r');
formatSpec = '%f';
maxd = fscanf(fileID,formatSpec);

k = size(greedy,1);
%k=221;
x = 0:1:k-1;
plot(x,greedy,'LineWidth',3);

hold on 
plot(x,maxd,'color',[0.6350    0.0780    0.1840],'LineWidth',2);

hold on 
plot(x,rand_up);



set(gca,'FontSize',14,'FontWeight','bold');
xlabel('k','FontSize',14,'FontWeight','bold');
ylabel('Expected Number of Infections','FontSize',14,'FontWeight','bold');
xlim([(min(x)-1) (max(x)+1)]);
ylim([(min(greedy-1)) (max(greedy)+1)]);
title('ER Network', 'FontSize',14);

legend('Greedy','GreedyRho','Random');