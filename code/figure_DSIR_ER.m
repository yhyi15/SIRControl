fileID = fopen('results/dsGreedyResult_bbc_13.txt','r');
formatSpec = '%f';
greedy = fscanf(fileID,formatSpec);

fileID = fopen('results/dsMaxDResult_bbc_13.txt','r');
formatSpec = '%f';
maxd = fscanf(fileID,formatSpec);

fileID = fopen('results/dsRandResult_bbc_13.txt','r');
formatSpec = '%f';
rand_up = fscanf(fileID,formatSpec);



k = size(greedy,1);
%k=221;
x = 1:1:k;
plot(x,greedy(1:k,:),'LineWidth',3)

hold on 
plot(x,maxd(1:k,:),'LineWidth',2,'color', [0.6350    0.0780    0.1840])

hold on 
plot(x,rand_up(1:k,:))

set(gca,'FontSize',14,'FontWeight','bold');
xlabel('k','FontSize',14,'FontWeight','bold');
ylabel('Upper Bound of Infections','FontSize',14,'FontWeight','bold');
title('Haslemere Network', 'FontSize',14);
xlim([(min(x)-10) (max(x)+10)]);
ylim([2 (max(greedy(1:k,:)+0.5))]);

legend('Greedy','Max-Degree','Random')