fileID = fopen('dsGreedyResult_2.txt','r');
formatSpec = '%f';
greedy = fscanf(fileID,formatSpec);

fileID = fopen('dsRandResult_2.txt','r');
formatSpec = '%f';
rand_up = fscanf(fileID,formatSpec);

fileID = fopen('dsMaxDResult_2.txt','r');
formatSpec = '%f';
maxd = fscanf(fileID,formatSpec);

k = size(A,1);
x = 1:1:k;
plot(x,greedy,'LineWidth',3)

hold on 
plot(x,maxd,'LineWidth',2)

hold on 
plot(x,rand_up,'color',[0.6350    0.0780    0.1840])

set(gca,'FontSize',14,'FontWeight','bold');
xlabel('k','FontSize',14,'FontWeight','bold')
ylabel('Upper Bound of Expected Infections','FontSize',14,'FontWeight','bold')
xlim([-10 540]) 
ylim([3 21])

legend('Greedy','Max-Degree','Random')