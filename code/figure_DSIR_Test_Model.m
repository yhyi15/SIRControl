fileID = fopen('testModel/ds-upper.txt','r');
formatSpec = '%f';
greedy = fscanf(fileID,formatSpec);

fileID = fopen('testModel/ds-real.txt','r');
formatSpec = '%f';
rand_up = fscanf(fileID,formatSpec);


fileID = fopen('testModel/ds-markov.txt','r');
formatSpec = '%f';
maxd = fscanf(fileID,formatSpec);

k = size(greedy,1);
x = 0:1:k-1;
plot(x,greedy,'LineWidth',3);

hold on 
plot(x,rand_up,'LineWidth',2);

hold on 
plot(x,maxd,'color',[0.6350    0.0780    0.1840]);

set(gca,'FontSize',14,'FontWeight','bold');
xlabel('k','FontSize',14,'FontWeight','bold');
ylabel('Expected Infections','FontSize',14,'FontWeight','bold')
xlim([-1 17]);
ylim([2 7.5]);

legend('Upper','D-SIR','G-SIR');