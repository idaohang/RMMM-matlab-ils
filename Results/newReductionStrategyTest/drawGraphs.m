st=6;
nd = 9;
for i = st:nd
    subplot(2,2,i-st+1);
    plot(noises,mean(squeeze(time1(i,:,:))'),'-x',noises,mean(squeeze(time2(i,:,:))'),'-o',noises,mean(squeeze(time3(i,:,:))'),'--');    
    leg = legend('LLL','LLL+Permu','LLL+Permu Babai','Location','NorthWest');
    xlab = xlabel('Noise Vector Sigma');
    ylab = ylabel('Average Time (s)');
    temp = strcat('Average Times over 200 Runs, Problem Size = ',num2str(sizes(i)));
    ttl = title(temp);
    set(xlab,'FontSize',16,'FontWeight','bold');
    set(ylab,'FontSize',16,'FontWeight','bold');
    set(ttl,'FontSize',16,'FontWeight','bold');
    set(leg,'FontSize',14,'FontWeight','bold');
end