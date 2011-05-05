for i = 1:length(noises)
    subplot(4,2,i);
    plot(sizes,mean(squeeze(expand1(:,i,:))'),sizes,mean(squeeze(expand2(:,i,:))'),sizes,mean(squeeze(expand3(:,i,:))'));
    legend('LLL','LLL+Permu','TestReduction','Location','NorthWest');
    xlabel('Problem Size');
    ylabel('Nodes Expanded');
    temp = strcat('Noise = ',num2str(noises(i)));
    title(temp);
end