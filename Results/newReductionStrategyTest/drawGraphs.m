for i = 6:length(sizes)
    subplot(2,2,i-6+1);
    plot(noises,mean(squeeze(time1(i,:,:))'),noises,mean(squeeze(time2(i,:,:))'));    
    legend('LLL','LLL+Permu','Location','NorthWest');
    xlabel('Noise Vector Sigma');
    ylabel('Average Time (s)');
    temp = strcat('Average Times over 200 Runs, Problem Size = ',num2str(sizes(i)));
    title(temp);
end