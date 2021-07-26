function PRshowbcd(iterInfo,probInfo)
    % This function plots the results from IRbcd
    % Inputs:
    %       iterInfo: output from IRbcd
    %       probInfo: output from PRtomo_var
    
    figure(1), clf
    PRshowx(probInfo.true.x, probInfo.ProbParams)
    title('True Solution','fontsize', 20)

    figure(2), clf
    PRshowx(iterInfo.xtrueparam,  probInfo.ProbParams)
    title('Solution with True A','fontsize', 20)

    figure(3), clf
    PRshowx(iterInfo.xsols(:,1),  probInfo.ProbParams)
    title('Solution with inital Guess Parameters','fontsize', 20)

    figure(4), clf
    PRshowx(iterInfo.xsols(:,end), probInfo.ProbParams)
    title('Solution After BCD','fontsize', 20)

    figure(5), clf
    plot(iterInfo.xErrors,'-o','LineWidth',1.5);
    hold on 
    ax = gca;
    ax.FontSize = 15;
    plot(iterInfo.pErrors,'-*','LineWidth',1.5);
    plot(iterInfo.RErrors,'-d','LineWidth',1.5);
    plot(iterInfo.angleErrors,'-^','LineWidth',1.5);
    hold off
    legend('x Error Norms','p Error Norms', 'R Error Norms','Angle Error Norms','Fontsize',15);
    xlabel('Number of Iterations','Fontsize',20);
    ylabel('Relative Error','Fontsize',20);
    title('Errors at each Iteration','Fontsize',20);
end