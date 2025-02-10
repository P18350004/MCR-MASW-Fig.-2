
function  MASWaves_plot_dispersion_curve(f_curve0,c_curve0,lambda_curve0,...
    f_curve0_up,c_curve0_up,lambda_curve0_up,f_curve0_low,c_curve0_low,lambda_curve0_low,...
    type,up_low_boundaries,FigWidth,FigHeight,FigFontSize)

%Frequency vs. Rayleigh wave phase velocity
if strcmp(type,'f_c')
    
    % With upper/lower boundaries
    if strcmp(up_low_boundaries,'yes')
        hold on
        plot(f_curve0,c_curve0,'ko-','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k')
        plot(f_curve0_up,c_curve0_up,'r+--','MarkerSize',3,'MarkerFaceColor','r','MarkerEdgeColor','r')
        plot(f_curve0_low,c_curve0_low,'r+--','MarkerSize',3,'MarkerFaceColor','r','MarkerEdgeColor','r')
        legend({'Exp.','Exp. up/low'},'location','northeast','FontSize',FigFontSize)
        hold on
    end
    
    % Without upper/lower boundaries
    plot(f_curve0,c_curve0,'ko-','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k')
    if strcmp(up_low_boundaries,'no')
        legend('Exp.','location','northeast','FontSize',FigFontSize)
    end
    
    % Axis labels and axis limits
    set(gca, 'FontSize', FigFontSize)
    grid off
    xlabel('Frequency [Hz]','FontSize',FigFontSize,'Fontweight','normal')
    ylabel('Rayleigh wave velocity [m/s]','FontSize',FigFontSize,'Fontweight','normal')
    
    % Size of figure
    set(gcf,'units','centimeters')
    pos = [2, 2, FigWidth, FigHeight];
    set(gcf,'Position',pos)
    box on
    set(gca,'TickDir','in')
end

% Rayleigh wave phase velocity vs. wavelength                     
if strcmp(type,'c_lambda')
    
    % With upper/lower boundaries
    if strcmp(up_low_boundaries,'yes')
        hold on
        plot(c_curve0,lambda_curve0,'ko-','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k')
        plot(c_curve0_up,lambda_curve0_up,'r+--','MarkerSize',3,'MarkerFaceColor','r','MarkerEdgeColor','r')
        plot(c_curve0_low,lambda_curve0_low,'r+--','MarkerSize',3,'MarkerFaceColor','r','MarkerEdgeColor','r')
        legend({'Exp.','Exp. up/low'},'location','southwest','FontSize',FigFontSize)
        hold on
    end
    
    % Without upper/lower boundaries
    plot(c_curve0,lambda_curve0,'ko-','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k')
        if strcmp(up_low_boundaries,'no')
        legend('Exp.','location','southwest','FontSize',FigFontSize)
    end
    
    % Axis labels and axis limits
    set(gca, 'FontSize', FigFontSize)
    axis ij
    grid off
    xlabel('Rayleigh wave velocity [m/s]','FontSize',FigFontSize,'Fontweight','normal')
    ylabel('Wavelength [m]','FontSize',FigFontSize,'Fontweight','normal')
    
    % Size of figure
    set(gcf,'units','centimeters')
    pos = [2, 2, FigWidth, FigHeight];
    set(gcf,'Position',pos)
    box off
    set(gca,'TickDir','in')
end
hold off
end
