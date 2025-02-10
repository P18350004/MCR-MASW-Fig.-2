
function MASWaves_plot_data(u,N,dx,x1,L,t,Tmax,du,FigWidth,FigHeight,FigFontSize)

for j = 1:N
    signal_plot(:,j) = u(:,j) + (j-1)*du*dx;
    plot(signal_plot(:,j),t,'k');
    hold on
end
set(gca,'YDir','reverse')
set(gca, 'FontSize', FigFontSize)

% Axis limits and axis labels
axis([-5*dx*du,L*du+5*dx*du,0,Tmax]) 
set(gca,'XTick',0:6*dx*du:L*du,'XTickLabel',x1:(6*dx):x1+L)
set(gca,'TickDir','out')
box off
xlabel('Distance from source [m]','FontSize',FigFontSize)
ylabel('Time [s]','FontSize',FigFontSize)

% Size of figure
set(gcf,'units','centimeters')
pos = [2, 2, FigWidth, FigHeight]; 
set(gcf,'Position',pos)

hold off
end
