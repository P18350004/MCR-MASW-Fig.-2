
function [fplot, cplot, Aplot] = MASWaves_plot_dispersion_image_2D(f,c,A,fmin,fmax,resolution,FigWidth,FigHeight,FigFontSize)

% Limits of frequency axis
[~,no_fmin] = (min(abs(f(:,1)-fmin)));
[~,no_fmax] = (min(abs(f(:,1)-fmax)));

% Select data corresponding to frequnecy range [fmin,fmax]
% Compute absolute value (length) of complex numbers
Aplot = A(no_fmin:no_fmax,:); 
fplot = f(no_fmin:no_fmax,:);
cplot = c(no_fmin:no_fmax,:);

% Plot the 2D dispersion image
[~,ch] = contourf(fplot,cplot,Aplot,resolution);
set(ch,'LineStyle','none');
set(ch,'edgecolor','none');
colormap(jet)
shading flat
grid off

% Axis limits and axis labels
set(gca,'XTick',0:10:fmax+0.01);
set(gca,'FontSize',FigFontSize,'FontName','Times New Roman');
xlim([fmin fmax])
xlabel('Frequency [Hz] ','FontSize',FigFontSize,'Fontweight','normal','color','k')
ylabel('Phase velocity [m/s] ','FontSize',FigFontSize,'Fontweight','normal','color','k')

% Size of figure
set(gcf,'units','centimeters')
pos = [2, 2, FigWidth, FigHeight]; 
set(gcf,'Position',pos)
box on
set(gca,'TickDir','in')

% Colorbar
c = colorbar('location','NorthOutside','color','k');
c.FontSize = FigFontSize;
c.Label.String = 'Normalized amplitude';
c.Label.FontSize = FigFontSize;

% Width and location of colorbar
axpos = get(gca,'Position');
cpos = get(c,'Position');
cpos(2) = cpos(2) - 0.5*cpos(4);
cpos(4) = 0.5*cpos(4);
set(c,'Position',cpos)
set(gca,'Position',axpos)

end