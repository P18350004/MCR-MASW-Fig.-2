
function [fplot,cplot,Aplot] = MASWaves_plot_dispersion_image_3D(f,c,A,fmin,fmax,FigWidth,FigHeight,FigFontSize)
        
% Limits of frequency axis
[~,no_fmin] = (min(abs(f(:,1)-fmin)));
[~,no_fmax] = (min(abs(f(:,1)-fmax)));

% Select data corresponding to frequnecy range [fmin,fmax]
% Compute absolute value (length) of complex numbers
Aplot = A(no_fmin:no_fmax,:); 
fplot = f(no_fmin:no_fmax,:);
cplot = c(no_fmin:no_fmax,:);

% Plot the 3D dispersion image
surf(fplot,cplot,Aplot);
colormap(jet)
shading interp
grid on

% Axis limits and axis labels
set(gca,'FontSize',FigFontSize);
xlabel('Frequency [Hz]','FontSize',FigFontSize,'Fontweight','normal')
ylabel('Phase velocity [m/s]','FontSize',FigFontSize,'Fontweight','normal')
zlabel('Normalized amplitude','FontSize',FigFontSize,'Fontweight','normal')

% Size of figure
set(gcf,'units','centimeters')
pos = [2, 2, FigWidth, FigHeight]; 
set(gcf,'Position',pos)
box off
set(gca,'TickDir','out')

% Colorbar
c = colorbar('location','NorthOutside','color','k');
c.FontSize = FigFontSize;

% Width and location of colorbar
axpos = get(gca,'Position');
cpos = get(c,'Position');
cpos(2) = cpos(2) - 0.5*cpos(4);
cpos(4) = 0.5*cpos(4);
set(c,'Position',cpos)
set(gca,'Position',axpos)

% View
az = 45;
el = 45;
view(az,el)
end