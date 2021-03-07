figure

line(P(end,:), param.z, 'Linewidth', 1.5, 'Color', 'b')
axis ij
ax1 = gca; % current axes
ax1.XColor = 'b';
ax1.YColor = 'b';
xlabel('Phytoplankton concentration [cells/ m3]')
ylabel( 'Depth [m]')



ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
ax2.XColor = 'r';
ax2.YColor = 'r';
line(N(end,:), param.z, 'Parent', ax2, 'Linewidth', 1.2, 'Color','r', 'LineStyle', '--')
axis ij
xlabel('Concentration of nutrients [mmol nutrient/m3]')




