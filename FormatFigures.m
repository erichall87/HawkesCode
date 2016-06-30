function out = FormatFigures(ax,save, name)
%FORMATFIGURES Takes as input an axis handle and formats the axis in a
%way consistent with our paper. The axis should be the only plot on the
%figure to work correctly. Additional optional arguments about whether
%or not to save the figure, and the name of the files to save
if nargin == 1
    save = false;
elseif nargin == 2
    name = 'TempFig';
end
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
if save
    print(fig,name,'-dpdf')
    savefig(strcat(name,'.fig'))
end