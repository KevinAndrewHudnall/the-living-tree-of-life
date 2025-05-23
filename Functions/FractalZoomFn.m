function FractalZoomFn(Point_0, Contraction_Factor)

% This function zooms through the multifractal tree. There is the option
% (commented out here) to save an image at each zoomed in level for
% purposes of making an animation.

% get the figure and axes handles
figure(1) = gcf;
hAx  = gca;

%fullscreen
set(gcf,'units','normalized','outerposition',[0 0 1 1])

% Suppress the axes
ax = gca; 
ax.Visible = 'off';

% % set the axes to full screen
% set(hAx,'Unit','normalized','Position',[0 0 1 1]);

%hide the toolbar
set(gcf,'menubar','none')

zoom on

node_coord = Point_0;

xscale = xlim;
yscale = ylim;

xcenter = (xscale(2) + xscale(1)) / 2;
ycenter = (yscale(2) + yscale(1)) / 2;

xmove = abs(node_coord(1) - xcenter);
ymove = abs(node_coord(2) - ycenter);

%Center the x coordinates
if ((0 < xcenter) && (xcenter < node_coord(1)))
    xlim([(xscale(1) + xmove), (xscale(2) + xmove)]);
elseif ((0 < node_coord(1)) && (node_coord(1) < xcenter))
    xlim([(xscale(1) - xmove), (xscale(2) - xmove)]);
elseif ((node_coord(1) < xcenter) && (xcenter < 0))
    xlim([(xscale(1) - xmove), (xscale(2) - xmove)]);
elseif ((xcenter < node_coord(1)) && (node_coord(1) < 0))
    xlim([(xscale(1) + xmove), (xscale(2) + xmove)]);
elseif ((xcenter < 0) && (0 < node_coord(1)))
    xlim([(xscale(1) + xmove), (xscale(2) + xmove)]);
elseif ((node_coord(1) < 0) && (0 < xcenter))
    xlim([(xscale(1) - xmove), (xscale(2) - xmove)]);
end    

%Center the y coordinates
if ((0 < ycenter) && (ycenter < node_coord(2)))
    ylim([(yscale(1) + ymove), (yscale(2) + ymove)]);
elseif ((0 < node_coord(2)) && (node_coord(2) < ycenter))
    ylim([(yscale(1) - ymove), (yscale(2) - ymove)]);
elseif ((node_coord(2) < ycenter) && (ycenter < 0))
    ylim([(yscale(1) - ymove), (yscale(2) - ymove)]);
elseif ((ycenter < node_coord(2)) && (node_coord(2) < 0))
    ylim([(yscale(1) + ymove), (yscale(2) + ymove)]);
elseif ((ycenter < 0) && (0 < node_coord(2)))
    ylim([(yscale(1) + ymove), (yscale(2) + ymove)]);
elseif ((node_coord(2) < 0) && (0 < ycenter))
    ylim([(yscale(1) - ymove), (yscale(2) - ymove)]);
end    

% %Capture the graph at current scale
% PictureName = sprintf('Random_Multifractal_Tree_Zoom_%d', 1);
% print(gcf, PictureName, '-dpng')

% pause for visualization. Comment out if capturing pictures instead
pause(0.1)

for k = 2:620 % adjust depth as needed

    %Scale the x coordinate
    xscale = Contraction_Factor * xlim;
    xlim([xscale(1), xscale(2)]);

    %Scale the y coordinate
    yscale = Contraction_Factor * ylim;
    ylim([yscale(1), yscale(2)]);

    xscale = xlim;
    yscale = ylim;

    xcenter = (xscale(2) + xscale(1)) / 2;
    ycenter = (yscale(2) + yscale(1)) / 2;

    xmove = abs(node_coord(1) - xcenter);
    ymove = abs(node_coord(2) - ycenter);
    
    %Center the x coordinates
    if ((0 < xcenter) && (xcenter < node_coord(1)))
        xlim([(xscale(1) + xmove), (xscale(2) + xmove)]);
    elseif ((0 < node_coord(1)) && (node_coord(1) < xcenter))
        xlim([(xscale(1) - xmove), (xscale(2) - xmove)]);
    elseif ((node_coord(1) < xcenter) && (xcenter < 0))
        xlim([(xscale(1) - xmove), (xscale(2) - xmove)]);
    elseif ((xcenter < node_coord(1)) && (node_coord(1) < 0))
        xlim([(xscale(1) + xmove), (xscale(2) + xmove)]);
    elseif ((xcenter < 0) && (0 < node_coord(1)))
        xlim([(xscale(1) + xmove), (xscale(2) + xmove)]);
    elseif ((node_coord(1) < 0) && (0 < xcenter))
        xlim([(xscale(1) - xmove), (xscale(2) - xmove)]);
    end    

    %Center the y coordinates
    if ((0 < ycenter) && (ycenter < node_coord(2)))
        ylim([(yscale(1) + ymove), (yscale(2) + ymove)]);
    elseif ((0 < node_coord(2)) && (node_coord(2) < ycenter))
        ylim([(yscale(1) - ymove), (yscale(2) - ymove)]);
    elseif ((node_coord(2) < ycenter) && (ycenter < 0))
        ylim([(yscale(1) - ymove), (yscale(2) - ymove)]);
    elseif ((ycenter < node_coord(2)) && (node_coord(2) < 0))
        ylim([(yscale(1) + ymove), (yscale(2) + ymove)]);
    elseif ((ycenter < 0) && (0 < node_coord(2)))
        ylim([(yscale(1) + ymove), (yscale(2) + ymove)]);
    elseif ((node_coord(2) < 0) && (0 < ycenter))
        ylim([(yscale(1) - ymove), (yscale(2) - ymove)]);
    end    
    
%     %Capture the graph at current scale
%     PictureName = sprintf('Random_Multifractal_Tree_Zoom_%d', k);
%     print(gcf, PictureName, '-dpng')

    % Pause for visualization. Comment out if capturing pictures instead.
    pause(0.01)

end
end