function [meanRMS] = plotSphRms(Y_in,scale,INDB,Y_ref,secDirs,BALLOON)
%PLOTSPHRMS Plotting helper showing SH input signal amplitude RMS.
%   Chris Hold, 2022

    if nargin < 6
        BALLOON = false;
    end
    if nargin < 5
        secDirs = [];  % in deg
    end
    if nargin < 4
        Y_ref = [];
    end
    if nargin < 3
        INDB = false;
    end
    if nargin < 2
        scale = 1;
    end  

    numSH = size(Y_in, 2);
    orderSH = sqrt(numSH) - 1;
    [grid_azi, grid_ele] = ndgrid(linspace(pi, -pi, 120), ...
                                  linspace(pi/2, -pi/2, 60));
    %Y_smp = getRSH(orderSH, [rad2deg(grid_azi(:)), rad2deg(grid_ele(:))] );
    Y_smp = getSH(orderSH, [grid_azi(:), pi/2 - grid_ele(:)], 'real').';
    discrete_in = (4*pi/(numSH)) * Y_in * Y_smp;
    rms_in = rms(discrete_in);

    if ~isempty(Y_ref)
        discrete_ref = (4*pi/(numSH)) * Y_ref * Y_smp;
        rms_ref = rms(discrete_ref);
        rms_plot = rms_in - rms_ref;
    else
        rms_plot = rms_in;
    end
    
    meanRMS = mean(rms_plot);

    if INDB
        rms_plot = db(rms_plot);
        meanRMS = db(meanRMS);
    end
 
    % 2D
    figure;
    p = pcolor(rad2deg(grid_azi), rad2deg(grid_ele),...
                scale * reshape(rms_plot, size(grid_azi)) );
    p.FaceColor = 'interp';
    p.EdgeColor = 'none';
    set(gca, 'XDir','reverse')
    xticks([-180:90:180])
    xlim([-180, 180])
    yticks([-90:45:90])
    ylim([-90, 90])
    xlabel('Azimuth')
    ylabel('Elevation')
    hold on
    line([0, 0], [-90, 90], 'color', 'black', 'lineWidth', 1, 'lineStyle', ':')
    line([-180, 180], [0, 0], 'color', 'black', 'lineWidth',1, 'lineStyle', ':')
    daspect([1 1 1])
    c = colorbar();
    c.Label.String = 'RMS';
    if INDB
        c.Label.String = c.Label.String + " in dB";
    end
    % sectors
    if ~isempty(secDirs)
        secPos = unitSph2cart(deg2rad(secDirs));
        for idxSec = 1:size(secPos, 1)
            scatter(secDirs(idxSec, 1), secDirs(idxSec, 2), 'k+')
            text(secDirs(idxSec, 1) -5, secDirs(idxSec, 2),...
                 num2str(idxSec));
        end
    end


    % 3D
    if BALLOON
    figure;
    [xplot, yplot, zplot] = sph2cart(grid_azi, grid_ele, ...
                                     scale * reshape(rms_plot, ...
                                     size(grid_azi)));
    mesh(xplot, yplot, zplot, reshape(rms_plot, size(grid_azi)), ...
         'FaceAlpha', '0.25', 'FaceColor', 'interp')
    c = colorbar();
    c.Label.String = 'RMS';
    if INDB
        c.Label.String = c.Label.String + " in dB";
    end
    
    hold on
    % axes
    line([-1, 1], [0, 0], [0, 0], 'Color', 'black', 'LineWidth', 1)
    line([0, 0], [-1, 1], [0, 0], 'Color', 'black', 'LineWidth', 1)
    line([0, 0], [0, 0], [-1, 1], 'Color', 'black', 'LineWidth', 1)
    % sectors
    if ~isempty(secDirs)
        secPos = unitSph2cart(deg2rad(secDirs));
        for idxSec = 1:size(secPos, 1)
            plot3([0, secPos(idxSec, 1)], ...
                  [0, secPos(idxSec, 2)], ...
                  [0, secPos(idxSec, 3)], ...
                  'LineStyle', ':', 'Color', '#838383', 'LineWidth', 1)
            text(secPos(idxSec, 1),...
                 secPos(idxSec, 2),...
                 secPos(idxSec, 3),...
                 num2str(idxSec));
        end
    end

    grid on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')

    axis equal
    axis vis3d

    view(-60, 25)
    end

end
