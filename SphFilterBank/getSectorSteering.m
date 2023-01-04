function [secDirs] = getSectorSteering(in_nm, method)
%GETSECTORSTEERING Sector steering directions and rotate based on 'method'.
% Inputs:
%    in_nm  : [t x nSH]         - real SH stream, ACN-N3D
%    method : ['front', 'pint', 'sphMUSIC']
%
% Outputs:
%    secDirs : [d x 2]          -  azi, ele in deg
%   Chris Hold, 2022
N_sph = sqrt(size(in_nm, 2)) - 1;
[~, secDirs] = getTdesign(2*N_sph);
secPos = unitSph2cart(secDirs);

if(strcmp(method, 'front'))
    v = [1, 0, 0];
elseif(strcmp(method, 'pint'))
    pint = [in_nm(:,1)' * in_nm(:,4), ...
            in_nm(:,1)' * in_nm(:,2), ...
            in_nm(:,1)' * in_nm(:,3)];
    v = pint ./ norm(pint);
elseif(strcmp(method, 'sphMUSIC'))
    sphCov = in_nm' * in_nm;
    grid_dirs = grid2dirs(1, 2, false);
    [~, est_dirs] = sphMUSIC(sphCov, grid_dirs, 1);
    v = unitSph2cart(est_dirs);
else
    v = secPos(1, :);
end

R = calculateRotationMatrix(secPos(1, :), v);
secPos = secPos * R;  % rotate to 'method'
secDirs = rad2deg(unitCart2sph(secPos));
end