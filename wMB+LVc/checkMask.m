function [xc_mm, yc_mm, Rc_mm, logicalMask ] = checkMask( LOAD_MASK_DIMS, mask , sigMat_pathName , Im_anatomy, im_w )

if LOAD_MASK_DIMS
    %   [xc, yc, Rc, logicalMask] = loadMaskParams( sigMat_pathName, mask );
    %TO DO: save others masks in a structure form, then delete 'loadMaskParams()'!!
    disp(['mask ' mask ' found on disk & loaded into workspace..'])
    temp = load([sigMat_pathName '\masks\mask_' mask]);
    logicalMask = temp.logicalMask;
    xc_mm = temp.xc_mm;
    yc_mm = temp.yc_mm;
    Rc_mm = temp.Rc_mm;
    
else              % draw image & choose ROI to get mask params and save them
    figure('name', ['please draw mask ' mask]), imagesc(Im_anatomy), title('anatomy'), colormap(bone), axis image off;
    %%% draw an ellipse to get approx. mask
    %   h = imellipse;
    %   ellipseVerts = wait(h);
    %   logicalMask = h.createMask;
    %   delete(h);
    %   [xc_pix, yc_pix, Rc_pix, ~] = circfit( ellipseVerts(:,1), ellipseVerts(:,2));     % fit drawn ellipse vertices to a circle (as we use a cicular mask for the weighting)
    h = drawcircle;
    logicalMask = h.createMask;
    xc_pix = h.Center(1); 
    yc_pix = h.Center(2);
    Rc_pix = h.Radius;
    delete(h);
    hold on; fnplt(rsmak('circle', Rc_pix, [xc_pix yc_pix]));                  % overlay mask on anatomy to confirm circle params
    
    [xc_mm, yc_mm, Rc_mm] = getCircleParams_mm( xc_pix, yc_pix, Rc_pix, size(Im_anatomy,1), im_w );    % convert fitted circle params to mm
    if ~exist([sigMat_pathName '\masks'], 'dir')
        mkdir([sigMat_pathName '\masks'])
    end
    save([sigMat_pathName '\masks\mask_' mask '.mat'],...
        'logicalMask', 'xc_mm', 'yc_mm', 'Rc_mm', 'im_w'); % next time change 'LOAD_CIRCLE_DIMS' in tryLuisWeighting.m to true
end

