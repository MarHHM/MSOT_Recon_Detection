function [I] = ImOverlay(im, umx, save, name)
    
    figure;
    im = im - min(im(:));im = im./max(im(:));
    rgb_im(:,:,1) = im;
    rgb_im(:,:,2) = im;
    rgb_im(:,:,3) = im;
    %axes('Units','pixels','Position',[10 10 size(im,2) size(im,1)]);
    imshow(rgb_im);     % Originalbild
    hold on;

    % Binärbild
    umx = umx - min(umx(:)); umx = umx./max(umx(:));
    rgb_umx(:,:,1) = zeros(size(im,1),size(im,2));
    rgb_umx(:,:,2) = umx;
    rgb_umx(:,:,3) = zeros(size(im,1),size(im,2));
    h = imagesc(rgb_umx); 
    set(h,'AlphaData',umx);
    
    if save == 1
        fr = getframe();
        imwrite(fr.cdata,[ name '.tif']);
    end
    I = fr.cdata;
end
