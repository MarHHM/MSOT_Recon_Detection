function [I] = ImOverlayR(im, umx, save, name)
    
    figure;
    im = im - min(im(:));im = im./max(im(:));
    rgb_im(:,:,1) = im;
    rgb_im(:,:,2) = im;
    rgb_im(:,:,3) = im;
    %axes('Units','pixels','Position',[10 10 size(im,2) size(im,1)]);
%     imshow(rgb_im);     % Originalbild
%     hold on;

    % Bin�rbild
    umx = umx - min(umx(:)); 
    if (max(umx(:))>0)
        umx = umx./max(umx(:));
    end
    rgb_umx(:,:,1) = zeros(size(im,1),size(im,2));
    rgb_umx(:,:,2) = abs(umx);
    rgb_umx(:,:,3) = zeros(size(im,1),size(im,2));
    
    I = rgb_im + rgb_umx;
%     I = imfuse(rgb_im,rgb_umx);
%     h = imagesc(rgb_umx); 
%     set(h,'AlphaData',umx);
%     
%     if save == 1
%         fr = getframe();
%         imwrite(fr.cdata,[ name '.tif']);
%     end
%     I = fr.cdata;
end

