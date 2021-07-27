function [overlay] = overlay(sat, im, background, Mask)
    
    %We use a visualization from green (0% sO2) to red (100% sO2) which is
    %half the total jet scale.
    
    sat(~Mask) = -1; sat = sat - min(sat(:)); sat = sat./2;
%     figure; imagesc(sat)
    rgb_SAT = ind2rgb(uint8(ceil(64*sat)),jet);
    
    % blue channel - remove values that correspond to the outside of the mask
    rgb_SAT_3 = rgb_SAT(:,:,3); rgb_SAT_3(find(rgb_SAT_3 ==0.5625)) = 0; rgb_SAT(:,:,3) = rgb_SAT_3;
    
    %Background
    background = background - min(background(:));
    background = background./max(background(:));
    
    im = im- min(im(:));
    im = im./max(im(:));
    
    %modify the intensity of all channels using the background intenity
    %image - merging of anatomical and functional information
    for i=1:3
        rgb_SAT(:,:,i) = rgb_SAT(:,:,i).*(sqrt(background));
    end
    rgb_BCK = ind2rgb(uint8(ceil(64*im)),gray); 
    
    overlay = rgb_BCK;
    idx = find(Mask);
    
    %combine the background anatomical image with the merged sO2 image
    rgb_SAT = permute(rgb_SAT,[3 1 2]);
    overlay = permute(rgb_BCK,[3 1 2]); overlay(:,idx) = rgb_SAT(:,idx); overlay = permute(overlay,[2 3 1]);

end
