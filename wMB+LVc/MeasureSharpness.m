% calc sharpness of brain dataset ROIs using "sharpness()" (part of k-wave - mainly Brenner's gradient used in Subhamoy 2014)

I = imread( 'C:\Users\marwan.muhammad\Desktop\ROI 1 - MB.tif' );
% figure, imagesc(I), colormap(gray), truesize

sh1 = sharpness(I);
sh1_norm = sh1/(size(I,1)*size(I,2));

I = imread( 'C:\Users\marwan.muhammad\Desktop\ROI 1 - wMB+LVc.tif' );
% figure, imagesc(I), colormap(gray), truesize

sh2 = sharpness(I);
sh2_norm = sh2/(size(I,1)*size(I,2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I = imread( 'C:\Users\marwan.muhammad\Desktop\ROI 2 - MB.tif' );
% figure, imagesc(I), colormap(gray), truesize

sh3 = sharpness(I);
sh3_norm = sh3/(size(I,1)*size(I,2));

I = imread( 'C:\Users\marwan.muhammad\Desktop\ROI 2 - wMB+LVc.tif' );
% figure, imagesc(I), colormap(gray), truesize

sh4 = sharpness(I);
sh4_norm = sh4/(size(I,1)*size(I,2));
