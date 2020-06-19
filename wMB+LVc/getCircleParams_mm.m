function [xc, yc, Rc] = getCircleParams_mm(xc_pix, yc_pix, Rc_pix, n, image_width)
% convert from pixel --> mm

xc = (xc_pix-n/2) *image_width/n;
yc = (yc_pix-n/2) *image_width/n;
Rc = Rc_pix *image_width/n;


