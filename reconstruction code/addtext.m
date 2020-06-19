function R=addtext(R,n,CurrRep,CurrSlice,Wavelengths,CurrWav,image_width,c,recon_method,MB_regu,truncated,noneg)
ratio_R=1/max(R(:));
   R=ratio_R*R;
position =  [n (n-12); n n];
text_str{2} = ['R',num2str(CurrRep),'_S',num2str(CurrSlice),'_W',num2str(Wavelengths(CurrWav))];

if strcmp (recon_method,'MB_Tik')
    if noneg
    text_str{1} = ['Tik_NN'];
    else
    text_str{1} = ['Tik'];    
    end
end
    
if strcmp(recon_method,'MB_TVL1')
    if noneg
   text_str{1} = ['L1_NN'];
    else
    text_str{1} = ['L1'];    
    end
end

if strcmp(recon_method,'WaveletPacket')
    if noneg
    text_str{1} = ['WP_NN'];
    else
    text_str{1} = ['WP'];    
    end
end

if strcmp(recon_method,'BackProjection')
    if noneg
    text_str{1} = ['BP_NN'];
    else
    text_str{1} = ['BP'];    
    end
end

R = (insertText(R ,position,text_str,'FontSize',9,'BoxColor','black','TextColor','white','AnchorPoint','RightBottom'));
R=squeeze(R(:,:,1));
R=R/ratio_R;
end