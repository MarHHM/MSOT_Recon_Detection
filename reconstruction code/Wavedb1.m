function y = Wavedb1(x,xsize,trans)
persistent s;
level =6;
wname = 'db4';

if ~trans;  % W*X
    [Y,s] = wavedec2(reshape(x,xsize),level,wname);
else        % W'*X
    if ~exist('s','var')
    [~,s] = wavedec2(reshape(x,xsize),level,wname);
    end
    Y = waverec2(x,s,wname);
    %Y=reshape(Y,xsize);
end
y = Y;