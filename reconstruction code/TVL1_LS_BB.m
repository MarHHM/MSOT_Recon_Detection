function [x,ERR,RN,OF,RMSD] = TVL1_LS_BB(x0,params)
%-----------------------------------------------------------------------
% 
% The function solves the following problem:
%  
% finds the image x that minimizes:
%
% Phi(x) = ||M* x - p||^2 + lambda1*|x|_1 + lambda2*TV(W'*x) 
%
%-------------------------------------------------------------------------
x = x0;


% line search parameters
maxlsiter = params.lineSearchItnlim ;
gradToll = params.gradToll ;
alpha = params.lineSearchAlpha;    beta = params.lineSearchBeta;
%t0 = params.lineSearchT0;
k = 1;
ERR=zeros(1,params.Itnlim);
RN=zeros(1,params.Itnlim);
OF=zeros(1,params.Itnlim);
RMSD=zeros(1,params.Itnlim);
% copmute g0  = grad(Phi(x))
g0 = wGradient(x,params);
aa=(params.XFM(g0));
% norm((aa(:)))
t=(norm((aa(:)))/norm((params.FT*(aa(:)))))^2;
dx = -g0;

% iterations
while(1)

     %k  
%      dx=params.XFMIN(max(params.XFM(dx),0));
%      [FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx] = preobjective(x, dx, params);
%      [f1, ERRobj, RMSerr]  =  objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,x,dx, t, params);
	
     x1 = (x + t*dx);
    
    if params.noneg
    x1=params.XFMIN(max(params.XFM(x1),0));
    end
    
        Recon = params.XFM(x1);
   Recon=Recon/max(max(Recon));
        
	%--------- uncomment for debug purposes ------------------------	
   	%disp(sprintf('%d   , obj: %f, RMS: %f, L-S: %d', k,f1,RMSerr,lsiter));

	%---------------------------------------------------------------
    	%shw=params.XFM(x1);%shw=shw/max(max(shw));
        %ERR(k+1)=sqrt(sum(sum(abs(shw-params.ori).^2))/length(params.ori).^2);
         ERR(k+1)=norm(x1-x)/norm(x1);
         %RN(k+1)=ERRobj;
        % OF(k+1)=f1;
    
%conjugate gradient calculation
    
	[g1,gradObj,gradXFM,gradTV]= wGradient(x1,params);

    t=((norm(x1-x))^2)/((x1-x)*((g1-g0)'));
	
	g0 = g1;
    x = x1;
	dx =  - g1 ;
    k = k + 1;
	
	%TODO: need to "think" of a "better" stopping criteria ;-)
	if (k > params.Itnlim) | (norm(dx(:)) < gradToll) 
		break;
	end

end

return;


function [FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx] = preobjective(x, dx, params)

% precalculates transforms to make line search cheap
shw=params.XFM(x);
% size(shw)
shw1=params.XFM(dx);
% size(shw1)
FTXFMtx = params.FT*(shw(:));
FTXFMtdx = params.FT*(shw1(:));

if params.TVWeight
    %DXFMtx = params.TV*(shw);
    DXFMtx = params.TV(shw);
    %DXFMtx=reshape(DXFMtx,size(x,1),size(x,2));
    %DXFMtdx = params.TV*(shw1);
    DXFMtdx = params.TV(shw1);
    %DXFMtdx=reshape(DXFMtdx,size(x,1),size(x,2));
else
    DXFMtx = 0;
    DXFMtdx = 0;
end

function [res, obj, RMS] = objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx, x,dx,t, params)
%calculated the objective function

p = params.pNorm;

obj = FTXFMtx + t*FTXFMtdx - params.data;
obj = obj(:)'*obj(:);

if params.TVWeight
    w = DXFMtx(:) + t*DXFMtdx(:);
    TV = (w.*conj(w)+params.l1Smooth).^(p/2); 
else
    TV = 0;
end

if params.xfmWeight
   w = x(:) + t*dx(:); 
   XFM = (w.*conj(w)+params.l1Smooth).^(p/2);
else
    XFM=0;
end

TV = sum(TV.*params.TVWeight(:));
XFM = sum(XFM.*params.xfmWeight(:));
RMS = sqrt(obj/sum(abs(params.data(:))>0));

res = obj + (TV) + (XFM) ;

function [grad,gradObj,gradXFM,gradTV] = wGradient(x,params)

gradXFM = 0;
gradTV = 0;

gradObj = gOBJ(x,params);
if params.xfmWeight
gradXFM = gXFM(x,params);
end
if params.TVWeight
gradTV = gTV(x,params);
end

grad = (gradObj +  params.xfmWeight.*gradXFM + params.TVWeight.*gradTV);

function gradObj = gOBJ(x,params)
% computes the gradient of the data consistency

	shw=params.XFM(x);
        gradObj = params.XFMIN(params.FT'*(params.FT*(shw(:)) -params.data));   

gradObj = 2*gradObj ;


function grad = gXFM(x,params)
% compute gradient of the L1 transform operator

p = params.pNorm;

grad = p*x.*(x.*conj(x)+params.l1Smooth).^(p/2-1);

function grad = gTV(x,params)
% compute gradient of TV operator

p = params.pNorm;
Dx =params.TV(params.XFM(x));
G = p*Dx.*(Dx.*conj(Dx) + params.l1Smooth).^(p/2-1);
%size(params.TV'*G)
grad = params.XFMIN(params.GTV(G));
