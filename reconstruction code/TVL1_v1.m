function x = TVL1_v1(x0,params)
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
alpha = params.lineSearchAlpha;
beta = params.lineSearchBeta;
t0 = params.lineSearchT0;
k = 0;

% copmute g0  = grad(Phi(x))
g0 = wGradient(x,params);
dx = -g0;


% iterations
while(1)

% backtracking line-search

	% pre-calculate values, such that it would be cheap to compute the objective
	% many times for efficient line-search
	[FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx] = preobjective(x, dx, params);
	f0 = objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,x,dx, 0, params);
	t = t0;
        [f1, ERRobj, RMSerr]  =  objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,x,dx, t, params);
	
	lsiter = 0;

	while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 & (lsiter<maxlsiter)
		lsiter = lsiter + 1;
		t = t * beta;
		[f1, ERRobj, RMSerr]  =  objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,x,dx, t, params);
	end

	if lsiter == maxlsiter
		disp('Reached max line search,.... not so good... might have a bug in operators. exiting... ');
		return;
	end

	% control the number of line searches by adapting the initial step search
	if lsiter > 2
		t0 = t0 * beta;
	end 
	
	if lsiter<1
		t0 = t0 / beta;
	end

	x = (x + t*dx);
    
    if params.noneg
    x=params.XFMIN(max(params.XFM(x),0));
    end

	%--------- uncomment for debug purposes ------------------------	
	%disp(sprintf('%d   , obj: %f, RMS: %f, L-S: %d', k,f1,RMSerr,lsiter));

	%---------------------------------------------------------------
	
    %conjugate gradient calculation
    
	[g1,gradObj,gradXFM,gradTV]= wGradient(x,params);
	bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
	g0 = g1;
	dx =  - g1 + bk* dx;
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
    
    DXFMtx = params.TV(shw);
    DXFMtdx = params.TV(shw1);
    
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
grad = params.XFMIN(params.GTV(G));






