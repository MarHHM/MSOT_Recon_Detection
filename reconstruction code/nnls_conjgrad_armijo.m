function [x, f0, t] = nnls_conjgrad_armijo(A,b,x,zetaRel,nIter,cgIter)
%performing NNLS with the projected conjugate gradient method
%Input: x--initial guess of the solution (usually set as vector 0)
%       zetaRel--relative infinity norm of the projected gradient,used as stopping criteria (e.g. 0.01)
%       nIter--max. number of outer iteration
%       cgIter--max. number of inner iteration (usually set below 5, depending on the problem size)
%Output: x--solution
%        f0--vector containing the value ||Ax-b||^2 after each outer iteration
%        t-- vector containing the execution time after each outer iteration
%terminates if nIter or zetaRel is reached or t>10 min.
tic
if ~exist('cgIter','var')
    cgIter = 5;
end
if ~exist('nIter','var')
    nIter = 1000;
end
r = b-A*x;
g = A'*(-r);
gOld = g;
f0 = r'*r;
t = toc;
sigma = 0.1;
tau = 0.5;

s = 1;
p = -g;
gPrev = g;

for iter=1:nIter
    I = x>0 | g<0;
    
    gProj = zeros(size(g));
    gProj(I) = g(I);
    if iter == 1
        zeta = max(abs(gProj));
    else if max(abs(gProj)) < zeta*zetaRel
            break
        end
    end
    
    p(~I) = 0;
    d = zeros(size(x));
    
    
    for i = 1:cgIter
        
        beta = max(0, (-g(I)'*(-g(I)+gOld(I)))/(gOld(I)'*gOld(I)));
        p(I) = -g(I) + beta*p(I);
        
        ap = A*p;
        alpha = -g'*p/(ap'*ap);
        
        d(I) = d(I) + alpha * p(I);
        if i == cgIter
            break;
        else
            r = r - alpha*ap;
            gOld = g;
            g=A'*(-r);
        end
        
    end
    
    xNew = x + s*d;
    nbz = sum(xNew<0);
    
    if nbz == 0
        r = r - alpha*ap;
        x = xNew;
        
        gOld = g;
        g = A'*(-r);
        f0(iter+1) = r'*r;
        t(iter+1) = toc;
        continue
    end
    
    xNew = max(xNew,0);
    r = (b-A*xNew);
    
    fun = r'*r;
    It = x>0 | d>0;
    d(~It) = 0;
    gLin = d'*gPrev;
    
    if fun - f0(iter) > sigma*s*gLin
        funPrev = fun;
        xPrev = xNew;
        rPrev = r;
        s = s*tau;
        while true
            xNew = max((x + s*d), 0);
            r = (b-A*xNew);
            fun = r'*r;
            
            if fun - f0(iter) <= sigma*s*gLin
                if funPrev < fun
                    xNew = xPrev;
                    r = rPrev;
                    fun = funPrev;
                end
                break
                
            end
            funPrev = fun;
            xPrev = xNew;
            rPrev = r;
            s = s*tau;
        end
    end
    %     else
    %         funPrev = fun;
    %         yPrev = y;
    %         rPrev = r;
    %         s = s/beta;
    %         while true
    %             y = y0 + s*d;
    %             y = max(y,0);
    %             r = (b-A(:,I)*y);
    %             fun = r'*r;
    %
    %             if fun - f0(iter) > sigma*s*gLin
    %                 if funPrev < fun
    %                     y = yPrev;
    %                     r = rPrev;
    %                     fun = funPrev;
    %                 end
    %                 break
    %
    %             end
    %             funPrev = fun;
    %             yPrev = y;
    %             rPrev = r;
    %             s = s/beta;
    %         end
    %     end
    
    x = xNew;
    gPrev = g;
    g=A'*(-r);
    f0(iter+1) = fun;
    t(iter+1) = toc;
    
    %terminates after the max. execution time of 10 min.
    if t(iter+1) > 600
        break
    end
end
end

