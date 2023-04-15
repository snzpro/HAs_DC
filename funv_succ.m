function [funv,suc] = funv_succ(A,b,x,xs,t)
lambda = 0.01;
hzt = sum(max(x-t,0)+max(0,-x-t));
funv = norm(A*x-b)^2+lambda/t*(norm(x,1)-hzt);
%suc = norm(x-xs,inf)<1e-2;
suc = norm(x-xs)/norm(xs)<1e-2;
end

