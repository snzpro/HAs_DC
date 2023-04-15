function x = fnHPDCA(A,b,lambda,t,x0)
%Algorithm 2
epsilon = 1e-4;
eta = 1e-1;
x=x0;
n=length(x);
s = lambda/t;
hzt = sum(max(x-t,0)+max(0,-x-t));
FO_rho = norm(A*x-b)^2/2 + s*(norm(x,1)-hzt);
M = 5;
Fv = zeros(M,1);
alpha = 1;
gamma = 1e-3;
alpha_min = 1e-8;
alpha_max = 1e+8;
gf = A'*(A*x-b);
rho=1;
for iter=1:1000     
    beta=alpha;
    dxt=max(x-t,0)+max(0,-x-t);
    II=zeros(n,1);
    II(x-t==dxt) = 1;
    II(0==dxt) = 0;
    II(-x-t==dxt) = -1;
    gx = II;   
    ii = rem(iter,M);
    Fv(ii+1) = FO_rho;  
    while 1        
        z=max(s*gx-gf+alpha*x-s,min(0,s*gx-gf+alpha*x+s))/alpha;
        hzt = sum(max(z-t,0)+max(0,-z-t));
        F_rho = norm(A*z-b)^2/2 + s*(norm(z,1)-hzt);
        if F_rho<max(Fv)-gamma*norm(z-x)^2
            break
        end
        alpha=2*alpha;    
    end           
    TERM = max(Fv) -F_rho;   
    if TERM >= rho*epsilon
        x_new=z;
    else
        I=zeros(n,3);
        dxt=max(x-t,0)+max(0,-x-t)-eta;
        I(x-t>=dxt,1) = 1;
        I(0>=dxt,2) = 1;
        I(-x-t>=dxt,3) = 1;
        Is = sum(I,2);
        num = sum(Is==2);
        Num = 2^num;
        while 1 
            lin_s=1;
            funv = Inf;
            for i=1:Num
                i22 = bitget(i-1,num:-1:1)';
                gx(Is==2) = sign(x(Is==2)).*i22;
                z=max(s*gx-gf+alpha*x-s,min(0,s*gx-gf+alpha*x+s))/alpha; 
                hzt = sum(max(z-t,0)+max(0,-z-t));
                F_rho = norm(A*z-b)^2/2 + s*(norm(z,1)-hzt);
                if F_rho<max(Fv)-max(rho*epsilon,gamma*norm(z-x)^2)
                    lin_s = 0;
                    x_new = z;
                    break
                end
                l_rho = (gf-s*gx)'*(z-x) + 1/2*norm(z-x)^2 + s*norm(z,1) - s*sum(gx.*z-gx.^2.*t);
                if l_rho<funv
                    x_new = z;
                    funv = l_rho;
                end   
            end
            if lin_s==0
                break
            else
                hzt = sum(max(x_new-t,0)+max(0,-x_new-t));
                F_rho = norm(A*x_new-b)^2/2 + s*(norm(x_new,1)-hzt);
                if F_rho<max(Fv)-gamma*norm(x_new-x)^2
                    break
                else
                    alpha=2*alpha;
                end
            end
        end  
        TERM = max(Fv) -F_rho;
        if TERM < rho*epsilon
            if epsilon<=1e-6
                break
            end
            epsilon = 0.1*epsilon;
        end    
    end
    gfn = A'*(A*x_new-b);  
    if norm(x-x_new)==0
        alpha = alpha_max;
    else
        alpha = abs((gfn-gf)'*(x_new-x))/norm(x_new-x)^2;
        alpha = max(alpha_min,min(alpha,alpha_max));
    end    
    gf = gfn;  
    FO_rho=F_rho;
    x = x_new;
end
