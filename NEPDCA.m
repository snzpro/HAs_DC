function x = NEPDCA(A,b,lambda,t,x0)
%proposed by Lu, Zhou (Algorithm 1)
epsilon = 1e-6;
eta=1e-3;
x=x0;
n=length(x);
s = lambda/t;
hzt = sum(max(x-t,0)+max(0,-x-t));
FO_rho = norm(A*x-b)^2/2 + s*(norm(x,1)-hzt);
M = 5;
Fv = zeros(M,1);
alpha = 1;
gamma = 1e-4;
alpha_min = 1e-8;
alpha_max = 1e+8;
gf = A'*(A*x-b);
for iter=1:1000  
    I=zeros(n,3);
    dxt=max(x-t,0)+max(0,-x-t)-eta;
    I(x-t>=dxt,1) = 1;
    I(0>=dxt,2) = 1;
    I(-x-t>=dxt,3) = 1;
    Is = sum(I,2);
    num = sum(Is==2);
    if num>20 
        break
    end
    Num = 2^num;
    Gx = [ones(n,1),zeros(n,1),-ones(n,1)]; 
    nzg = sum(I.*Gx,2);
    gx = zeros(n,1);
    gx(Is==1) = nzg(Is==1);
    
    remm = rem(iter,M);
    Fv(remm+1) = FO_rho;  
    while 1    
        funv = Inf;
        II=zeros(Num,2);
        II(:,1) = norm(A*x-b)^2/2+s*norm(x,1);
        for i=1:Num
            i22 = bitget(i-1,num:-1:1)';
            gx(Is==2) = sign(x(Is==2)).*i22;
            z=max(s*gx-gf+alpha*x-s,min(0,s*gx-gf+alpha*x+s))/alpha;         
            II(i,2) = norm(z-x)^2;
            hzt = sum(max(z-t,0)+max(0,-z-t));
            F_rho = norm(A*z-b)^2/2 + s*(norm(z,1)-hzt) + gamma/2*II(i,2);   
            II(i,1) = II(i,1) - s*sum(gx.*x-gx.^2.*t);%f_s(x^k)+f_n(x^k)-psi_i(x^k)
            if F_rho<funv
                x_new=z;
                funv=F_rho;
            end
        end
        Mum=0;
        hzt = sum(max(x_new-t,0)+max(0,-x_new-t));
        F_rho = norm(A*x_new-b)^2/2 + s*(norm(x_new,1)-hzt);
        for i=1:Num 
            if F_rho<=max(max(Fv),II(i,1))-gamma*norm(x_new-x)^2/2-gamma*II(i,2)/2
                Mum=Mum+1;
            end
        end 
        if Mum==Num
            break
        else
            alpha=2*alpha;
        end  
    end 
    TERM = FO_rho -F_rho;
    if abs(TERM) < epsilon
        break
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
