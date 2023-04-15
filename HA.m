function x = HA(A,b,lambda,t,x0)
%Algorithm 1
epsilon = 1e-4;
eta = 1.0e-1;
x=x0;
n=length(x);
s = lambda/t;
x_new = x;
alpha = 1;
hzt = sum(max(x-t,0)+max(0,-x-t));
FO_rho = norm(A*x-b)^2/2 + s*(norm(x,1)-hzt);
for iter=1:2000   
    dxt = max(x-t,0) + max(0,-x-t);
    II = zeros(n,1);
    II(x-t==dxt) = 1;
    II(0==dxt) = 0;
    II(-x-t==dxt) = -1;
    gf = A'*(A*x-b);
    gx = II;   
    z=max(s*gx-gf+alpha*x-s,min(0,s*gx-gf+alpha*x+s))/alpha; 
    hzt = sum(max(z-t,0)+max(0,-z-t));
    F_rho = norm(A*z-b)^2/2 + s*(norm(z,1)-hzt);
    TERM = FO_rho -F_rho;    
    if TERM >= epsilon       
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
        funv = Inf;
        for i=1:Num
            i22 = bitget(i-1,num:-1:1)';
            gx(Is==2) = sign(x(Is==2)).*i22;
            gf = A'*(A*x-b);
            z=max(s*gx-gf+alpha*x-s,min(0,s*gx-gf+alpha*x+s))/alpha;
            hzt = sum(max(z-t,0)+max(0,-z-t));
            F_rho = norm(A*z-b)^2/2 + s*(norm(z,1)-hzt);
            TERM = FO_rho -F_rho;                      
            if TERM>=epsilon
                x_new=z;
                break
            end           
            l_rho = (gf-s*gx)'*(z-x) + 1/2*norm(z-x)^2 + s*norm(z,1) - s*sum(gx.*z-gx.^2.*t);     
            if l_rho<funv
                x_new=z;
                funv=l_rho;
            end          
        end       
        hzt = sum(max(x_new-t,0)+max(0,-x_new-t));
        F_rho = norm(A*x_new-b)^2/2 + s*(norm(x_new,1)-hzt);
        TERM = FO_rho - F_rho;       
        if TERM < epsilon
            if epsilon<=1e-6
                break
            end
            epsilon = 0.1*epsilon;
        end   
    end   
    FO_rho = F_rho;
    x = x_new;
end
