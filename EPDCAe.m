function x = EPDCAe(A,b,lambda,t,x0)
%proposed by Lu, Zhou, Sun (Algorithm 2)
epsilon = 1e-6;
eta = 1e-3;
theta = zeros(200,1);
beta_v = zeros(200,1); 
theta(1) = (1+sqrt(5))/2;
beta_v(1) = 0;
for i = 2:200
    theta(i) = (1+sqrt(1+4*theta(i-1)^2))/2;
    beta_v(i) = 1*(theta(i-1)-1)/theta(i);
end
alpha=1;
c=1;
x=x0;
bx=x;
n=length(x);
s = lambda/t;
hzt = sum(max(x-t,0)+max(0,-x-t));
FO_rho = norm(A*x-b)^2/2 + s*(norm(x,1)-hzt);
gf = A'*(A*bx-b);
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
    funv = Inf;   
    for i=1:Num
        i22 = bitget(i-1,num:-1:1)';
        gx(Is==2) = sign(x(Is==2)).*i22;
        z=max(s*gx-gf+alpha*x+c*bx-s,min(0,s*gx-gf+alpha*x+c*bx+s))/(alpha+c);
        hzt = sum(max(z-t,0)+max(0,-z-t));
        F_rho = norm(A*z-b)^2/2 + s*(norm(z,1)-hzt)+(z-x)'*(z-x)/2;
        if F_rho<funv
            x_new=z;
            funv=F_rho;
        end
    end
    hzt = sum(max(x_new-t,0)+max(0,-x_new-t));
    F_rho = norm(A*x_new-b)^2/2 + s*(norm(x_new,1)-hzt);
    if abs(FO_rho -F_rho)<epsilon
        break
    end    
    FO_rho=F_rho;
    if rem(iter,200)==0
        beta = 0;
    else
        remm = rem(iter,200);
        beta = beta_v(remm);
    end
    bx=x_new+beta*(x_new-x);
    x = x_new;
    gf = A'*(A*bx-b);   
end
