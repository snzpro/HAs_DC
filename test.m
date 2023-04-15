% Solve
% min{1/2||Ax-b||^2+lambda*(||x||_1-h(x))/t}
% by HA, fnHPDCA, pnHPDCA, EPDCAe, NEPDCA and NEPDCAe.
clc
clear all
close all
filename = 'computational result.txt';
pathname = 'C:\code\Code_DC'; 
fid = fopen([pathname filename],'w');
n = 2^13;
m = 2^11;
Kr = 400:405;
num_K = length(Kr);
num_test = 20;
t = 0.11; 
lambda = 0.01;
r = 0.001;
for kk = 1:num_K
        K = Kr(kk);      
        funv = zeros(1,num_test);
        succ = zeros(1,num_test);
        time = zeros(1,num_test);
        for nnt=1:num_test
            fprintf('\n\nExperiment on m = %d, n = %d, K = %d, \t No. test = %d.\n', m,n,K,nnt);
            fprintf('----------------------------------------------------------\n');
            randn('seed',nnt);
            rand('seed',nnt);
            xs = zeros(n,1);
            q = randperm(n);
            xs(q(1:K)) = sign(randn(K,1));
            A = randn(m,n);
            A = orth(A')';
            b = A*xs + r*randn(m,1);
            xinit = zeros(size(xs));
            
            st = tic;
            xs_f = HA(A,b,lambda,t,xinit);
            %xs_f = pnHPDCA(A,b,lambda,t,xinit);
            %xs_f = fnHPDCA(A,b,lambda,t,xinit);
            %xs_f = EPDCAe(A,b,lambda,t,xinit);
            %xs_f = NEPDCA(A,b,lambda,t,xinit);
            %xs_f = NEPDCAe(A,b,lambda,t,xinit);
            time(nnt) = toc(st);
            
            [funv(nnt),succ(nnt)] = funv_succ(A,b,xs_f,xs,t);
            fprintf('funv=%1.4e succ=%1.4e time=%1.4f\n',funv(nnt),succ(nnt),time(nnt));
            fprintf('----------------------------------------------------------\n');
        end
        filename = ['m=' num2str(m) 'n=' num2str(n) 'K=' num2str(K)];
        save([pathname filename]);     
        fprintf(fid,'\n\nAveraged result for m = %d, n = %d, K = %d\n', m,n,K);
        fprintf(fid,'----------------------------------------------------------\n');
        fprintf(fid,'SOLVER:  funv=%1.4e succ=%1.4e time=%1.4e\n',sum(funv)/num_test,sum(succ)/num_test,sum(time)/num_test);
        fprintf(fid,'----------------------------------------------------------\n');
end
fclose(fid);

