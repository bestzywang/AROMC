function [Out_X, Y1_RMSE,Y2_RMSE,peaksnr,U1,V1,NRE,PMD] = AROMC(M, M_Omega, array_Omega,maxiter)
%AUTOMATIC_RANK Summary of this function goes here
%   Detailed explanation goes here
% M - true data
% M_Omega - observed data with missing entries
% array_Omega - observed index, which is a matrix, '1' means observed
% data,while '0' refers to missing data. For example, per =0.6,
% array_Omega = binornd( 1, per, [ m, n ] ), which implies 60% observations
% and 40% missing entries,
% M_Omega = M.*array_Omega;
% maxiter - maximum iteration number

[r,c] = size(M_Omega);
Y1_RMSE = [];
Y2_RMSE = [];
NRE(1:maxiter)=zeros;
X = zeros(r, c);
peaksnr = [];
U = randn(r,1);
V = randn(1,c);
M_2 = M_Omega;
U1 = [];

PMD = [];
SPMD = [];
for iter = 1 : maxiter
    M_1 = M_2;
    M_Omega = M_1;
    rak = iter;
    in2 = 0;

    while (in2 <2)
%         a = norm((M_Omega - U * V).*array_Omega,'fro')^(2)/norm(M_Omega.*array_Omega,'fro')^(2);
        a = norm((M_Omega - U * V).*array_Omega,'fro')^(2)/norm(M_Omega.*array_Omega,'fro')^(2);               
        
        if iter>1
            M_1= M_Omega - U1*V(1:iter-1,:);
        end
        for i = 1:r
            U_bf = U(i,end);
%             U_bf = 0;
            col = find(array_Omega(i,:) == 1);
            V_I = V(iter,col);
            b_I = M_1(i,col);
%             U(i,iter) = (2*b_I*V_I'+balta*U_bf)/(balta+2*V_I(end,:)*V_I(end,:)');
             U(i,iter) = b_I*V_I'/(V_I(end,:)*V_I(end,:)');
        end  
        
        for j = 1:c
%             V_bf = V(:,j);
%             V_bf = 0;
            row = find(array_Omega(:,j) == 1);
            U_I =  U(row,:);
            b_I = M_Omega(row,j);
%             V(:,j) = inv(2*U_I'*U_I+BALTA)*(2*U_I'*b_I+BALTA*V_bf);
%             V(:,j) = inv(2*U_I'*U_I)*(2*U_I'*b_I);
            V(:,j) = pinv(U_I)*b_I;
            
        end
        
%         b = norm((M_Omega - U * V).*array_Omega,'fro')^(2)/norm(M_Omega.*array_Omega,'fro')^(2);
        b = norm((M_Omega - U * V).*array_Omega,'fro')^2/norm(M_Omega.*array_Omega,'fro')^2;
        if a-b < 0.000001
            in2 = in2 + 1;
        end
    end
    U1 = U;
    V1 = V;
    U = [U1 randn(r,1)];
%     V = randn(iter+1,c);
    V = [V1;randn(1,c)];
    X = U1*V1;
    Y = X;
    iter
    Rank_X=rank(X)
    peaksnr = [peaksnr psnr(Y, M)];
    Y1_RMSE = [Y1_RMSE norm((M-X),'fro')/sqrt(r*c)];
    Y2_RMSE = [Y2_RMSE norm((M_Omega-X).*array_Omega,'fro')/sqrt(r*c)];
%     Y2_RMSE = [Y2_RMSE norm((M_Omega-X).*array_Omega,'fro')/norm((M_Omega).*array_Omega,'fro')];%/sqrt(sum(sum(array_Omega)))];%/norm((M_Omega).*array_Omega,'fro');%/sqrt(sum(sum(array_Omega)))
     NRE(iter) = (norm((M_Omega).*array_Omega,'fro'))^2/(norm((M_Omega-X).*array_Omega,'fro'))^2;
     
%      SPMD = [SPMD norm(U1(:,end)*V1(end,:),'fro')]
    if iter>1
        
        PMD = [PMD (norm(U1(:,end-1)*V1(end-1,:),'fro')^2 - norm(U1(:,end)*V1(end,:),'fro')^2)/norm(U1(:,end-1)*V1(end-1,:),'fro')^2]
        if iter>2
%             SPMD = [SPMD abs((abs(PMD(end))+abs(PMD(end-1)))/2)]
            if abs((abs(PMD(end))+abs(PMD(end-1)))/2)<0.5
%             if PMD(end)<xi
                break
            end
        end
    end

end
Out_X = X;
end


