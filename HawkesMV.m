% This code implements the method of Dassios and Zhao to synthetically
% generate a Multivariate Hawkes process with exponentially decaying
% intensity rates
%
% Paper: Exact simulation of Hawkes process with exponentially decaying
% intensity - Angelos Dassios and Hongbiao Zhao
%
% The code has been written to conform to the notation of our paper instead
% of the notation of the Dassios and Zhao paper.  The important
% differences:
%
% Us     Them
% mu     lambda
% delta  (N/A)
% r      delta
% mu_bar a
%
% Code created: 10/09/2013
% Last Modified: 1/06/2016

function [T, T_index, mu_T] = HawkesMV(W,mu_bar,alpha,mu_0,T_horizon)

K = 5000;
D = length(mu_bar);
r= -log(alpha);
if numel(r) == 1
    r = r*ones(D,1);
end
mu_T=mu_0;

T=zeros(1,K+1);
T_index=zeros(size(T));
N=zeros(D,K+1);
k = 1;
fprintf('Generating Data: ')
%disp_space = 5000;
disp_space = floor(.2*T_horizon);
while max(T) < T_horizon
    if k>1 && mod(T(k),disp_space) < mod(T(k-1),disp_space)
        if ceil(T(k)/disp_space)*disp_space ~= T_horizon
            fprintf('%.0f%%, ',100*T(k)/T_horizon)
        else
            fprintf('%.0f%%', 100*T(k)/T_horizon)
        end
    end
    S=zeros(1,D);
    for j=1:D
        D_j=1+r(j)*log(rand)/(mu_T(j)-mu_bar(j));
        S_1=-log(D_j)/r(j);
        S_2=-log(rand)/mu_bar(j);
        
        if D_j>0
            S(j)=min(S_1,S_2);
        elseif D_j<0
            S(j)=S_2;
        end
    end
    
    [W_k, ell]=min(S);
    
    T_index(k+1)=ell;
    T(k+1)=T(k)+W_k;
    
    delta_N=zeros(D,1);
    delta_N(ell)=1;
    N(:,k+1)=N(:,k)+delta_N;
    mu_T=(mu_T-mu_bar).*exp(-r.*(T(k+1)-T(k)))+mu_bar+W(:,ell);
    
    k = k+1;
end
fprintf(', 100%%\n')


T=T(2:k-1);
T_index=T_index(2:k-1);



