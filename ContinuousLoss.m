function out = ContinuousLoss(T,T_index,alpha,mu_bar,W)
% CONTINUOUSLOSS   Calculates the continuous time loss for a Hawkes
% Process with exponential decay
%
% This function takes in event times and nodes, as well with parameters of
% an exponential decay influence function and calculates the continuous
% time, exact value of the negative log likelihood

p=length(W);
N=length(T);
T_horizon=ceil(T(end));
mu_true=zeros(p,1); %The true rate of the process at the N event times, T_horizon
loss_true=zeros(1,N+1);
e_kn=zeros(p,1);

for n=1:N+1
    if n==1
        mu_true=mu_bar;
        delta=T(1);
        e_kn(T_index(n))=1;
        loss_true(n)=delta*sum(mu_bar)-e_kn(:,n)'*log(mu_bar);
    elseif n==N+1
        delta=T_horizon-T(n-1);
        mu_true_prev=mu_true;
        mu_true=alpha^delta*mu_true_prev+W*e_kn*(alpha^delta)+(1-alpha^delta)*mu_bar;
        a_delta=(alpha^delta-1)/log(alpha);
        loss_true(n)=a_delta*sum(mu_true_prev+W*e_kn)+(delta-a_delta)*sum(mu_bar);
    else
        delta=T(n)-T(n-1);
        mu_true_prev=mu_true;
        mu_true=alpha^delta*(mu_true+W*e_kn)+(1-alpha^delta)*mu_bar;
        e_prev=e_kn;
        e_kn=zeros(p,1);
        e_kn(T_index(n))=1;
        a_delta=(alpha^delta-1)/log(alpha);
        loss_true(n)=a_delta*sum(mu_true_prev+W*e_prev)+(delta-a_delta)*sum(mu_bar)-e_kn'*log(mu_true);
    end
end

out=loss_true;