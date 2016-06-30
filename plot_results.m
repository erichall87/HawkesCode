load mainSynthResults_temp
savefigs = false;

%Set up filter to do moving average over 'zeta' time points
%(length zeta*delta time window)
zeta=5000;
filt=ones(1,zeta)/zeta;
filt=repmat(fft(filt,N_delta),no_iters,1);

%Calculate moving average losses
ma_loss_DMD=ifft(fft(loss_DMD,[],2).*filt,[],2);
ma_loss_DMD=ma_loss_DMD(:,zeta:end);
ma_loss_DMD_W0=ifft(fft(loss_DMD_W0,[],2).*filt,[],2);
ma_loss_DMD_W0=ma_loss_DMD_W0(:,zeta:end);
ma_loss_DMD_exp=ifft(fft(loss_DMD_exp,[],2).*filt,[],2);
ma_loss_DMD_exp=ma_loss_DMD_exp(:,zeta:end);
ma_loss_OGD=ifft(fft(loss_OGD,[],2).*filt,[],2);
ma_loss_OGD=ma_loss_OGD(:,zeta:end);

ma_loss_true = ifft(fft(loss_true,[],2).*filt,[],2);
ma_loss_true = ma_loss_true(:,zeta:end);

ma_loss_DMD_beta=ifft(fft(loss_DMD_beta,[],2).*filt,[],2);
ma_loss_DMD_beta=ma_loss_DMD_beta(:,zeta:end);
ma_loss_DMD_W0_beta=ifft(fft(loss_DMD_W0_beta,[],2).*filt,[],2);
ma_loss_DMD_W0_beta=ma_loss_DMD_W0_beta(:,zeta:end);
ma_loss_DMD_exp_beta=ifft(fft(loss_DMD_exp_beta,[],2).*filt,[],2);
ma_loss_DMD_exp_beta=ma_loss_DMD_exp_beta(:,zeta:end);
ma_loss_OGD_beta=ifft(fft(loss_OGD_beta,[],2).*filt,[],2);
ma_loss_OGD_beta=ma_loss_OGD_beta(:,zeta:end);

J=1000;

figure(1)
plot(delta*(zeta:J:N_delta),mean(ma_loss_DMD_W0(:,1:J:end),1),...
    delta*(zeta:J:N_delta), mean(ma_loss_DMD(:,1:J:end),1),...
    delta*(zeta:J:N_delta),mean(ma_loss_DMD_exp(:,1:J:end),1),...
    delta*(zeta:J:N_delta),mean(ma_loss_OGD(:,1:J:end),1),'r--',...
    delta*(zeta:J:N_delta),mean(ma_loss_true(:,1:J:end),1),'k--','linewidth',3)
curr_leg=legend('Algorithm 1:W=0','Alg 1:True W','Algorithm 2','OGD','True Rate','location','east');
set(curr_leg,'fontsize',18);
xlabel('Time','fontsize',18)
ylabel('Loss','fontsize',18)
title('Moving Average Loss','fontsize',16)
hold off
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
if savefigs
    print(fig,'MALossAlg2','-dpdf')
    savefig('MALossAlg2.fig')
end

figure(2)
plot(delta*(zeta:J:N_delta),median(ma_loss_OGD(:,1:J:end) - ma_loss_DMD_exp(:,1:J:end)),'b-',...
    delta*(zeta:J:N_delta),prctile(ma_loss_OGD(:,1:J:end) - ma_loss_DMD_exp(:,1:J:end),95),'k-',...
    delta*(zeta:J:N_delta),prctile(ma_loss_OGD(:,1:J:end) - ma_loss_DMD_exp(:,1:J:end),5),'k-','linewidth',3)
xlabel('Time','fontsize',18)
ylabel('Loss','fontsize',18)
curr_leg = legend('Median','95th Percentile','5th Percentile','location','southeast');
set(curr_leg,'fontsize',18);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
if savefigs
    print(fig,'MADiffAlg2','-dpdf')
    savefig('MADiffAlg2.fig')
end

figure(3)
plot(delta*(zeta:J:N_delta),mean(ma_loss_DMD_W0_beta(:,1:J:end),1),...
    delta*(zeta:J:N_delta),mean(ma_loss_DMD_beta(:,1:J:end),1),...
    delta*(zeta:J:N_delta),mean(ma_loss_DMD_exp_beta(:,1:J:end),1),...
    delta*(zeta:J:N_delta),mean(ma_loss_OGD_beta(:,1:J:end),1),'r--',...
    delta*(zeta:J:N_delta),mean(ma_loss_true(:,1:J:end)),'k--','linewidth',3)
curr_leg = legend('Algorithm 1:W=0','Algorithm 1:True W','Algorithm 2','OGD','True Rate','location','East');
set(curr_leg,'fontsize',18)
xlabel('Time','fontsize',18)
ylabel('Loss','fontsize',18)
title('Moving Average Loss, Mismatch','fontsize',16)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
if savefigs
    print(fig,'MALossAlg2Mismatch','-dpdf')
    savefig('MALossAlg2Mismatch.fig')
end

figure(4)
plot(delta*(zeta:J:N_delta),median(ma_loss_OGD_beta(:,1:J:end) - ma_loss_DMD_exp_beta(:,1:J:end)),'b-',...
    delta*(zeta:J:N_delta),prctile(ma_loss_OGD_beta(:,1:J:end) - ma_loss_DMD_exp_beta(:,1:J:end),95),'k-',...
    delta*(zeta:J:N_delta),prctile(ma_loss_OGD_beta(:,1:J:end) - ma_loss_DMD_exp_beta(:,1:J:end),5),'k-','linewidth',3)
xlabel('Time','fontsize',18)
ylabel('Loss','fontsize',18)
curr_leg = legend('Median','95th Percentile','5th Percentile','location','southeast');
set(curr_leg,'fontsize',18);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
if savefigs
    print(fig,'MADiffAlg2Mismatch','-dpdf')
    savefig('MADiffAlg2Mismatch.fig')
end

figure(5)
imagesc(W,[0 1*max(W(:))])
axis off; axis square; 
title('True W Matrix','fontsize',18)
colormap hot; colorbar;
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
if savefigs
    print(fig,'TrueW','-dpdf')
    savefig('TrueW.fig')
end


figure(6)
imagesc(W_DMD,[0 1*max(W_DMD(:))])
axis off; axis square; 
title('Alg 2 Estimate of W Matrix','fontsize',16)
colormap hot; colorbar;
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
if savefigs
    print(fig,'Alg2Estimate','-dpdf')
    savefig('Alg2Estimate.fig')
end

figure(7)
imagesc(W_OGD,[0 1*max(W_DMD(:))])
axis off; axis square; 
title('OGD Estimate of W Matrix','fontsize',16)
colormap hot; colorbar;
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
if savefigs
    print(fig,'OGDEstimate','-dpdf')
    savefig('OGDEstimate.fig')
end

figure(8)
imagesc(W_DMD_beta,[0 1*max(W_DMD(:))])
axis off; axis square; 
title('Alg 2 Estimate of W Matrix, mismatch','fontsize',16)
colormap hot; colorbar;
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
if savefigs
    print(fig,'Alg2EstimateMismatch','-dpdf')
    savefig('Alg2EstimateMismatch.fig')
end

figure(9)
imagesc(W_OGD_beta,[0 1*max(W_DMD(:))])
axis off; axis square; 
title('OGD Estimate of W Matrix, mismatch','fontsize',16)
colormap hot; colorbar;
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
if savefigs
    print(fig,'OGDEstimateMismatch','-dpdf')
    savefig('OGDEstimateMismatch.fig')
end

figure(10)
plot(delta*((0:no_saves)*save_iters),mean(batch_loss_DMD),'b-',...
    delta*((0:no_saves)*save_iters),mean(batch_loss_OGD),'r--',...
    delta*((0:no_saves)*save_iters),mean(batch_loss_DMD_beta),'k-',...
    delta*((0:no_saves)*save_iters),mean(batch_loss_OGD_beta),'c-',...
    ...%[0 T_horizon],repmat(sum(loss_cont_true),1,2),
    'linewidth',3)
hold on
errorbar(delta*((0:no_saves)*save_iters),mean(batch_loss_DMD),std(batch_loss_DMD),'b.')
errorbar(delta*((0:no_saves)*save_iters),mean(batch_loss_OGD),std(batch_loss_OGD),'r.')
errorbar(delta*((0:no_saves)*save_iters),mean(batch_loss_DMD_beta),std(batch_loss_DMD_beta),'k.')
errorbar(delta*((0:no_saves)*save_iters),mean(batch_loss_OGD_beta),std(batch_loss_OGD_beta),'c.')
hold off
xlabel('Time','fontsize',18)
ylabel('Loss','fontsize',18)
curr_leg = legend('Alg 2','OGD','Alg 2, Mismatch','OGD, Mismatch',0);
set(curr_leg,'fontsize',18)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
if savefigs
    print(fig,'BatchLoss','-dpdf')
    savefig('BatchLoss.fig')
end

no_points = 500;
thresh = linspace(0,max([W_final_DMD(:);W_final_DMD_beta(:); W_final_OGD(:);W_final_OGD_beta(:)]),no_points);

P_d_DMD=zeros(1,no_points);
P_f_DMD=zeros(1,no_points);
P_d_OGD = zeros(1,no_points);
P_f_OGD = zeros(1,no_points);
P_d_DMD_beta=zeros(1,no_points);
P_f_DMD_beta=zeros(1,no_points);
P_d_OGD_beta=zeros(1,no_points);
P_f_OGD_beta=zeros(1,no_points);

W_pos = W_true > 0;
W_sum = sum(W_pos(:));
W_neg = sum(1 - W_pos(:));
for t=1:length(thresh);
    W_DMD_temp = W_final_DMD > thresh(t);
    W_DMD_beta_temp = W_final_DMD_beta > thresh(t);
    W_OGD_temp = W_final_OGD > thresh(t);
    W_OGD_beta_temp = W_final_OGD_beta > thresh(t);
    
    P_d_DMD(t)=sum(W_DMD_temp(:).*W_pos(:))/W_sum;
    P_f_DMD(t)=sum(W_DMD_temp(:).*(1-W_pos(:)))/(W_neg);
    
    P_d_DMD_beta(t)=sum(W_DMD_beta_temp(:).*W_pos(:))/W_sum;
    P_f_DMD_beta(t)=sum(W_DMD_beta_temp(:).*(1-W_pos(:)))/(W_neg);
    
    P_d_OGD(t)=sum(W_OGD_temp(:).*W_pos(:))/W_sum;
    P_f_OGD(t)=sum(W_OGD_temp(:).*(1-W_pos(:)))/(W_neg);
    
    P_d_OGD_beta(t)=sum(W_OGD_beta_temp(:).*W_pos(:))/W_sum;
    P_f_OGD_beta(t)=sum(W_OGD_beta_temp(:).*(1-W_pos(:)))/(W_neg);
    
end

figure(11)
plot(P_f_DMD,P_d_DMD,'b-',P_f_OGD,P_d_OGD,'r--',...
    P_f_DMD_beta,P_d_DMD_beta,'k-',...
    P_f_OGD_beta,P_d_OGD_beta,'c-','linewidth',3)
xlabel('P_f','fontsize',18)
ylabel('P_d','fontsize',18)
curr_leg=legend('Alg 2','OGD','Alg 2, Mismatch','OGD, Mismatch','location','southeast');
set(curr_leg,'fontsize',18);
axis square
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
if savefigs
    print(fig,'ROC','-dpdf')
    savefig('ROC.fig')
end

W_true_top = zeros(size(W_true));
for k = 1:100
    vals = sort(reshape(W_true(:,:,k),[10000,1]));
    W_true_top(:,:,k) = W_true(:,:,k) > vals(9001);
end

no_points = 500;
thresh = linspace(0,max([W_final_DMD(:);W_final_DMD_beta(:); W_final_OGD(:);W_final_OGD_beta(:)]),no_points);

P_d_DMD=zeros(1,no_points);
P_f_DMD=zeros(1,no_points);
P_d_OGD = zeros(1,no_points);
P_f_OGD = zeros(1,no_points);
P_d_DMD_beta=zeros(1,no_points);
P_f_DMD_beta=zeros(1,no_points);
P_d_OGD_beta=zeros(1,no_points);
P_f_OGD_beta=zeros(1,no_points);


W_sum = sum(W_true_top(:));
W_neg = sum(1 - W_true_top(:));
for t=1:length(thresh);
    W_DMD_temp = W_final_DMD > thresh(t);
    W_DMD_beta_temp = W_final_DMD_beta > thresh(t);
    W_OGD_temp = W_final_OGD > thresh(t);
    W_OGD_beta_temp = W_final_OGD_beta > thresh(t);
    
    P_d_DMD(t)=sum(W_DMD_temp(:).*W_true_top(:))/W_sum;
    P_f_DMD(t)=sum(W_DMD_temp(:).*(1-W_true_top(:)))/(W_neg);
    
    P_d_DMD_beta(t)=sum(W_DMD_beta_temp(:).*W_true_top(:))/W_sum;
    P_f_DMD_beta(t)=sum(W_DMD_beta_temp(:).*(1-W_true_top(:)))/(W_neg);
    
    P_d_OGD(t)=sum(W_OGD_temp(:).*W_true_top(:))/W_sum;
    P_f_OGD(t)=sum(W_OGD_temp(:).*(1-W_true_top(:)))/(W_neg);
    
    P_d_OGD_beta(t)=sum(W_OGD_beta_temp(:).*W_true_top(:))/W_sum;
    P_f_OGD_beta(t)=sum(W_OGD_beta_temp(:).*(1-W_true_top(:)))/(W_neg);
    
end

figure(12)
plot(P_f_DMD,P_d_DMD,'b-',P_f_OGD,P_d_OGD,'r--',...
    P_f_DMD_beta,P_d_DMD_beta,'k-',...
    P_f_OGD_beta,P_d_OGD_beta,'c-','linewidth',3)
xlabel('P_f','fontsize',18)
ylabel('P_d','fontsize',18)
curr_leg = legend('Alg 2','OGD','Alg 2, Mismatch','OGD, Mismatch','location','southeast');
set(curr_leg,'fontsize',18)
axis square
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
if savefigs
    print(fig,'ROCTop','-dpdf')
    savefig('ROCTop.fig')
end