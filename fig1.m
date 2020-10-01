function fig1

a = 0.2;
b = 1;
d = 0.5;
q = 1e-3;
alpha = 0.4;
beta = 0.2;
gamma = 0.2;
sigtau = 5;
betamult = 1;
tau = betamult*beta;
sigma = sigtau*tau;

c = linspace(1,20,1001);
ymax = 100;
c(c==0)=[];

figure(1)
clf
set(gcf,'color','w')
set(gcf,'PaperUnits','centimeters')
xSize = 8; ySize = 8;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 100 xSize*50 ySize*50])

labs1 = {'(a)','(b)'};

%%% Approximation
Bcc = max(0,1-(sigma./(tau.*c.^2)));
R0_I_approx = (tau.*c.^2)/sigma;
R0_D_approx = (beta*c.^2)./(d*(1-(1-a)*Bcc)+alpha+gamma);

%%% Full model
% Starting from disease-free and information-free states
R0_I_DF = (tau.*c.^2)/(a*d+sigma);
R0_D_IF = (beta.*c.^2)/(d+alpha+gamma);

% Starting from disease or information present
SU = (d + alpha + gamma).*(alpha.^2 + (-beta.*c.^2 + d + gamma).*alpha + c.^2.*beta.*(b - d))./(beta.^2.*c.^4.*q);
IU = max(0,(alpha.^2 + (-beta.*c.^2 + d + gamma).*alpha + c.^2.*beta.*(b - d)).*(beta.*c.^2 - alpha - d - gamma)./(beta.^2.*c.^4.*q));
list = SU<=0 | IU<=0;
SU(list) = (b-d)/q;
IU(list) = 0;
N = SU+IU;
R0_I = tau.*c.^2.*(((a.*d + gamma + sigma).*IU + SU.*(a.*d + alpha + gamma + sigma)).*N + c.^2.*IU.*beta.*(IU + SU))./(N.*((a.*d + sigma).*(a.*d + alpha + gamma + sigma).*N + c.^2.*IU.*beta.*(a.*d + alpha + sigma)));
R0_I(IU==0) = R0_I_DF(IU==0);

SU = (a.*d + sigma).*((a.^2 - a).*d.^2 + ((-c.^2.*tau + sigma).*a - sigma).*d + b.*c.^2.*tau)./(c.^4.*q.*tau.^2);
SA = max(0,-(-c.^2.*tau + a.*d + sigma).*((a.^2 - a).*d.^2 + ((-c.^2.*tau + sigma).*a - sigma).*d + b.*c.^2.*tau)./(c.^4.*q.*tau.^2));
list = SU<=0 | SA<=0;
SU(list) = (b-d)/q;
SA(list) = 0;
N = SU+SA;
R0_D = beta.*(((sigma + d + alpha + gamma).*SA + SU.*(a.*d + alpha + gamma + sigma)).*N + c.^2.*SA.*tau.*(SU + SA)).*c.^2./(((d + alpha + gamma).*(a.*d + alpha + gamma + sigma).*N + c.^2.*SA.*tau.*(a.*d + alpha + gamma)).*N);
R0_D(SA==0) = R0_D_IF(SA==0);

subplot(2,2,1)
hold on
plot(c,R0_I,'k','linewidth',1.5)
plot(c,R0_I_DF,'b:','linewidth',1.5)
plot(c,R0_I_approx,'r--','linewidth',1.5)
plot(c,R0_I_DF,'b:','linewidth',1.5)
box on
set(gca,'yscale','log')
set(gca,'fontsize',10)
ylim([1,ymax])
text(0,140,labs1{1},'interpreter','latex','fontsize',12)
ylabel('Basic reproductive ratio','interpreter','latex','fontsize',16);
title('Information','interpreter','latex','fontsize',14);

L = legend('$R_0^k$','$\hat{R}_0^k$','$\tilde{R}_0^k$');
set(L,'interpreter','latex')
temp = get(L,'position');
temp(1) = 0.15;
temp(2) = 0.82;
set(L,'position',temp)
x1=xlabel('Contact initiation rate, $c$','interpreter','latex','fontsize',16);
temp=get(x1,'position');
temp(1) = temp(1) + 13;
set(x1,'position',temp);

subplot(2,2,2)
hold on
plot(c,R0_D,'k','linewidth',1.5)
plot(c,R0_D_approx,'r--','linewidth',1.5)
plot(c,R0_D_IF,'b:','linewidth',1.5)
box on
set(gca,'fontsize',10)
set(gca,'yscale','log')
ylim([1,ymax])
text(0,140,labs1{2},'interpreter','latex','fontsize',12)
title('Disease','interpreter','latex','fontsize',14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% if(~exist('fig1.mat','file'))
    clear
    
    % Parameters
    t_max = 1e4;
    a = 0.2;
    b = 1;
    c = 4;
    d = 0.5;
    q = 1e-3;
    alpha = 0.4;
    beta = 0.2;
    gamma = 0.2;
    eqtol = 1e-3;
    sigtau = 5;
    init_pop1 = ((b-d)/(4*q)*ones(1,4));
    init_pop2 = ((b-d)/(2*q)*ones(1,2));    
    RES = 1001;
    
    TAU = beta*logspace(-1,1,RES);
    SIGMA = sigtau*TAU;
    OUT1 = zeros(RES,4);
    OUT2 = zeros(RES,4);
    
    for i=1:RES
        tau=TAU(i);
        sigma=SIGMA(i);
        
        Ac =  max(0,1 - sigma/(c^2*tau));
        
        [~,x1] = sociality_slowinfo_mex(t_max,a,b,c,d,q,alpha,beta,gamma,sigma,tau,eqtol,init_pop1);
        [~,x2] = sociality_fastinfo_mex(t_max,a,b,c,d,q,alpha,beta,gamma,sigma,tau,eqtol,init_pop2);
        
        OUT1(i,:) = x1(end,:)./sum(x1(end,:));
        OUT2(i,:) = [x2(end,1)*[1-Ac,Ac],x2(end,2)*[1-Ac,Ac]]./sum(x2(end,:));
    end
    
    DISPREV1 = sum(OUT1(:,3:4),2);
    DISPREV2 = sum(OUT2(:,3:4),2);
    INFOPREV1 = sum(OUT1(:,[2,4]),2);
    INFOPREV2 = sum(OUT2(:,[2,4]),2);
    
    cols = [166,206,227
        31,120,180
        178,223,138
        51,160,44]/255;
    
    cols2 = [27,158,119
        117,112,179]/255;
    
    clear ans i x1 x2 tau sigma Ac
%     save('fig1.mat')
% else
%     load('fig1.mat')
% end

labs2 = {'(c)','(d)'};
subplot(2,2,3)
hold on
for i=1:4
    plot(TAU/beta,OUT1(:,i),'color',cols(i,:),'linewidth',1.5)
end
for i=1:4
    plot(TAU/beta,OUT2(:,i),'--','color',cols(i,:),'linewidth',1.5)
end
set(gca,'ytick',0:0.2:1)
set(gca,'xscale','log')
set(gca,'fontsize',10)
ylim([0,0.8])
box on
ylabel('Proportion','interpreter','latex','fontsize',16)
text(TAU(1)/beta,1.07*0.8,labs2{1},'interpreter','latex','fontsize',12)
text(TAU(end)/beta,OUT2(end,1)-0.02,'$\frac{S_U}{N}$','interpreter','latex','fontsize',9)
text(TAU(end)/beta,OUT2(end,2)-0.01,'$\frac{S_A}{N}$','interpreter','latex','fontsize',9)
text(TAU(end)/beta,OUT2(end,3)+0.02,'$\frac{I_U}{N}$','interpreter','latex','fontsize',9)
text(TAU(end)/beta,OUT2(end,4),'$\frac{I_A}{N}$','interpreter','latex','fontsize',9)
x1 = xlabel('Relative transmission rate, $\tau/\beta$','interpreter','latex','fontsize',16);
temp=get(x1,'position');
temp(1) = temp(1) + 17;
set(x1,'position',temp);

subplot(2,2,4)
hold on
plot(TAU/beta,DISPREV1,'color',cols2(1,:),'linewidth',1.5)
plot(TAU/beta,INFOPREV1,'color',cols2(2,:),'linewidth',1.5)
plot(TAU/beta,DISPREV2,'--','color',cols2(1,:),'linewidth',1.5)
plot(TAU/beta,INFOPREV2,'--','color',cols2(2,:),'linewidth',1.5)
set(gca,'ytick',0:0.2:1)
set(gca,'xscale','log')
set(gca,'fontsize',10)
ylim([0,1])
box on
text(TAU(1)/beta,1.07,labs2{2},'interpreter','latex','fontsize',12)
text(TAU(end)/beta,DISPREV2(end)+0.04,'$\frac{I_U+I_A}{N}$','interpreter','latex','fontsize',10)
text(TAU(end)/beta,INFOPREV2(end)-0.04,'$\frac{S_A+I_A}{N}$','interpreter','latex','fontsize',10)
legend('infected','aware','location','southeast')

% save2pdf('fig1.pdf');