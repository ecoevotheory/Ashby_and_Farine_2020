function fig6

a = 0.2;
d = 0.5;
gamma = 0.2;
sigtau = 6;
tau = 0.2;
TAU = [tau,tau/10];
sigma = sigtau.*tau;
SIGMA = sigtau*TAU;
kappa = 1;

c = linspace(0,10,1001);
c(c==0)=[];

figure(6)
clf
set(gcf,'color','w')
set(gcf,'PaperUnits','centimeters')
xSize = 5; ySize = 4.5;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 100 xSize*50 ySize*50])

alpha = (a.*c.^2.*d.*tau + c.^2.*gamma.*tau - a.*d.*sigma + d.*sigma)./(c.^2.*tau);
beta = kappa*sqrt(alpha);
disprev = (c.^4.*tau.*beta - tau.*(a.*d + alpha + gamma).*c.^2 + d.*sigma.*(-1 + a))./(c.^4.*tau.*beta);
Bcc = max(0,1-(sigma./(tau.*c.^2)));

alpha(Bcc==0)=gamma+d;
beta(Bcc==0)=kappa*(sqrt(alpha(Bcc==0)));
disprev(Bcc==0)=(c(Bcc==0).^2.*kappa.*sqrt(alpha(Bcc==0)) - gamma - alpha(Bcc==0) - d)./(kappa.*c(Bcc==0).^2.*sqrt(alpha(Bcc==0)));

R0_D_approx = (beta.*c.^2)./(d*(1-(1-a)*Bcc)+alpha+gamma);

alpha(R0_D_approx<1)=NaN;
disprev(R0_D_approx<1)=NaN;

yyaxis left
plot(c,alpha,'k','linewidth',1.5)
hold on
plot(c(find(isnan(alpha),1,'last'))*[1,1],[0,ceil(10*max(alpha))/10],'k:','linewidth',1.5)
plot(c(find(Bcc==0,1,'last'))*[1,1],[0,ceil(10*max(alpha))/10],'k:','linewidth',1.5)
set(gca,'fontsize',10)
ylabel('Optimal virulence, $\alpha^*$','interpreter','latex','fontsize',16);
text(c(end)*0.75,alpha(end)*0.9,'$\alpha^*$','interpreter','latex','fontsize',12);
text(0.7,max(alpha)*0.75,'$\tilde{R}_0^D<1$','interpreter','latex','fontsize',12,'rotation',90);
text(1.9,max(alpha)*0.75,'$\tilde{R}_0^I<1$','interpreter','latex','fontsize',12,'rotation',90);
set(gca,'ycolor','k')

yyaxis right
plot(c,disprev,'k--','linewidth',1.5)
text(c(end)*0.75,disprev(end)*0.9,'$\frac{I_U+I_A}{N}$','interpreter','latex','fontsize',12);
set(gca,'ycolor','k')

xlabel('Contact initiation rate, $c$','interpreter','latex','fontsize',16);
ylabel('Disease prevalence, $\frac{I_U+I_A}{N}$','interpreter','latex','fontsize',16);
box on

C = c;
alpha = linspace(0,2,1001);
alpha(alpha==0)=[];
for k=1:2
    tau = TAU(k);
    sigma = SIGMA(k);
    for i=1:length(C)
        c = C(i);
        if(tau.*c.^2/(a*d+sigma)<1)
            alpha_ss_slow(i,1) = gamma+d;
        else
            R0_D = -((a.^2 - a).*d.^2 - (c.^2.*tau - sigma).*(a - 1).*d + c.^2.*tau.*(c.^2.*tau + alpha + gamma)).*kappa.*sqrt(alpha)./(tau.*((a.^2 - a).*d.^2 + (-a.*c.^2.*tau + a.*sigma - alpha - gamma - sigma).*d - (alpha + gamma).*(c.^2.*tau + alpha + gamma)));
            alpha_ss_slow(i,1) = alpha(find(R0_D==max(R0_D),1));
        end
    end
    yyaxis left
    if(k==1)
        plot(C,alpha_ss_slow,'r--','linewidth',1.5)
    else
        plot(C,alpha_ss_slow,'r-.','linewidth',1.5)
    end
end

% save2pdf('fig6.pdf');