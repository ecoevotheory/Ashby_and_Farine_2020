function fig7

if(~exist('fig7.mat','file'))
    fig7_data;
end
load('fig7.mat')
gamma = 0.2;


PARAM = SIGTAU;
cols = [251,154,153
    171,217,233
    200,10,10
    69,117,180
    0,0,0
    0,0,0
    0,0,0
    0,0,0
    0,0,0]/255;
% 1) pseudo-unstable (pink)
% 2) psuedo-stable (light blue)
% 3) repeller (red)
% 4) CSS (dark blue)

figure(7)
clf
set(gcf,'color','w')
set(gcf,'PaperUnits','centimeters')
xSize = 18; ySize = 4.5;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 100 xSize*50 ySize*50])

labs = {'(a)','(b)','(c)'};
ymax = 5;
for i2=1:length(A)
    subplot(1,length(A),i2)
    singstrat = permute(SINGSTRAT_SIGTAU_FAST(:,i2,:),[1,3,2]);
    outcome = permute(OUTCOME_SIGTAU_FAST(:,i2,:),[1,3,2]);
    
    outcome2 = NaN*zeros(xres,9);
    singstrat2 = NaN*zeros(xres,9);
    for j=1:9
        list = find(sum(outcome==j,2));
        outcome2(list,j)=j;
        for k=1:length(list)
            singstrat2(list(k),j) = singstrat(list(k),find(outcome(list(k),:)==j,1));
        end
    end
    ALPHA = (singstrat2.^2.*(a.*d + gamma).*tau - d.*repmat(SIGTAU',[1,9]).*(-1 + a))./(singstrat2.^2.*tau);
    ALPHA(:,[1,3,5:9])=0;
    ALPHA(~isnan(ALPHA(:,2)),2) = gamma+d;
    ALPHA(isnan(ALPHA))=0;
    ALPHA=sum(ALPHA,2);
    ALPHA(ALPHA==0)=NaN;
    BETA = kappa*sqrt(ALPHA);
    
    yyaxis left
    hold on
    plot(SIGTAU,sqrt(SIGTAU),'-','color',0.8*[1,1,1],'linewidth',1.5) % Fast information threshold (no disease)
    plot(SIGTAU,sqrt((d+ALPHA+gamma)./BETA),'-','color',0.8*[1,1,1],'linewidth',1.5) % Disease threshold (no information)
    
    % Plot each outcome
    for j=1:9
        list = find(outcome2(:,j)==j);
        while(~isempty(list))
            list_diff = diff(list);
            seg_length = find(list_diff>1);
            if(isempty(seg_length))
                seg_length=length(list);
            end
            list_seg = list(1:seg_length);
            yyaxis left
            hold on
            plot(PARAM(list_seg),singstrat2(list_seg,j),'-','color',cols(j,:),'linewidth',1.5)
            list(1:seg_length)=[];
        end
    end
    xlim([floor(PARAM(1)),ceil(PARAM(end))])
    set(gca,'xtick',[0.1,1,10])
    set(gca,'xscale','log')
    yyaxis left
    set(gca,'fontsize',10)
    set(gca,'ycolor','k')
    ylim([0,ymax])
    text(SIGTAU(1),ymax*1.05,labs{i2},'fontsize',12)
    yyaxis right
    plot(SIGTAU,ALPHA,'k-','linewidth',1.5)
    set(gca,'ycolor','k')
    set(gca,'fontsize',10)
    ylim([0,1])
    box on
    title(strcat('$a=',num2str(A(i2)),'$'),'interpreter','latex','fontsize',12);    
end

for i=1:3
subplot(1,3,i)
yyaxis left
ylabel('Contact initiation rate, $c$','interpreter','latex','fontsize',16);
yyaxis right
ylabel('Virulence, $\alpha$','interpreter','latex','fontsize',16);
end

subplot(1,3,2)
xlabel({'Information expiry rate/transmission rate, $\sigma/\tau$'},'interpreter','latex','fontsize',16);

% save2pdf('fig7.pdf')

function fig7_data

% Fixed and default parameter values
cmin = 0;
cmax = 5;
a = 0.2;
d = 0.5;
kappa = 1;
gamma = 0.2;
xres = 10001;
res1 = 1001;
tau = 1;

% Variables
SIGTAU = logspace(-1,log10(20),xres);
A = [0,0.2,0.5];

% OUTPUTS
SINGSTRAT_SIGTAU_FAST1 = NaN*zeros(xres,length(A));
SINGSTRAT_SIGTAU_FAST2 = NaN*zeros(xres,length(A));
SINGSTRAT_SIGTAU_FAST3 = NaN*zeros(xres,length(A));
SINGSTRAT_SIGTAU_FAST4 = NaN*zeros(xres,length(A));
SINGSTRAT_SIGTAU_FAST5 = NaN*zeros(xres,length(A));
OUTCOME_SIGTAU_FAST1 = NaN*zeros(xres,length(A));
OUTCOME_SIGTAU_FAST2 = NaN*zeros(xres,length(A));
OUTCOME_SIGTAU_FAST3 = NaN*zeros(xres,length(A));
OUTCOME_SIGTAU_FAST4 = NaN*zeros(xres,length(A));
OUTCOME_SIGTAU_FAST5 = NaN*zeros(xres,length(A));
NUM_OUTCOMES_SIGTAU_FAST = NaN*zeros(xres,length(A));

% Fast info version depends on SIGTAU but not on TAUBETA, so can fix tau
COUNT = 0;
TOTAL = length(A)
for i2=1:length(A)
    % Sweep over SIGTAU
    tic;
    parfor i1=1:length(SIGTAU)
        [singstrat,outcome] = singstrat_fastinfo_coevo_approx(A(i2),cmin,cmax,d,kappa,gamma,SIGTAU(i1)*tau,tau,res1);
        
        NUM_OUTCOMES_SIGTAU_FAST(i1,i2) = length(singstrat);
        SINGSTRAT_SIGTAU_FAST1(i1,i2) = singstrat(1);
        OUTCOME_SIGTAU_FAST1(i1,i2) = outcome(1);
        if(NUM_OUTCOMES_SIGTAU_FAST(i1,i2)>1)
            SINGSTRAT_SIGTAU_FAST2(i1,i2) = singstrat(2);
            OUTCOME_SIGTAU_FAST2(i1,i2) = outcome(2);
            if(NUM_OUTCOMES_SIGTAU_FAST(i1,i2)>2)
                SINGSTRAT_SIGTAU_FAST3(i1,i2) = singstrat(3);
                OUTCOME_SIGTAU_FAST3(i1,i2) = outcome(3);
                if(NUM_OUTCOMES_SIGTAU_FAST(i1,i2)>3)
                    SINGSTRAT_SIGTAU_FAST4(i1,i2) = singstrat(4);
                    OUTCOME_SIGTAU_FAST4(i1,i2) = outcome(4);
                    if(NUM_OUTCOMES_SIGTAU_FAST(i1,i2)>4)
                        SINGSTRAT_SIGTAU_FAST5(i1,i2) = singstrat(5);
                        OUTCOME_SIGTAU_FAST5(i1,i2) = outcome(5);
                    end
                end
            end
        end
    end
    toc;
    COUNT=COUNT+1;
    PROGRESS = COUNT/TOTAL
end

% Combine results
SINGSTRAT_SIGTAU_FAST = cat(3,SINGSTRAT_SIGTAU_FAST1,SINGSTRAT_SIGTAU_FAST2,SINGSTRAT_SIGTAU_FAST3,SINGSTRAT_SIGTAU_FAST4,SINGSTRAT_SIGTAU_FAST5);
OUTCOME_SIGTAU_FAST = cat(3,OUTCOME_SIGTAU_FAST1,OUTCOME_SIGTAU_FAST2,OUTCOME_SIGTAU_FAST3,OUTCOME_SIGTAU_FAST4,OUTCOME_SIGTAU_FAST5);

clear i1 i2 i3 i4 singstrat outcome ans COUNT TOTAL PROGRESS
clear SINGSTRAT_SIGTAU_FAST1 SINGSTRAT_SIGTAU_FAST2 SINGSTRAT_SIGTAU_FAST3 SINGSTRAT_SIGTAU_FAST4 SINGSTRAT_SIGTAU_FAST5 OUTCOME_SIGTAU_FAST1 OUTCOME_SIGTAU_FAST2 OUTCOME_SIGTAU_FAST3 OUTCOME_SIGTAU_FAST4 OUTCOME_SIGTAU_FAST5 SINGSTRAT_SIGTAU_SLOW1 SINGSTRAT_SIGTAU_SLOW2 SINGSTRAT_SIGTAU_SLOW3 SINGSTRAT_SIGTAU_SLOW4 SINGSTRAT_SIGTAU_SLOW5 OUTCOME_SIGTAU_SLOW1 OUTCOME_SIGTAU_SLOW2 OUTCOME_SIGTAU_SLOW3 OUTCOME_SIGTAU_SLOW4 OUTCOME_SIGTAU_SLOW5 SINGSTRAT_TAUBETA_SLOW1 SINGSTRAT_TAUBETA_SLOW2 SINGSTRAT_TAUBETA_SLOW3 SINGSTRAT_TAUBETA_SLOW4 SINGSTRAT_TAUBETA_SLOW5 OUTCOME_TAUBETA_SLOW1 OUTCOME_TAUBETA_SLOW2 OUTCOME_TAUBETA_SLOW3 OUTCOME_TAUBETA_SLOW4 OUTCOME_TAUBETA_SLOW5
save('fig7.mat')
