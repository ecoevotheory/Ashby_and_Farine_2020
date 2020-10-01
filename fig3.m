function fig3

% Note: arrows require arrow.m from Mathworks File Exchange

if(~exist('fig3.mat','file'))
    fig3_data;
end
load('fig3.mat')
beta=0.2;
gamma=0.2;
tau=0.2;
i1 = 1;

simstart = [1.25,3.8,4.2,4.5];
simsigtau = [1,3,3,15];

% Plot outcome if necessary
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


figure(3)
clf
set(gcf,'color','w')
set(gcf,'PaperUnits','centimeters')
xSize = 14; ySize = 14;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 100 xSize*50 ySize*50])
colormap([1,1,1;0.5,0.5,0.5;0.75,0.75,0.75])

labs = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'};
sub=0;
for i3=1:length(A)
    for i2=1:length(ALPHA)
        singstrat = permute(SINGSTRAT_SIGTAU_FAST(:,i2,i3,:),[1,4,3,2]);
        outcome = permute(OUTCOME_SIGTAU_FAST(:,i2,i3,:),[1,4,3,2]);
        sub=sub+1;
        subplot(length(A),length(ALPHA),sub)
        
        hold on
        plot(SIGTAU,sqrt(SIGTAU),':','color',0.8*[1,1,1],'linewidth',1.5) % Fast information threshold (no disease)
        plot(SIGTAU,sqrt(A(i3)*d/tau + SIGTAU),'color',0.8*[1,1,1],'linewidth',1.5) % Slow information threshold (no disease)
        plot(SIGTAU,sqrt((d+ALPHA(i2)+gamma)/beta) + 0*SIGTAU,'color',0.8*[1,1,1],'linewidth',1.5) % Disease threshold (no information)
        
        % rebuild arrays
        outcome2 = NaN*zeros(xres,9);
        singstrat2 = NaN*zeros(xres,9);
        for j=1:9
            list = find(sum(outcome==j,2));
            outcome2(list,j)=j;
            for k=1:length(list)
                singstrat2(list(k),j) = singstrat(list(k),outcome(list(k),:)==j);
            end
        end
        
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
                plot(SIGTAU(list_seg),singstrat2(list_seg,j),':','color',cols(j,:),'linewidth',1.5)
                list(1:seg_length)=[];
            end
        end
        xlim([floor(SIGTAU(1)),ceil(SIGTAU(end))])
        
        singstrat = permute(SINGSTRAT_SIGTAU_SLOW(:,i1,i2,i3,:),[1,5,3,4,2]);
        outcome = permute(OUTCOME_SIGTAU_SLOW(:,i1,i2,i3,:),[1,5,3,4,2]);
        
        % rebuild arrays
        outcome2 = NaN*zeros(xres,9);
        singstrat2 = NaN*zeros(xres,9);
        for j=1:9
            list = find(sum(outcome==j,2));
            outcome2(list,j)=j;
            for k=1:length(list)
                singstrat2(list(k),j) = singstrat(list(k),outcome(list(k),:)==j);
            end
        end
        
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
                plot(SIGTAU(list_seg),singstrat2(list_seg,j),'-','color',cols(j,:),'linewidth',1.5)
                list(1:seg_length)=[];
            end
        end
        xlim([SIGTAU(1),SIGTAU(end)])
        set(gca,'xtick',[0.1,1,10])
        set(gca,'xscale','log')
        ymax = 10;
        ylim([0,ymax])
        box on
        set(gca,'fontsize',10)
        text(SIGTAU(1),ymax*1.05,labs{sub},'fontsize',12)
        if(sub==4)
            ylabel('Contact initiation rate, $c$','interpreter','latex','fontsize',16);
        elseif(sub==8)
            xlabel('Information expiry rate/transmission rate, $\sigma/\tau$','interpreter','latex','fontsize',16);
        end
        title(strcat('$a=',num2str(A(i3)),', \alpha=',num2str(ALPHA(i2)),'$'),'interpreter','latex','fontsize',10);
        
        if(sub==4)
            plot(2,1.6,'k*')
        elseif(sub==5)
            plot(simsigtau,simstart,'kx','markersize',5)
        end
        if(exist('arrow.m','file'))
            
            if(sub==1)
                arrow([0.15,2.7],[0.15,1.8],'length',6)
                arrow([0.15,0.7],[0.15,1.6],'length',6)
                arrow([0.15,7],[0.15,7.9],'length',6)
                
                arrow([1.6,1.5],[1.6,3],'length',6)
                
                arrow([18,4],[18,2.5],'length',6)
                arrow([18,5],[18,6.5],'length',6)
            elseif(sub==2)
                arrow([0.5,3.2],[0.5,2.2],'length',6)
                arrow([0.5,0.9],[0.5,1.9],'length',6)
                
                arrow([18,4.5],[18,3],'length',6)
                arrow([18,5.5],[18,7],'length',6)
            elseif(sub==3)
                arrow([0.5,3.4],[0.5,2.4],'length',6)
                arrow([0.5,0.9],[0.5,1.9],'length',6)
                arrow([18,4.9],[18,3.4],'length',6)
                arrow([18,5.9],[18,7.4],'length',6)
            elseif(sub==4)
                arrow([0.5,3.2],[0.5,2.3],'length',6)
                arrow([0.5,1.2],[0.5,2.1],'length',6)
                arrow([0.5,7],[0.5,7.9],'length',6)
                
                plot(2,1.6,'k*')
                
                arrow([2,1.75],[2,3],'length',6)
                
                arrow([18,4],[18,2.5],'length',6)
                arrow([18,5],[18,6.5],'length',6)
            elseif(sub==5)
                arrow([0.5,3.2],[0.5,2.4],'length',6)
                arrow([0.5,1.1],[0.5,1.9],'length',6)
                
                arrow([18,4.5],[18,3],'length',6)
                arrow([18,5.5],[18,7],'length',6)
                
                plot(simsigtau,simstart,'kx','markersize',5)
            elseif(sub==6)
                arrow([0.5,3.6],[0.5,2.6],'length',6)
                arrow([0.5,1.1],[0.5,2.1],'length',6)
                
                arrow([18,5.2],[18,3.7],'length',6)
                arrow([18,6.2],[18,7.7],'length',6)
            elseif(sub==7)
                arrow([0.5,4],[0.5,2.5],'length',6)
                arrow([18,4.5],[18,3],'length',6)
                arrow([18,5.5],[18,7],'length',6)
            elseif(sub==8)
                arrow([0.5,3.4],[0.5,2.5],'length',6)
                arrow([0.5,1.3],[0.5,2.2],'length',6)
                
                arrow([18,5],[18,3.5],'length',6)
                arrow([18,6],[18,7.5],'length',6)
            elseif(sub==9)
                arrow([0.5,3.5],[0.5,2.6],'length',6)
                arrow([0.5,1.4],[0.5,2.3],'length',6)
                
                arrow([18,5.9],[18,4.4],'length',6)
                arrow([18,6.9],[18,8.4],'length',6)
            end
        end
    end
end

% save2pdf('fig3.pdf')

function fig3_data

% Fixed and default parameter values
t_max = 1000;
b = 1;
cmin = 0;
cmax = 20;
d = 0.5;
q = 1e-3;
beta = 0.2;
gamma = 0.2;
xres = 101;
res1 = 1001;

% Variables
SIGTAU = logspace(-1,log10(20),xres);
TAUBETA = logspace(-1,2,xres);
A = [0,0.2,0.5];
ALPHA = [0.3,0.4,0.5];
SIGTAU_DEFAULT = 1;%[0.1,1,10];
TAUBETA_DEFAULT = 1;% [0.1,1,10];

% OUTPUTS
SINGSTRAT_SIGTAU_FAST1 = NaN*zeros(xres,length(ALPHA),length(A));
SINGSTRAT_SIGTAU_FAST2 = NaN*zeros(xres,length(ALPHA),length(A));
SINGSTRAT_SIGTAU_FAST3 = NaN*zeros(xres,length(ALPHA),length(A));
SINGSTRAT_SIGTAU_FAST4 = NaN*zeros(xres,length(ALPHA),length(A));
SINGSTRAT_SIGTAU_FAST5 = NaN*zeros(xres,length(ALPHA),length(A));
OUTCOME_SIGTAU_FAST1 = NaN*zeros(xres,length(ALPHA),length(A));
OUTCOME_SIGTAU_FAST2 = NaN*zeros(xres,length(ALPHA),length(A));
OUTCOME_SIGTAU_FAST3 = NaN*zeros(xres,length(ALPHA),length(A));
OUTCOME_SIGTAU_FAST4 = NaN*zeros(xres,length(ALPHA),length(A));
OUTCOME_SIGTAU_FAST5 = NaN*zeros(xres,length(ALPHA),length(A));
NUM_OUTCOMES_SIGTAU_FAST = NaN*zeros(xres,length(ALPHA),length(A));

SINGSTRAT_SIGTAU_SLOW1 = NaN*zeros(xres,length(TAUBETA_DEFAULT),length(ALPHA),length(A));
SINGSTRAT_SIGTAU_SLOW2 = NaN*zeros(xres,length(TAUBETA_DEFAULT),length(ALPHA),length(A));
SINGSTRAT_SIGTAU_SLOW3 = NaN*zeros(xres,length(TAUBETA_DEFAULT),length(ALPHA),length(A));
SINGSTRAT_SIGTAU_SLOW4 = NaN*zeros(xres,length(TAUBETA_DEFAULT),length(ALPHA),length(A));
SINGSTRAT_SIGTAU_SLOW5 = NaN*zeros(xres,length(TAUBETA_DEFAULT),length(ALPHA),length(A));
OUTCOME_SIGTAU_SLOW1 = NaN*zeros(xres,length(TAUBETA_DEFAULT),length(ALPHA),length(A));
OUTCOME_SIGTAU_SLOW2 = NaN*zeros(xres,length(TAUBETA_DEFAULT),length(ALPHA),length(A));
OUTCOME_SIGTAU_SLOW3 = NaN*zeros(xres,length(TAUBETA_DEFAULT),length(ALPHA),length(A));
OUTCOME_SIGTAU_SLOW4 = NaN*zeros(xres,length(TAUBETA_DEFAULT),length(ALPHA),length(A));
OUTCOME_SIGTAU_SLOW5 = NaN*zeros(xres,length(TAUBETA_DEFAULT),length(ALPHA),length(A));
NUM_OUTCOMES_SIGTAU_SLOW = NaN*zeros(xres,length(TAUBETA_DEFAULT),length(ALPHA),length(A));
SINGSTRAT_TAUBETA_SLOW1 = NaN*zeros(xres,length(SIGTAU_DEFAULT),length(ALPHA),length(A));
SINGSTRAT_TAUBETA_SLOW2 = NaN*zeros(xres,length(SIGTAU_DEFAULT),length(ALPHA),length(A));
SINGSTRAT_TAUBETA_SLOW3 = NaN*zeros(xres,length(SIGTAU_DEFAULT),length(ALPHA),length(A));
SINGSTRAT_TAUBETA_SLOW4 = NaN*zeros(xres,length(SIGTAU_DEFAULT),length(ALPHA),length(A));
SINGSTRAT_TAUBETA_SLOW5 = NaN*zeros(xres,length(SIGTAU_DEFAULT),length(ALPHA),length(A));
OUTCOME_TAUBETA_SLOW1 = NaN*zeros(xres,length(SIGTAU_DEFAULT),length(ALPHA),length(A));
OUTCOME_TAUBETA_SLOW2 = NaN*zeros(xres,length(SIGTAU_DEFAULT),length(ALPHA),length(A));
OUTCOME_TAUBETA_SLOW3 = NaN*zeros(xres,length(SIGTAU_DEFAULT),length(ALPHA),length(A));
OUTCOME_TAUBETA_SLOW4 = NaN*zeros(xres,length(SIGTAU_DEFAULT),length(ALPHA),length(A));
OUTCOME_TAUBETA_SLOW5 = NaN*zeros(xres,length(SIGTAU_DEFAULT),length(ALPHA),length(A));
NUM_OUTCOMES_TAUBETA_SLOW = NaN*zeros(xres,length(SIGTAU_DEFAULT),length(ALPHA),length(A));

% Fast info version depends on SIGTAU but not on TAUBETA, so can fix tau
COUNT = 0;
TOTAL = length(A)*length(ALPHA)*(1 + length(TAUBETA_DEFAULT) + length(SIGTAU_DEFAULT))
tau = beta;
for i3=1:length(A)
    for i2=1:length(ALPHA)
        % Sweep over SIGTAU        
        tic;
        parfor i1=1:length(SIGTAU)
            [singstrat,outcome] = singstrat_fastinfo(A(i3),cmin,cmax,d,ALPHA(i2),beta,gamma,SIGTAU(i1)*tau,tau,res1);
            
            NUM_OUTCOMES_SIGTAU_FAST(i1,i2,i3) = length(singstrat);
            SINGSTRAT_SIGTAU_FAST1(i1,i2,i3) = singstrat(1);
            OUTCOME_SIGTAU_FAST1(i1,i2,i3) = outcome(1);
            if(NUM_OUTCOMES_SIGTAU_FAST(i1,i2,i3)>1)
                SINGSTRAT_SIGTAU_FAST2(i1,i2,i3) = singstrat(2);
                OUTCOME_SIGTAU_FAST2(i1,i2,i3) = outcome(2);
                if(NUM_OUTCOMES_SIGTAU_FAST(i1,i2,i3)>2)
                    SINGSTRAT_SIGTAU_FAST3(i1,i2,i3) = singstrat(3);
                    OUTCOME_SIGTAU_FAST3(i1,i2,i3) = outcome(3);
                    if(NUM_OUTCOMES_SIGTAU_FAST(i1,i2,i3)>3)
                        SINGSTRAT_SIGTAU_FAST4(i1,i2,i3) = singstrat(4);
                        OUTCOME_SIGTAU_FAST4(i1,i2,i3) = outcome(4);
                        if(NUM_OUTCOMES_SIGTAU_FAST(i1,i2,i3)>4)
                            SINGSTRAT_SIGTAU_FAST5(i1,i2,i3) = singstrat(5);
                            OUTCOME_SIGTAU_FAST5(i1,i2,i3) = outcome(5);
                        end
                    end
                end
            end
        end
        toc;
        COUNT=COUNT+1;
        PROGRESS = COUNT/TOTAL
    end
end

for i4=1:length(A)
    for i3=1:length(ALPHA)
         
        % Fix TAUBETA and sweep over SIGTAU
        for i2=1:length(TAUBETA_DEFAULT)
            tic;
            tau = TAUBETA_DEFAULT(i2)*beta;
            parfor i1=1:length(SIGTAU)
                [singstrat,outcome] = singstrat_slowinfo(t_max,A(i4),b,cmin,cmax,d,q,ALPHA(i3),beta,gamma,SIGTAU(i1)*tau,tau,res1);
                
                NUM_OUTCOMES_SIGTAU_SLOW(i1,i2,i3,i4) = length(singstrat);
                SINGSTRAT_SIGTAU_SLOW1(i1,i2,i3,i4) = singstrat(1);
                OUTCOME_SIGTAU_SLOW1(i1,i2,i3,i4) = outcome(1);
                if(NUM_OUTCOMES_SIGTAU_SLOW(i1,i2,i3,i4)>1)
                    SINGSTRAT_SIGTAU_SLOW2(i1,i2,i3,i4) = singstrat(2);
                    OUTCOME_SIGTAU_SLOW2(i1,i2,i3,i4) = outcome(2);
                    if(NUM_OUTCOMES_SIGTAU_SLOW(i1,i2,i3,i4)>2)
                        SINGSTRAT_SIGTAU_SLOW3(i1,i2,i3,i4) = singstrat(3);
                        OUTCOME_SIGTAU_SLOW3(i1,i2,i3,i4) = outcome(3);
                        if(NUM_OUTCOMES_SIGTAU_SLOW(i1,i2,i3,i4)>3)
                            SINGSTRAT_SIGTAU_SLOW4(i1,i2,i3,i4) = singstrat(4);
                            OUTCOME_SIGTAU_SLOW4(i1,i2,i3,i4) = outcome(4);
                            if(NUM_OUTCOMES_SIGTAU_SLOW(i1,i2,i3,i4)>4)
                                SINGSTRAT_SIGTAU_SLOW5(i1,i2,i3,i4) = singstrat(5);
                                OUTCOME_SIGTAU_SLOW5(i1,i2,i3,i4) = outcome(5);
                            end
                        end
                    end
                end
            end
            toc;
            COUNT=COUNT+1;
            PROGRESS = COUNT/TOTAL
        end
        
        % Fix SIGTAU and sweep over TAUBETA
        for i2=1:length(SIGTAU_DEFAULT)
            tic;
            parfor i1=1:length(TAUBETA)
                tau = TAUBETA(i1)*beta;
                [singstrat,outcome] = singstrat_slowinfo(t_max,A(i4),b,cmin,cmax,d,q,ALPHA(i3),beta,gamma,SIGTAU_DEFAULT(i2)*tau,tau,res1);
                
                NUM_OUTCOMES_TAUBETA_SLOW(i1,i2,i3,i4) = length(singstrat);
                SINGSTRAT_TAUBETA_SLOW1(i1,i2,i3,i4) = singstrat(1);
                OUTCOME_TAUBETA_SLOW1(i1,i2,i3,i4) = outcome(1);
                if(NUM_OUTCOMES_TAUBETA_SLOW(i1,i2,i3,i4)>1)
                    SINGSTRAT_TAUBETA_SLOW2(i1,i2,i3,i4) = singstrat(2);
                    OUTCOME_TAUBETA_SLOW2(i1,i2,i3,i4) = outcome(2);
                    if(NUM_OUTCOMES_TAUBETA_SLOW(i1,i2,i3,i4)>2)
                        SINGSTRAT_TAUBETA_SLOW3(i1,i2,i3,i4) = singstrat(3);
                        OUTCOME_TAUBETA_SLOW3(i1,i2,i3,i4) = outcome(3);
                        if(NUM_OUTCOMES_TAUBETA_SLOW(i1,i2,i3,i4)>3)
                            SINGSTRAT_TAUBETA_SLOW4(i1,i2,i3,i4) = singstrat(4);
                            OUTCOME_TAUBETA_SLOW4(i1,i2,i3,i4) = outcome(4);
                            if(NUM_OUTCOMES_TAUBETA_SLOW(i1,i2,i3,i4)>4)
                                SINGSTRAT_TAUBETA_SLOW5(i1,i2,i3,i4) = singstrat(5);
                                OUTCOME_TAUBETA_SLOW5(i1,i2,i3,i4) = outcome(5);
                            end
                        end
                    end
                end
            end            
            toc;
            COUNT=COUNT+1;
            PROGRESS = COUNT/TOTAL
        end
    end
end

% Combine results
SINGSTRAT_SIGTAU_FAST = cat(4,SINGSTRAT_SIGTAU_FAST1,SINGSTRAT_SIGTAU_FAST2,SINGSTRAT_SIGTAU_FAST3,SINGSTRAT_SIGTAU_FAST4,SINGSTRAT_SIGTAU_FAST5);
OUTCOME_SIGTAU_FAST = cat(4,OUTCOME_SIGTAU_FAST1,OUTCOME_SIGTAU_FAST2,OUTCOME_SIGTAU_FAST3,OUTCOME_SIGTAU_FAST4,OUTCOME_SIGTAU_FAST5);
SINGSTRAT_SIGTAU_SLOW = cat(5,SINGSTRAT_SIGTAU_SLOW1,SINGSTRAT_SIGTAU_SLOW2,SINGSTRAT_SIGTAU_SLOW3,SINGSTRAT_SIGTAU_SLOW4,SINGSTRAT_SIGTAU_SLOW5);
OUTCOME_SIGTAU_SLOW = cat(5,OUTCOME_SIGTAU_SLOW1,OUTCOME_SIGTAU_SLOW2,OUTCOME_SIGTAU_SLOW3,OUTCOME_SIGTAU_SLOW4,OUTCOME_SIGTAU_SLOW5);
SINGSTRAT_TAUBETA_SLOW = cat(5,SINGSTRAT_TAUBETA_SLOW1,SINGSTRAT_TAUBETA_SLOW2,SINGSTRAT_TAUBETA_SLOW3,SINGSTRAT_TAUBETA_SLOW4,SINGSTRAT_TAUBETA_SLOW5);
OUTCOME_TAUBETA_SLOW = cat(5,OUTCOME_TAUBETA_SLOW1,OUTCOME_TAUBETA_SLOW2,OUTCOME_TAUBETA_SLOW3,OUTCOME_TAUBETA_SLOW4,OUTCOME_TAUBETA_SLOW5);

clear i1 i2 i3 i4 singstrat outcome tau ans COUNT TOTAL PROGRESS
clear SINGSTRAT_SIGTAU_FAST1 SINGSTRAT_SIGTAU_FAST2 SINGSTRAT_SIGTAU_FAST3 SINGSTRAT_SIGTAU_FAST4 SINGSTRAT_SIGTAU_FAST5 OUTCOME_SIGTAU_FAST1 OUTCOME_SIGTAU_FAST2 OUTCOME_SIGTAU_FAST3 OUTCOME_SIGTAU_FAST4 OUTCOME_SIGTAU_FAST5 SINGSTRAT_SIGTAU_SLOW1 SINGSTRAT_SIGTAU_SLOW2 SINGSTRAT_SIGTAU_SLOW3 SINGSTRAT_SIGTAU_SLOW4 SINGSTRAT_SIGTAU_SLOW5 OUTCOME_SIGTAU_SLOW1 OUTCOME_SIGTAU_SLOW2 OUTCOME_SIGTAU_SLOW3 OUTCOME_SIGTAU_SLOW4 OUTCOME_SIGTAU_SLOW5 SINGSTRAT_TAUBETA_SLOW1 SINGSTRAT_TAUBETA_SLOW2 SINGSTRAT_TAUBETA_SLOW3 SINGSTRAT_TAUBETA_SLOW4 SINGSTRAT_TAUBETA_SLOW5 OUTCOME_TAUBETA_SLOW1 OUTCOME_TAUBETA_SLOW2 OUTCOME_TAUBETA_SLOW3 OUTCOME_TAUBETA_SLOW4 OUTCOME_TAUBETA_SLOW5
save('fig3.mat')