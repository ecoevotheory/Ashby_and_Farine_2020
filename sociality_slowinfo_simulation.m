function [SOCIALITY,DISPREV,INFOPREV] = sociality_slowinfo_simulation(t_max,a,b,cmin,cmax,cstart,d,q,alpha,beta,gamma,sigma,tau,res0,nevol,plotflag)

% Note: ensure sociality_slowinfo_competition_mex.c is compiled first

eqtol = 1e-3;
exttol = 1e-5;

C = linspace(cmin,cmax,res0);
SOCIALITY = zeros(nevol,res0);
DISPREV = zeros(nevol,1);
INFOPREV = zeros(nevol,1);

% initial conditions
initial = find(C>=cstart,1);
c_current = C(initial);
index_current = initial;
strain_total = 1;
init_pop = ((b-d)/(4*q)*ones(1,4));
for ievol=1:nevol
    
    % Run solver
    [~,SA,SB,IA,IB,~] = sociality_slowinfo_competition_mex(t_max,a,b,c_current,d,q,alpha,beta,gamma,sigma,tau,eqtol,init_pop,strain_total);
    SA = SA(end,:);
    SB = SB(end,:);
    IA = IA(end,:);
    IB = IB(end,:);
    N = SA+SB+IA+IB;
    
    % Remove extinct classes
    extinct = (N/sum(N))<exttol;
    strain_total = strain_total-sum(extinct);
    SA(extinct) = [];
    SB(extinct) = [];
    IA(extinct) = [];
    IB(extinct) = [];
    N(extinct) = [];
    index_current(extinct) = [];
    c_current(extinct) = [];
    
    % Update tracker
    SOCIALITY(ievol,index_current) = N/sum(N);
    DISPREV(ievol) = (sum(IA)+sum(IB))/sum(N);
    INFOPREV(ievol) = (sum(SB)+sum(IB))/sum(N);
    
    % Introduce rare mutant & update init_pop
    weightedprob = N/sum(N);
    cumsum1 = cumsum(weightedprob);
    r1 = rand*cumsum1(end);
    mutator_loc = (find(r1<cumsum1,1));
    mutator = index_current(mutator_loc);
    
    if(mutator==1) % Mutate up
        mutant = mutator+1;
    elseif(mutator==res0) % Mutate down
        mutant = mutator-1;
    else
        if(rand>0.5) % Mutate up
            mutant = mutator+1;
        else % Mutate down
            mutant = mutator-1;
        end
    end
    if(~ismember(mutant,index_current)) % New strain
        strain_total = strain_total+1;
        c_current(end+1) = C(mutant);
        index_current(end+1) = mutant;
        
        
        SA(end+1) = SA(mutator_loc)/10;
        SB(end+1) = SB(mutator_loc)/10;
        IA(end+1) = IA(mutator_loc)/10;
        IB(end+1) = IB(mutator_loc)/10;
    end
    
    % Update initial population
    init_pop = [];
    for i=1:strain_total
        init_pop = [init_pop,SA(i),SB(i),IA(i),IB(i)];
    end
end
N = SA+SB+IA+IB;

if(plotflag>0)
    figure(plotflag)
    clf
    a0=3;
    SOCIALITY0=log10(SOCIALITY);
    SOCIALITY0(SOCIALITY0<-a0)=-a0;
    SOCIALITY0=(SOCIALITY0+a0)/a0;
    subplot(1,3,1)
    imagesc(SOCIALITY0);set(gca,'ydir','normal')
    map=colormap('gray');
    map=flipud(map);
    colormap(map);
    ylabel('Evo time')
    xlabel('sociality, c')
    set(gca,'xtick',1:20:res0,'xticklabel',C(1:20:res0))
    subplot(1,3,2)
    plot(DISPREV)
    subplot(1,3,3)
    plot(INFOPREV)
end