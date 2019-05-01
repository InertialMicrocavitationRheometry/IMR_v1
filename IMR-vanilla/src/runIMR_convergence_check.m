% IMR solver convergence check
% ST = 0; numsim4: Neo-Hookean Kevin-Voigt; G = 7.69e3Pa; mu = 0.101 Pas;

%%
clear all; close all;
%%
% ====== Time duration ======
tspan = 2.25e-4; deltat = 1e-10; tList = 0:deltat:tspan; R0 = 225e-6; 

% ====== Material model ======
% G = 7.69e3; mu = 0.101; alpha = 1; model = 'fung'; simNo = 5;
% G = 9.5499e3; mu = 0.001; model = 'neoHook'; simNo = 4;
% G = 7.69e3; mu = 0.101; model = 'neoHook'; simNo = 4;
G = 2.12e3; mu = 0.118; model = 'neoHook'; simNo = 4;

expt = 1;
folderNamePrefix = ['E:\Jin\Franck\IMR-vanilla\data\numsim\',num2str(simNo),'\'];
fp = [folderNamePrefix,num2str(expt),'\'];
 
NT = 3000; IMRsolver_RelTolX = 1e-10;
tempfilename = ['numsim',num2str(simNo),'_NHKV_Exact_G',num2str(G),'_mu',num2str(mu),'_NT',num2str(NT),'_logRelTol',num2str(-log10(IMRsolver_RelTolX)),'.mat'];

load([fp tempfilename]); 
Rfit_pp_0 = interp1(t2',R2','pchip','pp');
RList0 = ppval(Rfit_pp_0,tList);

figure, plot(tList,RList0/R0);
% legend('Stiff PA: $G=7.69$kPa, $\mu=0.101$Pa$\cdot$s','Soft PA: $G=2.12$kPa, $\mu=0.118$Pa$\cdot$s','interpreter','latex')
a=gca; a.TickLabelInterpreter = 'latex';
set(gca,'fontsize',16); axis([0,2.25e-4,0,1]);

%%
NTList = [100,200,300,500,1000,1250,1500,1750,2000]; 
IMRsolver_RelTolXList = [1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9];
LSQErr = zeros(length(NTList),length(IMRsolver_RelTolXList));
ReqErr = LSQErr; ReqErrExact = LSQErr;
for tempi = 1:length(NTList)
    for tempj = 1:length(IMRsolver_RelTolXList)
        
        NT = NTList(tempi); IMRsolver_RelTolX = IMRsolver_RelTolXList(tempj);
        tempfilename = ['numsim',num2str(simNo),'_NHKV_Exact_G',num2str(G),'_mu',num2str(mu),'_NT',num2str(NT),'_logRelTol',num2str(-log10(IMRsolver_RelTolX)),'.mat'];

        load([fp tempfilename]); 
        Rfit_pp = interp1(t2',R2','pchip','pp');
        RList = ppval(Rfit_pp,tList);
        
        LSQErr(tempi,tempj) = norm(RList-RList0,2);
        ReqErrExact(tempi,tempj) = mean(RList(end-20:end)) - R0/4; % mean(RList0(end-20:end)) ;
        ReqErr(tempi,tempj) = mean(RList(end-20:end)) -  mean(RList0(end-20:end)) ;
        
    end
end

[NTGrid,RelTolXGrid] = ndgrid(NTList,IMRsolver_RelTolXList);

% ======
figure, surf(log10(NTGrid),log10(RelTolXGrid),log(LSQErr/sqrt(length(RList))/R0))
set(gca,'fontsize',16); axis tight;
xlabel('log10(NT)','interpreter','latex','Rotation',16); 
ylabel('log10(RelTol)','interpreter','latex','Rotation',-48);
%ylabel('log10(ode23tb solver "RelTol" parameter)','interpreter','latex','Rotation',-22);
title('log(RMS relative error of whole $R$-$t$ curve)','fontweight','normal','interpreter','latex');
a=gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
 
% ======
figure, surf(log10(NTGrid),log10(RelTolXGrid), log(abs(ReqErrExact)/R0));
set(gca,'fontsize',16);  axis tight;
xlabel('log10(NT)','interpreter','latex','Rotation',14); 
ylabel('log10(RelTol)','interpreter','latex','Rotation',-45);
%ylabel('log10(ode23tb solver "RelTol" parameter)','interpreter','latex','Rotation',-22);
title('log(Relative error of $R_{eq}$)','fontweight','normal','interpreter','latex');
a=gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

% ======
figure, surf(log10(NTGrid),log10(RelTolXGrid), log(abs(ReqErr)/R0));
set(gca,'fontsize',16);  axis tight;
xlabel('log10(NT)','interpreter','latex','Rotation',14); 
ylabel('log10(RelTol)','interpreter','latex','Rotation',-45);
%ylabel('log10(ode23tb solver "RelTol" parameter)','interpreter','latex','Rotation',-22);
title('log(Relative error of $R_{eq}$)','fontweight','normal','interpreter','latex');
a=gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';