%% Runfile for inertial microcavitation rheometry (IMR)

% Authors:
% Jin Yang: jyang526@wisc.edu
% Jon Estrada: jonathan_estrada@alumni.brown.edu; Brown Solid Mechanics, PhD '17
% Carlos Barajas: carlobar@umich.edu; Umich Mechanical Engineering, BS '16

% clear; close all;
warning('off','all')

%% Set up the file paths

%For CCV, file path is different
%general folder
% cd('/Users/yangjin/Documents/MATLAB/Franck/IMR-master');
% addpath(genpath('./data/')); % addpath(genpath('/gpfs/data/cf5/jbestrad/FL_Cav_Data/'));
%cd('E:\Jin\Franck\IMR-master');
%addpath(genpath('.\data\'));

%data folder
% fp = ['/Users/yangjin/Documents/MATLAB/Franck/IMR-master/data/numsim/',num2str(simNo),'/'];
% fp = './data/11kPa_PA/'; % fp = '/gpfs/data/cf5/jbestrad/FL_Cav_Data/160420/11kPa_PA/';
% fp = './data/1pt3kPa_PA/';
% fp = './data/collagen/1/'
% fp = ['./data/numsim/',num2str(simNo),'/'];
 
%simNo = 4; SampleNo = 100; GaussNo = 5; % fp = ['.\data\numsim\',num2str(simNo),'\'];
                               
%fp = '/gpfs/data/cf5/jbestrad/FL_Cav_Data/160420/water/';
%fp = '/gpfs/data/cf5/jbestrad/FL_Cav_Data/160511/collagen/1/';
%fp = '/gpfs/data/cf5/jbestrad/FL_Cav_Data/170403/Collagen/';
%fp = '/gpfs/data/cf5/jbestrad/FL_Cav_Data/170411/Collagen/20170302p1/';

%For Local comp usage (i.e. single runs)
%general folder
%addpath(genpath('V:\data\jbestrad\FL_Cav_Data\'));
%data folder
%fp = 'V:\data\jbestrad\FL_Cav_Data\160420\water\';
%fp = 'V:\data\jbestrad\FL_Cav_Data\170411\Collagen\20170302p1/';
%fp = 'C:\Users\Jon\Brown Google Drive\RESEARCH DATA\FL Cav Data\160420\11kPa_PA\';

%Load the file RofTdata.mat, which contains vars Rnew and t
%Rnew has size [num_expts num_video_frames]
% load([fp 'RofTdata.mat']);
% tempfilename  =  ['numsim',num2str(simNo),'_RofTdata_100.mat'];
% load([fp tempfilename]); Rnew = Rnew*1e6;
 

%There are 4 choices of model:
%linkv (linear Kelvin-Voigt)
%neoHook (neo-Hookean Kelvin-Voigt)
%sls (Standard Linear Solid)
%nhzen (neo-Hookean Standard Solid)
%model = 'neoHook'; % JY!!!
%savename = '190220_sweep_JY_';
%'170727_softsweep_coarse';%'170413_collagenKVcoarse';

%Parallel pool creation
% curCluster = parcluster('local');
% curCluster.NumWorkers = 8;
% saveProfile(curCluster);
% pool = parpool(8);
 
%Set the index choice of experiments +-

%expts = [12 14:19]; %water
%expts = [2,3,4,5,8,10,14,15,16,18,20,23,24]; %collagen
%expts = 1:1:4; PickNo = 7; % GaussNo = 5; SampleNo = 1; % numsim4
timestep = 1e-7; % by default
% expts = [1,2,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]; % soft

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== Comment G,mu,etc ======
%Set the range of G and mu (most conveniently done as powers of 10)
% temp_G_ooms_List = [10,100,500,1000,2000,3000,4000,5000,7000]; G_ooms = log10(temp_G_ooms_List); 
% G_ooms_step = 0.01; G_ooms = 3.2 : G_ooms_step : 3.8; % 3.9:0.02:4.3; % JY!!! % 1:0.2:5;%3.0:0.1:4.0; %soft PA %3.65:0.05:4.15 stiff PA
% temp_mu_ooms_List = [0.3, 0.2, 0.10, 0.05, 0.01, 0.001, 1e-4]; mu_ooms = log10(temp_mu_ooms_List); 
% mu_ooms_step = 0.02; mu_ooms = -2.7 : mu_ooms_step : -1.3; % -1.5:0.02:-0.7; % JY!!! 0.0533; %-1.4:0.05:-0.9;%[-inf -4:0.25:-1.25 -1.05];%[-1.65:0.05:-0.9]%-2.25:0.25:-0.5%[-inf -1.65:0.05:-0.9];
% ====== End of comment ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if expts==11
%     G_ooms = 2:0.1:4.3; mu_ooms =  -4 : 0.1: -1;
% elseif expts==1
%     G_ooms = 3.3:0.02:4.2; mu_ooms=-3:0.1:-1;
% elseif expts==6
%     G_ooms=3.5:0.02:4.3; mu_ooms=-2:0.05:-0.2;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% temp_alpha_ooms_List = [0.001,0.01,0.10,0.33,0]; alpha_ooms = log10(temp_alpha_ooms_List);
% alpha_ooms = -Inf;
% temp_lambda_nu_oooms_List = [0, 3.3e-3, 1e-2, 3.3e-2, 1e-1, 1e0, 1e1, 1e2]; lambda_nu_ooms = log10(temp_lambda_nu_oooms_List);
% lambda_nu_ooms = 0;
% G1_ooms = Inf;
%Note, mu/G1 should go to 0 in the limit, so G1 for K-V models should be infinite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% whetherToSolveP = 0; DefineLocsMaxPOrNot = 0; DefineTimetolListOrNot = 0; PlotDefineTimetolListOrNot = 0; 

% TimeRes = 10;


%% %%%%%%%%%%%%%%% First part %%%%%%%%%%%%%%%%%%
%% Run IMR to solve {P,C,T}
if whetherToSolveP > 0
for tempk = 1:length(expts)
    
    expt = expts(tempk)
    fp = [folderNamePrefix,num2str(expt),'/'];
    %cd([fp]);
    %fp = ['./data/numsim/',num2str(simNo),'/',num2str(expt),'/'];
    % fp = ['E:\Jin\Franck\IMR-master\data\numsim\',num2str(simNo),'\',num2str(expt),'\'];
    % tempfilename  =  ['numsim',num2str(simNo),'_expt',num2str(expt),'_RofTdata_sample',num2str(SampleNo),'_gauss',num2str(GaussNo),'.mat'];
    tempfilename  =  [fileNamePrefix,'.mat'];
    load([fp tempfilename]); Rnew  =  Rnew*1e6;
    
    %Set time duration for the simulation (s)
    %tspan = 2.25E-4; %1.3E-4;
    %allRmax = max(Rnew,[],2);
    
    % mkdir([fp num2str(expt)]);
    
     % tempfilename  =  ['Rfit_numsim',num2str(simNo),'_exp',num2str(expt),'_sample',num2str(SampleNo),'_gauss',num2str(GaussNo),'.mat']; load( tempfilename);  
    tempfilename  =  ['Rfit_numsim',num2str(simNo),'_exp',num2str(expt),'_',num2str(TimeRes),'.mat']; 
    load(tempfilename);
    
    %%
    %soln_mx = cell(length(G_ooms),length(mu_ooms),length(alpha_ooms),length(lambda_nu_ooms));
    % JY!!! % soln_mx = cell(length(G_ooms),length(mu_ooms),length(G1_ooms));
    % P_inf = 101325; % (Pa) Atmospheric Pressure
    % rho = 998.2; % JY!!! % rho = 998.2; % (Kg/m^3) 1060 Material Density
    % Uc = sqrt(P_inf/rho); % Characteristic velocity
    % 
    % Tgrad = 1; %1 Temperature transfer bt bubble and material
    % Cgrad = 1; %1;  %1 Vapor-non-condensible gas diffusion
    % Tmgrad = 0; %0 Off means cold liquid assumption
    
    % Load parameters: 
    Pmt = IMRcall_parameters(1,1,1,1); % The first input as "1" doesn't matter.
    P_inf = Pmt(19); T_inf = Pmt(20); rho = Pmt(23); Uc = Pmt(24); ST = Pmt(28);
    
    Tgrad = Pmt(25); Cgrad = Pmt(26); Tmgrad = Pmt(27);
       
    %%
    R0 = Rfit(1); %[R0,t0] =  calcR0(Rnew(1,:)*1E-6,t); %need to get the inital radius from a fit (R)
    eqR = mean(Rnew(1, end-20:end )) *10^-6; % Solve for equilibirium R_inf
    R_eq = eqR/R0;
                 
    %% Bisection method to get initial partial pressure of the non-condensible gas
     P_guess = (P_inf+2*ST./eqR-Pvsat(T_inf)).*((eqR/R0).^3);
    
%     % tic
%     % Initial guess (purely empirical) and range (e.g. ~226 for
%     % 11kPa, 92 for 1.3kPa, 20.6 for water)
%     % JY!!! which for collagen?
%     P_guess = 200;
%     P = [8, 400];
%     
%     % temp G, mu, G1
%     G = 2970; mu = 0.01; G1 = Inf;
% 
%     R_eq_lims = [IMRCalc_Req(R0,Tgrad,Cgrad,P(1),G,G1,mu), IMRCalc_Req(R0,Tgrad,Cgrad,P(2),G,G1,mu)];
%     R_eqf = IMRCalc_Req(R0,Tgrad,Cgrad,P_guess,G,G1,mu);
%     error = 1; BisectionIter = 1;
%     while abs(error) > 0.00001 && BisectionIter<1e6
%         if R_eq > R_eqf
%             P(1) = P_guess;
%             R_eq_lims(1) = IMRCalc_Req(R0,Tgrad,Cgrad,P_guess,G,G1,mu);
%         else
%             P(2) = P_guess;
%             R_eq_lims(2) = IMRCalc_Req(R0,Tgrad,Cgrad,P_guess,G,G1,mu);
%         end
%         P_guess = mean(P);
%         R_eqf  = IMRCalc_Req(R0,Tgrad,Cgrad,P_guess,G,G1,mu);
%         error = abs(R_eq-R_eqf);
%         BisectionIter = BisectionIter+1;
%     end
%     % toc
%     if abs(error) > 1e-5 || P_guess < 8 || P_guess > 400
%         disp(['Fail to find initial partial pressure of gas in ', num2str(toc),' seconds.']);
%     else   
%         disp(['P_guess = ',num2str(P_guess)]);
%         disp(['Find initial partial pressure of gas in ', num2str(toc),' seconds.']);
%     end
 %%
    NT = 500; % Mesh points in bubble, resolution should be >=500
    NTM = 10; % Mesh points in the material, should be 10
    Pext_type = 'IC'; %'IC' for Flynn, 'ga' for gaussian bubble growth

    if strcmp(Pext_type,'IC')
        Pext_Amp_Freq = [P_guess; 0];%[226; 0]; % Tune first number to match equi radii
    elseif strcmp(Pext_type,'ga')
        Pext_Amp_Freq = [P_guess; dt; tw;];
    end

    disptime = 0; % 1 = Displays time to complete simulation
    Dim = 1;  % 1 = displays results in dimensional form
    comp = 1; % 0 uses Rayleigh-Plesset, 1 uses Keller-Miksis

    %%                
    if strcmp(Pext_type,'ga')
        Rmax = R0;
        R0 = eqR;
    end
    
    %% ===== IMRsolver2 solve {P,C,T} with measured R, Rdot, and Rddot =====
     
    
    P0Old = -1; P0New = 0; P0Iter = 0; eqPTheory = P_inf+2*ST/eqR; % based on bubble equilibrium state. 
    P_upperlimit = 400; P_lowerlimit = 80; % JY!!! pay attention to the assigned limit values of P.
     
    while (abs(P_upperlimit-P_lowerlimit)/P_upperlimit > 1e-3) && (P0Iter < 100)
        P0Iter = P0Iter+1;
        % ====== IMR decoupled system solver + tune P_initial ======
        tic;
        [t2_num, P, T, C, tdel, Tdel, Cdel] = IMRsolver2_DecouplePCT(...
            model, tspan, R0, NT, NTM,...
            Pext_type, Pext_Amp_Freq, disptime, Tgrad, Tmgrad, Cgrad, Dim, comp,...
            Rfit_pp,Rfit_pp_d1,Rfit_pp_d2,t2Start,eqR);
        toc
        if whetherToSolveP == 1
            eqP = mean(P(end-20:end));
            if eqP > eqPTheory % we should decrease P_initial
                P_upperlimit = Pext_Amp_Freq(1);
                Pext_Amp_Freq(1) = Pext_Amp_Freq(1)*0.5 + P_lowerlimit*0.5;
            else % eqP < eqPTheory, we should increase P_initial
                P_lowerlimit = Pext_Amp_Freq(1);
                Pext_Amp_Freq(1) = Pext_Amp_Freq(1)*0.5 + P_upperlimit*0.5;
            end
            disp(['IterNo: ',num2str(P0Iter),'; Current val: ',num2str(Pext_Amp_Freq(1)),'; Upper limit: ',num2str(P_upperlimit),'; Lower limit: ',num2str(P_lowerlimit)]);
        else % whetherToSolveP == 2
            P0Iter = 100;
        end
        % ==========================================================
    end
    file_name = ['numsim',num2str(simNo),'_IMRSolver2_exp',num2str(expt),'_',num2str(TimeRes)];
    
    
    cd([fp]);
    save([file_name '.mat'],'t2_num','C','T','P','tdel','Tdel','Cdel','Pext_Amp_Freq');
    %cd('../../../');
    cd('/Users/yangjin/Documents/MATLAB/Franck/IMR_code');
    
    % *********** Modify P errors at t_inf ***************
    % eqP = median(P(end-20:end)); exactP = P_inf+2*0.056/eqR; DeltaP = exactP-eqP;
    % % P = P+(exp(log(1.2)/(1.5e-4)*t2_num)-1)/(1.2-1)*DeltaP; 
    % [roweqP,~]=find(t2_num>1.5e-4); P(roweqP)=exactP;

end
end


%% Run IMR
if whetherToSolveP == 0
    
    
matPropVarList = cell(length(expts),1);

for tempk = 1:length(expts)
    
     expt = expts(tempk); disp(['====== Expt #: ',num2str(expt),' ======']);
    
    %fp = ['./data/numsim/',num2str(simNo),'/',num2str(expt),'/'];
    fp = [folderNamePrefix,num2str(expt),'/'];
    % fp = ['E:\Jin\Franck\IMR-master\data\numsim\',num2str(simNo),'\',num2str(expt),'\'];
    % tempfilename  =  ['numsim',num2str(simNo),'_expt',num2str(expt),'_RofTdata_sample',num2str(SampleNo),'_gauss',num2str(GaussNo),'.mat'];
    %tempfilename  =  ['numsim',num2str(simNo),'_RofTdata_',num2str(TimeRes),'.mat'];
    tempfilename  =  [fileNamePrefix,'.mat']; 
    load([fp tempfilename]); Rnew = Rnew*1e6;
    
    %Set time duration for the simulation (s)
    %tspan = 2.25E-4; %1.3E-4;
    %allRmax = max(Rnew,[],2);
    
    % mkdir([fp num2str(expt)]);
    %cd([fp]);
    
    %%
    % soln_mx = cell(length(G_ooms),length(mu_ooms),length(alpha_ooms),length(lambda_nu_ooms));
    % JY!!! % soln_mx = cell(length(G_ooms),length(mu_ooms),length(G1_ooms));
    % P_inf = 101325; % (Pa) Atmospheric Pressure
    % rho = 998.2; % JY!!! % rho = 998.2; % 1060 (Kg/m^3) Material Density
    % Uc = sqrt(P_inf/rho); % Characteristic velocity
    % 
    % Tgrad = 1; %1 Temperature transfer bt bubble and material
    % Cgrad = 1; %1;  %1 Vapor-non-condensible gas diffusion
    % Tmgrad = 0; %0 Off means cold liquid assumption
    
    % Load parameters: 
    Pmt = IMRcall_parameters(1,1,1,1); % The first input as "1" doesn't matter.
    P_inf = Pmt(19); T_inf = Pmt(20); rho = Pmt(23); Uc = Pmt(24); ST = Pmt(28);
    
    Tgrad = Pmt(25); Cgrad = Pmt(26); Tmgrad = Pmt(27);
    
       
    %%
    [R0,t0] = calcR0(Rnew(1,:)*1E-6,t); %need to get the inital radius from a fit (R)
    eqR = mean(Rnew(1, end-20:end )) *10^-6; % Solve for equilibirium R_inf
    R_eq = eqR/R0;
    
    %% Bisection method to get initial partial pressure of the non-condensible gas
    % tic
    % Initial guess (purely empirical) and range (e.g. ~226 for
    % 11kPa, 92 for 1.3kPa, 20.6 for water)
    % JY!!! which for collagen?
%     P_guess = 200;
%     P = [8, 400];
%     
%     % temp G, mu, G1
%     G = 2970; mu = 0.01; G1 = Inf;
% 
%     R_eq_lims = [IMRCalc_Req(R0,Tgrad,Cgrad,P(1),G,G1,mu), IMRCalc_Req(R0,Tgrad,Cgrad,P(2),G,G1,mu)];
%     R_eqf = IMRCalc_Req(R0,Tgrad,Cgrad,P_guess,G,G1,mu);
%     error = 1; BisectionIter = 1;
%     while abs(error) > 0.000001 && BisectionIter<1e6
%         if R_eq > R_eqf
%             P(1) = P_guess;
%             R_eq_lims(1) = IMRCalc_Req(R0,Tgrad,Cgrad,P_guess,G,G1,mu);
%         else
%             P(2) = P_guess;
%             R_eq_lims(2) = IMRCalc_Req(R0,Tgrad,Cgrad,P_guess,G,G1,mu);
%         end
%         P_guess = mean(P);
%         R_eqf  = IMRCalc_Req(R0,Tgrad,Cgrad,P_guess,G,G1,mu);
%         error = abs(R_eq-R_eqf);
%         BisectionIter = BisectionIter+1;
%     end
%     % toc
%     if abs(error) > 1e-6 || P_guess < 8 || P_guess > 400
%         disp(['Fail to find initial partial pressure of gas in ', num2str(toc),' seconds.']);
%     else   
%         P_guess % disp(['Find initial partial pressure of gas in ', num2str(toc),' seconds.']);
%     end
 
 %% ===== Load IMRsolver2 solved {P,C,T} =====
    tempfilename = ['Rfit_numsim',num2str(simNo),'_exp',num2str(expt),'_',num2str(TimeRes),'.mat']; 
    load(tempfilename);  
    
    %****** This part already solved before ******
    file_name = ['numsim',num2str(simNo),'_IMRSolver2_exp',num2str(expt),'_',num2str(TimeRes),'.mat'];
    load(file_name); 
   
    
    
    NT = 500; % Mesh points in bubble, resolution should be >=500
    NTM = 10; % Mesh points in the material, should be 10
    Pext_type = 'IC'; %'IC' for Flynn, 'ga' for gaussian bubble growth

%     if strcmp(Pext_type,'IC')
%         Pext_Amp_Freq = [P_guess; 0];%[226; 0]; % Tune first number to match equi radii
%     elseif strcmp(Pext_type,'ga')
%         Pext_Amp_Freq = [P_guess; dt; tw;];
%     end

    disptime = 0; % 1 = Displays time to complete simulation
    Dim = 1;  % 1 = displays results in dimensional form
    comp = 1; % 0 uses Rayleigh-Plesset, 1 uses Keller-Miksis
                 
   
	%% Compute R2_num, R2dot, R2ddot Pdot 
    
    %R2_num = ppval(Rfit_pp,t2_num+t2Start); 
    %figure; hold on; plot(t2_num+t2Start,R2_num );
    tspanEndIndex = min(find(t2'>t2_num(end)));
    if isempty(tspanEndIndex), tspanEndIndex = length(t2); end
    figure; plot(t2(1:tspanEndIndex),Rfit(1:tspanEndIndex)); set(gca,'fontsize',18);  axis([0,tspan,0,max(Rfit)]);
    xlabel('Time($s$)','interpreter','latex'); ylabel('Fitted bubble radius $R$(m)','interpreter','latex');
    a=gca; a.TickLabelInterpreter = 'latex';

    % R2dot = ppval(Rfit_pp_d1,t2_num+t2Start); % figure; plot(t2_num,R2dot,'ro');
    %deltat = timestep; t2_temp = 0+t2Start : deltat : tspan+t2Start ;
    %R2dot_temp = ppval(Rfit_pp_d1,t2_temp);  R2dot_temp(1) = 0;
    %figure; plot(t2_temp,R2dot_temp);
    figure; plot(t2(1:tspanEndIndex),d1Rfit(1:tspanEndIndex)); set(gca,'fontsize',18);  axis([0,tspan,min(d1Rfit),max(d1Rfit)]);
    xlabel('Time($s$)','interpreter','latex'); ylabel('$Fitted {\partial}{R}/{\partial}{t}$','interpreter','latex');
    set(gca,'fontsize',18); a=gca; a.TickLabelInterpreter = 'latex';
    
    figure; plot(t2(1:tspanEndIndex),d2Rfit(1:tspanEndIndex)); set(gca,'fontsize',18);  axis([0,tspan,min(d2Rfit),max(d2Rfit)]);
    xlabel('Time($s$)','interpreter','latex'); ylabel('$Fitted {\partial^2}{R}/{\partial}{t^2}$','interpreter','latex');
    set(gca,'fontsize',18); a=gca; a.TickLabelInterpreter = 'latex';
    

%     %PickNo = 7; % prompt='How many peaks in radius?'; PickNo = input(prompt); 
%     tNaNs = zeros(3*PickNo,2);
% 
%     % ****** Filter R2dot ******
%     % for tempjj = 1:5, R2dot_tempNew  = conv(R2dot_temp,[1 2 1]/4);  R2dot_temp=R2dot_tempNew(2:end-1); end
%     for tempjj = 1:30, R2dot_tempNew = conv(R2dot_temp,[1 2 5 2 1]/11); R2dot_temp = R2dot_tempNew(3:end-2); end
%     % hold on; plot(t2_temp,R2dot_temp);
%     R2dotfit_pp = interp1(t2_temp,R2dot_temp,'pchip','pp');
%     R2dot = ppval(R2dotfit_pp,t2_num); hold on; plot(t2_num,R2dot);
% 
%     % ****** Find peaks of R2dot ******
%     [pksUp_R2dot, locsUp_R2dot] = findpeaks( (R2dot), t2_num, 'MinPeakDist',1e-6);
%     [pksDown_R2dot, locsDown_R2dot] = findpeaks( (-R2dot), t2_num, 'MinPeakDist',1e-6);
%     % To vanish points between [locsDown_R2dot,locsUp_R2dot].
%     for tempi = 1:PickNo, tNaNs(tempi,1:2) = [locsDown_R2dot(tempi), locsUp_R2dot(tempi)]; end
% 
%     % ****** Compute R2ddot ******
%     imgradientMatrix = [-1/60 3/20 -3/4 0 3/4 -3/20 1/60]';
%     R2ddot_temp = imfilter(R2dot_temp,imgradientMatrix')/deltat; 
%     figure; plot(t2_temp,R2ddot_temp); axis([0,tspan,-3/deltat,3/deltat]);
% 
%     % ****** Filter R2ddot ******
%     % for tempjj = 1:30 , R2ddot_tempNew  = conv(R2ddot_temp,[1 2 1]/4);  R2ddot_temp=R2ddot_tempNew(2:end-1); end
%     for tempjj = 1:30 , R2ddot_tempNew  = conv(R2ddot_temp,[1 2 5 2 1]/11); R2ddot_temp=R2ddot_tempNew(3:end-2); end
%     R2ddotfit_pp = interp1(t2_temp,R2ddot_temp,'pchip','pp');
%     R2ddot = ppval(R2ddotfit_pp,t2_num); hold on; plot(t2_num,R2ddot);  
% 
%     % ****** Find peaks of R2ddot ******
%     [pksUp_R2ddot, locsUp_R2ddot] = findpeaks( (R2ddot), t2_num, 'MinPeakDist',1e-6);
%     [pksDown_R2ddot, locsDown_R2ddot] = findpeaks( (-R2ddot), t2_num, 'MinPeakDist',1e-6);
%     % To vanish points between [locsDown_R2ddot(2*i-1),locsDown_R2ddot(2*i)].
%     % JY!!! haven't finished.
%     tempindex=[1:2*PickNo];
%     for tempi = 1:PickNo  
%         tNaNs(PickNo+tempi,1:2) = [locsDown_R2ddot(tempindex(2*(tempi)-1)),locsDown_R2ddot(tempindex(2*(tempi)))];
%     end

    % ****** Compute Pdot ******
%     figure; plot(t2_num+t2Start,P); set(gca,'yscale','log'); axis auto; set(gca,'fontsize',20); title('P-t','fontweight','normal');
%     deltat = timestep; t2_temp = 0+t2Start:deltat:tspan+t2Start;
%     Pfit_pp = interp1(t2_num(1:end)'+t2Start,P(1:end)','pchip','pp');
%     Pfit_pp_d1 = fnder(Pfit_pp,1); Pdot = ppval(Pfit_pp_d1,t2_num+t2Start); 
%     figure; plot(t2_num+t2Start,Pdot); set(gca,'yscale','log'); axis auto; set(gca,'fontsize',20); title('Pdot-t','fontweight','normal');

 
    Pfit_pp = interp1(t2_num(1:end)'+t2Start,P(1:end)','pchip','pp');
    Pfit = ppval(Pfit_pp,t2);
    Pfit_pp_d1 = fnder(Pfit_pp,1); 
    Pfitdot = ppval(Pfit_pp_d1,t2); 
     
    figure; plot(t2(1:tspanEndIndex),Pfit(1:tspanEndIndex)); set(gca,'fontsize',18);  axis([0,tspan,min(Pfit),max(Pfit)]); set(gca,'yscale','log');
     xlabel('Time($s$)','interpreter','latex'); ylabel('Bubble insie pressure $p$ (Pa)','interpreter','latex');
       a=gca; a.TickLabelInterpreter = 'latex';
    
    %figure; plot(t2(1:tspanEndIndex),Pfitdot(1:tspanEndIndex)); set(gca,'fontsize',18);  axis([0,tspan,min(Pfitdot),max(Pfitdot)]); set(gca,'yscale','log');
    %xlabel('Time($s$)','interpreter','latex'); ylabel('Fitted ${\partial}{p}/{\partial}{t}$ (Pa/s)','interpreter','latex');
    %a=gca; a.TickLabelInterpreter = 'latex';
    
    
    % ****** Find peaks of Pdot ******
    %[pksUp_P, locsUp_P] = findpeaks( (P), t2_num+t2Start, 'MinPeakDist', 5e-6);
    [pksMaxP, tspanLocsMaxP] = findpeaks( (Pfit), t2 , 'MinPeakDist', 5e-6);
    % Not totally finished. Need to modify code to do this automatically.
    % locsUp_P=1e-3*[0.0202, 0.0394, 0.0533, 0.0648, 0.0750, 0.0846, 0.0939, 0.1032,  0.1125, 0.1218]'; % numsim4
%     if DefineLocsMaxPOrNot==0
%         if (expt==1) && (GaussNo == 1) && (SampleNo == 4)
%             locsUp_P = 1e-5*[2.024 3.799 5.05 6.075 6.975 7.816 8.774];
%         end
%     end
%     if expt==1
%         locsUp_P = 1e-3*[0.0738, 0.1073, 0.1296, 0.1481, 0.1666, 0.1851]; % stiff_exp1
%     elseif expt==3
%         locsUp_P = 1e-3*[0.07582 0.1110 0.1369 0.1554 0.1772 0.1960];
%     elseif expt==7
%         locsUp_P = 1e-3*[0.0772 0.1110 0.1309 0.1517 0.1703  0.1905]; 
%     elseif expt==14
%         locsUp_P = 1e-3*[0.07395 0.1073 0.1296 0.1481 0.1666 0.1851]; 
%     elseif expt==16
%         locsUp_P = 1e-3*[0.0766 0.1073 0.1294 0.1480 0.1666 0.1828]; 
%     elseif expt==20
%         locsUp_P = 1e-3*[0.0738 0.1108 0.1331 0.1517 0.1735 0.1926]; 
%     else
%         locsUp_P = locsUp_P(1:6);  % 1e-5*[3.549, 6.881, 9.101, 10.95, 12.8, 14.66 ]; % stiff_exp1
%     end
    
    if DefineLocsMaxPOrNot==1
       fprintf('First three collapses where p values are large happens at time:  \n') 
            disp(num2str(tspanLocsMaxP(1:PickNo)));
            fprintf('Double these these time values are correct, and press "Enter".  \n') 
            pause;
    end
  
    % locsUp_P = locsUp_P(1:10);
    %[pksDown_Pdot, locsDown_Pdot] = findpeaks( (-Pfitdot), t2_num+t2Start, 'MinPeakDist',1e-6);
    % To vanish points between [locsDown_R2dot,locsUp_R2dot].
    
    %% ===== Compute A1,A2,B1,B2,E =======
    % comment this section
    for tempii = 1
%     [A1,A2,B1,B2,E] = IMRsolver2_Comp_ABE(model,G,G1,mu,tspan,R0,NT,NTM, ...
%         Pext_type,Pext_Amp_Freq,disptime,Tgrad,Tmgrad,Cgrad,Dim,comp, t2_num,R2_num,R2dot,R2ddot,P,Pdot,T,C,eqR);
%      
%     figure; plot(t2_num,A1); axis auto; set(gca,'fontsize',20); title('A1(t)-t','fontweight','normal');
%     figure; plot(t2_num,A2); axis([0,tspan,-10,30]);  set(gca,'fontsize',20); title('A2(t)-t','fontweight','normal');
%     figure; plot(t2_num,B1); axis([0,tspan,-8,15]);  axis auto; set(gca,'fontsize',20);  title('B1(t)-t','fontweight','normal');
%     figure; plot(t2_num,B2); axis([0,tspan,-1e2,2e2]); set(gca,'fontsize',20); title('B2(t)-t','fontweight','normal');
%     figure; plot(t2_num,E); axis([0,tspan,-1,2]); % set(gca,'fontsize',20); title('E(t)-t','fontweight','normal');
%     E_temp=E; for tempjj = 1:10 , E_tempNew  = conv(E_temp,[1 2 1]/4);  E_temp=E_tempNew(2:end-1); end
%         figure; plot(t2_num,E_temp); axis([0,tspan,-1,2]);  set(gca,'fontsize',20); title('E(t)-t','fontweight','normal');
% 
%     rowRemove = []; 
%     for tempi = 1:size(tNaNs,1)
%     [rowRemovetemp,~] = find(t2_num>tNaNs(tempi,1)-1e-6 & t2_num<tNaNs(tempi,2)+1e-6);
%     rowRemove = unique([rowRemove;rowRemovetemp]);
%     end
%     B2Neat = B2; B2Neat(rowRemove) = NaN; ENeat = E_temp; ENeat(rowRemove) = NaN; 
%     figure; plot(t2_num,ENeat); axis([0,tspan,-1,2]);  set(gca,'fontsize',20); title('ENeat(t)-t','fontweight','normal');
%     figure; plot(t2_num,B2Neat); axis([0,tspan,-1e2,2e2]); set(gca,'fontsize',20); title('B2Neat(t)-t','fontweight','normal');
% 
%     % Remove points where B2<1e-4
%     [rowRemovetemp,~] = find(abs(B2)<1e-1);
%     rowRemove=unique([rowRemove;rowRemovetemp]);
%     % Remove points where before first peak in E(t)
%     [pksUp_ENeat, locsUp_ENeat] = findpeaks( (ENeat), t2_num, 'MinPeakDist',1e-6);
%     [rowRemovetemp,~] = find(t2_num<locsUp_ENeat(1));
%     rowRemove=unique([rowRemove;rowRemovetemp]);
% 
%     % ****** Find peaks of B2(t) ******
%     % Delete points where abs(B2(t))<tol 
%     % ReTildeGuess = -E ./B2 ; figure; plot(t2_num,ReTildeGuess,'r.-');
%     ReTildeGuess = -ENeat./B2Neat; ReTildeGuess(rowRemove)=NaN; 
%     figure; plot(t2_num,ReTildeGuess,'r.-'); title('-E(t)/B2(t)','fontweight','normal'); set(gca,'fontsize',20);
%     % For mu=0.1; mean(ReTildeGuess(36:55)); std(ReTildeGuess(36:55))
%     % For mu = 0.01; mean(ReTildeGuess(39:48)); std(ReTildeGuess(39:48))
% 
%     % *****************************************************************
%     % *** This part originally to solve realtime Re# and Ca# ***
%     % deltat = timestep; t4 = t0:deltat:tspan+t0;
%     % Rfit_pp_num = interp1(t0+t2_num,R2_num,'pchip','pp'); 
%     % Rfit_num = ppval(Rfit_pp_num,t4); d = Rfit_num;
%     % Ts = (t4(end)-t0)/length(d); Fs = 1/Ts; Fn = Fs/2; tv_num = (0:length(d)-1)*Ts;
%     % st = find(d > 1e-5, 1, 'first'); tvst_num = tv_num(st);
%     % d = d(st:end); tv_num = tv_num(st:end); L = length(d); % figure; plot(tv_num,d);
%     % [pks, locs] = findpeaks( (d), tv_num, 'MinPeakDist',1e-6);
%     % pks = [d(1) pks]; locs = [0 locs]; % hold on; plot(locs,pks); 
%     % peakVal_num = pks(1:min([10,length(pks)]))';
%     % peakLoc_num = locs(1:min([10,length(pks)]))'+tvst_num ;
%     % 
%     % ReTildefit_pp_num = interp1(t0+t2_num,ReTildeGuess,'pchip','pp'); 
%     % ReTildefit = ppval(ReTildefit_pp_num,peakLoc_num);  
%     % 
%     % hold on; plot(peakLoc_num,ReTildefit,'ro'); axis([0,tspan,-1,1])
%     % % hold on; plot(peakLoc_num,0,'ro');
% 
%     % [t3,ReTilde] = IMRsolver2_Comp_ReTilde(G,G1,mu,tspan,R0,t2_num,B1,B2,E, ReTildeGuess(1));
%     % 
%     % P_inf = 101325;  Ca = P_inf/G;  CaTilde = 1/Ca;
%     % Uc = sqrt(P_inf/rho); Re = P_inf*R0/(mu*Uc); 
%     % 
%     % mu_num = ReTilde*rho* Uc*R0; figure; plot(t3,mu_num);
    end

    %% Tune material parameters
    close all; tic;
    % JY!!! dangerous
    % PickNo = 7; % already defined before
    % prompt='How many peaks in radius?'; PickNo = input(prompt); tic; % 8 for numsim4
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for i = 1 :length(G_ooms)
    % for j = 1 :length(mu_ooms) %mu
        % for k = 1 :length(alpha_ooms) % for j = 1:length(mu_ooms)
        % for l = 1 % 1:length(lambda_nu_ooms) % for i=1:length(G_ooms)
                
            % % waitbar((j-1+(i-1)*length(mu_ooms))/length(G_ooms)/length(mu_ooms)); % [i,j,k,l]
            % waitbar((k-1+(j-1+(i-1)*length(mu_ooms))*length(alpha_ooms))/length(alpha_ooms)/length(G_ooms)/length(mu_ooms)); % [i,j,k,l]

            % soln_mx{i,j,k,l} = struct('G',10^G_ooms(i),'mu',10^mu_ooms(j),'G1',10^G1_ooms(1),'alpha',10^alpha_ooms(k),'lambda_nu',10^lambda_nu_ooms(l),...
            %    'tcs_star',[],'Rratios',[],'tmaxs_star',[],'t2_num',[],'R2_num',[],'U',[],'P',[],'J',[],'T',[],'C',[],'Cdel',[],'tdel',[],'Tdel',[]);

            % G = soln_mx{i,j,k,l}.G;
            % mu = soln_mx{i,j,k,l}.mu;
            % G1 = soln_mx{i,j,k,l}.G1;
            % alpha = soln_mx{i,j,k,l}.alpha;
            % lambda_nu = soln_mx{i,j,k,l}.lambda_nu;
            
            % ====== expt 1 ======
            if PlotDefineTimetolListOrNot==1
                %G = 2.97e3; mu=0.01; 
                PlotExact=1;
                % G = 2.0959e3; mu=0.0110; PlotExact=0; openfig('numsim4_expt1_100.fig');
                % G = 2.8188e3; mu=0.0032; PlotExact=0; % openfig('numsim4_expt1_sample3_gauss5.fig');
                % ====== expt 6 ======
                % G = 7.690e3; mu=0.1; PlotExact=1;
                % G = 8.2773e3; mu=0.1005; PlotExact=0; openfig('numsim4_expt6_sample4_gauss5.fig');
                % G = 8.371e3; mu=0.0849; PlotExact=0; openfig('numsim4_expt6_sample3_gauss5.fig');
            end
             
           
           %%
           d_exp = Rfit(1:tspanEndIndex); %tNew =  0+t2Start:timestep:tspan+t2Start; Rfit_exp = ppval(Rfit_pp,tNew); d_exp =  Rfit_exp;
           Ts = tspan/length(d_exp); Fs = 1/Ts; Fn = Fs/2; tv_exp = (0:length(d_exp)-1)*Ts+t2Start; % figure; plot(tv,d,'rs');
           st = find(d_exp > 1e-5, 1, 'first'); tvst_exp = tv_exp(st);
           d_exp = d_exp(st:end); tv_exp = tv_exp(st:end); L = length(d_exp);
           [pks, locs, widthpks] = findpeaks( (d_exp), tv_exp, 'MinPeakDist',1E-6);
           pksMaxR = pks(1:PickNo)'; % pksMaxR = [pksMaxR;eqR];
           tspanLocsMaxR = locs(1:PickNo)'; tspanLocsMaxROld = tspanLocsMaxR; % tspanLocsMaxR=[tspanLocsMaxR;tspan];
            if tspanLocsMaxR(1)-t2Start>1e-5
                tspanLocsMaxR=[t2Start;tspanLocsMaxR];
                pksMaxR=[d_exp(1);pksMaxR];
            end
            for tempi = 2:PickNo 
                [~,temptCen] = find(tv_exp > tspanLocsMaxR(tempi)); temptCen = min(temptCen)-1;
                tempwidth = 0.4*widthpks(tempi)/timestep;
                tempt = tv_exp(temptCen-tempwidth:temptCen+tempwidth);
                tempd_exp = d_exp(temptCen-tempwidth:temptCen+tempwidth);
                % figure; plot(tempt,tempd_exp,'r.');
                tempp = polyfit(tempt,tempd_exp,2);
                tspanLocsMaxR(tempi) = -tempp(2)/2/tempp(1);
                pksMaxR(tempi) = (4*tempp(1)*tempp(3)-tempp(2)^2)/4/tempp(1);
            end
            
            % *************** plot ******************
            if DefineTimetolListOrNot==1
                q = [locs' ones(size(locs'))]\log(abs(pks))';          % Initial Parameter Estimaes
                fitfcn = @(b,t) b(1) .* exp(b(2).*t) + b(3);           % Fit Function
                SSECF = @(b) sum((pks' - fitfcn(b,locs')).^2);
                [B_exp,SSE] = fminsearch(SSECF, [pks(1); q(1); 0]);
                figure; plot(t,Rnew(1,:)*1e-6,'rs'); hold on; 
                %plot(tv_exp , fitfcn(B_exp,tv_exp),'k-.'); hold on;
                plot(tspanLocsMaxR,pksMaxR,'rd') ;hold on;
                plot(tv_exp,d_exp,'r'); hold off; 
                grid off; grid; set(gca,'fontsize',18); grid minor; 
                for tempi = 1:PickNo
                    hold on; plot(tspanLocsMaxP(tempi)*ones(1,101), linspace(0,2.5e-4,101),'k--');
                end
                % title(['G=',num2str(10^G_ooms),' \mu=',num2str(10^mu_ooms),', \alpha=',num2str(alpha)],'FontWeight','Normal');
                title('Experiment data','FontWeight','Normal'); set(gca,'fontsize',18)
                axis([0,tspan,0,1.12*R0]);
            end
            
            
            % *************** plot ******************
            timetol = timestep; timetolList = 1e-6*ones(max(expts),PickNo*2-1); % by default
            % **************** Tune timetolList *******************
%             if DefineTimetolListOrNot == 0
%                 if (expt==1) && (GaussNo == 5) && (SampleNo == 100) 
%                     % timetolList(1,:) = 1e-6*[];
%                 end
%             end
            
            %%
           
            % cd('E:\Jin\Franck\IMR-master\');
                
            % alpha_ooms = -Inf; lambda_nu_ooms = 0; G1_ooms = Inf; 
            % alpha = 10^alpha_ooms; lambda_nu = 10^lambda_nu_ooms; G1 = 10^G1_ooms;
            %NH_KV_varList0 = [2, -2];

            %lb = log10([500, 5e-3]);
            %ub = log10([1.5e5, 0.1]);
            % options = optimset('Display','off' );
            % options = optimset('Display','iter','PlotFcns','optimplotfval');
            options = optimset('Display','iter','PlotFcns','optimplotfval','TolX',NelderMead_TolX);
          
%             [NH_KV_varList_best,fval,exitflag,output] = simulannealbnd( @(NH_KV_varList) ...
%                 funIMR2LSQErr(expt,PickNo,model,NH_KV_varList,G1_ooms,alpha_ooms,lambda_nu_ooms,tspanLocsMaxR,pksMaxR,NT,NTM,...
%                 Pext_type,Pext_Amp_Freq,disptime,Tgrad,Tmgrad,Cgrad,Dim,comp,Pfit_pp,Pfit_pp_d1,eqR,timestep,timetol,timetolList,locsUp_P,Rfit_pp ), ...
%                 NH_KV_varList0, lb, ub, options);
            
           
            % lgd=legend('Envelope','Peaks of R','Exact','IMR Solver2 soln');
            % 
            % title(['G=',num2str(G),'; \mu=',num2str(mu),'; \lambda_\nu=',num2str(lambda_nu), ...
            %      '; \alpha=',num2str(alpha)],'FontWeight','Normal');
            % set(gca,'fontsize',18)
            % 
            % % fig_name = ['fig_exp',num2str(expt),'_Fung_mu',num2str(k),'_a',num2str(j),'_ln',num2str(i),'_pkloc.fig'];
            % fig_name = ['fig_simnum',num2str(simNo),'_NHKV_G',num2str(i),'_mu',num2str(j),'_a',num2str(k),'_ln',num2str(l),'_pkloc.fig'];
            % savefig(fig_name);

            % ==========================
            %if PlotDefineTimetolListOrNot==1, title(''); pause; end
            fprintf('================================================ \n')
       fprintf('Start defining time segments. After clicking all data points, press "Enter". \n');
            if DefineTimetolListOrNot==1
                [row1, col1] = ginput;
%                 if PickNo==7
%                     tempdata = 1e6*[tspanLocsMaxP(1)-row1(1);row1(2)-tspanLocsMaxP(1);tspanLocsMaxP(2)-row1(3);
%                         row1(4)-tspanLocsMaxP(2);tspanLocsMaxP(3)-row1(5);row1(6)-tspanLocsMaxP(3);tspanLocsMaxP(4)-row1(7);
%                         row1(8)-tspanLocsMaxP(4);tspanLocsMaxP(5)-row1(9);row1(10)-tspanLocsMaxP(5);tspanLocsMaxP(6)-row1(11);]'
%                     pause;
%                 elseif PickNo==4
%                     tempdata = 1e6*[tspanLocsMaxP(1)-row1(1);row1(2)-tspanLocsMaxP(1);tspanLocsMaxP(2)-row1(3);
%                         row1(4)-tspanLocsMaxP(2);tspanLocsMaxP(3)-row1(5);]'
%                     pause;
%                 elseif PickNo==11
%                     tempdata = 1e6*[tspanLocsMaxP(1)-row1(1);row1(2)-tspanLocsMaxP(1);tspanLocsMaxP(2)-row1(3);
%                         row1(4)-tspanLocsMaxP(2);tspanLocsMaxP(3)-row1(5);row1(6)-tspanLocsMaxP(3);tspanLocsMaxP(4)-row1(7);
%                         row1(8)-tspanLocsMaxP(4);tspanLocsMaxP(5)-row1(9);row1(10)-tspanLocsMaxP(5);tspanLocsMaxP(6)-row1(11);
%                         row1(12)-tspanLocsMaxP(6);tspanLocsMaxP(7)-row1(13);row1(14)-tspanLocsMaxP(7);tspanLocsMaxP(8)-row1(15);
%                         row1(16)-tspanLocsMaxP(8);tspanLocsMaxP(9)-row1(17);row1(18)-tspanLocsMaxP(9);tspanLocsMaxP(10)-row1(19);]'
%                     pause;
%                 else
                     for tempkk = 1:length(row1)
                         if mod(tempkk,2) == 1
                            tempdata(tempkk) = tspanLocsMaxP(floor(tempkk/2)+1)-row1(tempkk);
                         else
                            tempdata(tempkk) = row1(tempkk)-tspanLocsMaxP(floor(tempkk/2));
                         end
                     end
                     fprintf('Finish defining time segments. \n');
                  
%                 end
                
                timetolList(expt,1:PickNo*2-1) = reshape(tempdata(1:PickNo*2-1),1,PickNo*2-1);
                
                %timetolList = 1e-6*ones(max(expts),PickNo*2-1);
                
            end
        
             % if DefineTimetolListOrNot==0
            %     if (expt==1) && (GaussNo == 5) && (SampleNo == 4)
            %         timetolList(1,:) = 1e-6*[0.2887    1.0554    0.5418    0.4663 -0.0127    1.0207   -0.1131    1.6252 0.8191    1.1971    0.7202];
            %     elseif (expt==1) && (GaussNo == 5) && (SampleNo == 3)
            %         timetolList(1,:) = 1e-6*[0.9139    0.9006    0.7388    0.7733 0.6703    0.5394    0.7266    0.9367 0.9836    0.8310    1.1405];
            %     elseif (expt==1) && (GaussNo == 5) && (SampleNo == 2)
            %         timetolList(1,:) = 1e-6*[2.1440    1.1826    0.8973    1.0684 0.9766    1.8964    1.1983    1.2211  1.2353    1.6377    1.6538];
            %     elseif (expt==1) && (GaussNo == 5) && (SampleNo == 1)
            %         timetolList(1,:) = 1e-6*[1.7710    2.3117    2.0756    2.0071 2.0132    2.0694    1.7035    2.5304 2.0354    1.7449    1.8866];
            %     elseif (expt==1) && (GaussNo == 1) && (SampleNo == 1)
            %         timetolList(1,:) = 1e-6*[2.2855    2.3228    1.9499 2.0824    1.2597    2.1965 1.4776    1.4026    1.1113 1.7689    1.0668];
            %     elseif (expt==6) && (GaussNo == 5) && (SampleNo == 4)
            %         timetolList(6,:) = 1e-6*[0.6720    0.5376    0.8611    0.9534 2.112];
            %     elseif (expt==6) && (GaussNo == 5) && (SampleNo == 3)
            %         timetolList(6,:) = 1e-6*[0.9188    0.7445    1.0148    1.1022 2.0774];
            %     elseif (expt==6) && (GaussNo == 5) && (SampleNo == 2)
            %         timetolList(6,:) = 1e-6*[1.0165    1.2516    1.3216    1.0978 2.1584];
            %     elseif (expt==6) && (GaussNo == 5) && (SampleNo == 1)
            %         timetolList(6,:) = 1e-6*[1.7887    2.2940    2.7625    1.6226 3.4296];
            %     elseif (expt==11) && (GaussNo == 5) && (SampleNo == 4)
            %         timetolList(11,:) = 1e-6*[0.5639    0.6458    0.5664    0.7945 0.8611    0.9534    0.6559    0.7050 ...
            %             0.5755    0.4830    0.3981    0.8116  0.5700    0.6397    0.4039    0.6546 0.7067    0.8054    0.8169];
            %     elseif (expt==11) && (GaussNo == 5) && (SampleNo == 3)
            %         timetolList(11,:) = 1e-6*[0.9831    1.1339    0.7122    0.9511 0.9524    0.7109    0.8169    0.8464 ...
            %             0.7352    0.9281    0.8445    0.9700 0.6592    1.1553    0.8597    0.8036 1.2368    0.8802    0.6508];
            %     elseif (expt==11) && (GaussNo == 5) && (SampleNo == 2)
            %         timetolList(11,:) = 1e-6*[1.1249    1.1432    1.5122    1.2095 0.9104    2.2650    1.1220    0.9950 ...
            %             1.1386    2.1880    1.0080    1.1090 0.9460    1.4733    1.3148    1.2558   1.6879    0.8827    2.5182];
            %     elseif (expt==11) && (GaussNo == 5) && (SampleNo == 1 )
            %         timetolList(11,:) = 1e-6*[2.8703    1.9684    1.9896    1.9419 2.2022    2.0317    3.4673    1.3714 ...
            %             2.4286    1.5028    3.0525    0.8789 1.0677    1.2004    2.9092    0.8710 1.8112    0.9106    1.0658  ];
            %         %timetolList(1,:) = 1e-6*[1.7247    1.7315    1.6708    1.7854  ...
            %         % 1.8305    1.6257    1.4642    1.9921   0.6912    2.1890    0.4806 ];
            %     end
            % end
            fprintf('================================================\n');
            fprintf('Start tuning material properties using Nelder-Mead method.\n');
            fprintf('Please wait.\n');
            fprintf('================================================\n');
            
            if PlotDefineTimetolListOrNot == 1 && PlotExact == 1
              
                figure; % plot(t,Rnew(1,:)*1e-6,'rs'); hold on;
                plot(tspanLocsMaxR,pksMaxR,'kd') ;hold on;
                plot(tspanLocsMaxR,pksMaxR,'rd') ;hold on;
                plot(tv_exp,d_exp,'r'); hold off;
                grid off; grid; set(gca,'fontsize',18); grid minor;
                % title(['G=',num2str(10^G_ooms),' \mu=',num2str(10^mu_ooms),', \alpha=',num2str(alpha)],'FontWeight','Normal');
                title('Experiment data','FontWeight','Normal'); set(gca,'fontsize',18)
                axis([0,tspan,0,1.12*R0]);
                for tempi = 1:PickNo
                    hold on; plot(tspanLocsMaxP(tempi)*ones(1,101), linspace(0,2.5e-4,101),'k--');
                end
  
                PlotIMR2FitResults(model,matPropVar_best,tspanLocsMaxP(1:PickNo),pksMaxR(1:PickNo),tspanLocsMaxR(1:PickNo+1),NT/10,NTM,...
                    Pext_type,Pext_Amp_Freq,disptime,Tgrad,Tmgrad,Cgrad,Dim,comp,Pfit_pp,Pfit_pp_d1,eqR,timestep,...
                    timetolList(expt,:),Rfit_pp);
            else

            [matPropVar_best,fval,exitflag,output] = fminsearchbnd( @(matPropVar) ...
                funIMR2LSQErr(model,matPropVar,tspanLocsMaxP(1:PickNo),pksMaxR(1:PickNo),tspanLocsMaxR(1:PickNo+1),NT/10,NTM,...
                Pext_type,Pext_Amp_Freq,disptime,Tgrad,Tmgrad,Cgrad,Dim,comp,Pfit_pp,Pfit_pp_d1,eqR,timestep,timetolList(expt,:),Rfit_pp), ...
                matPropVar0, lb,ub,options);
            
            end
              
            disp(['expt = ',num2str(expt),'; G_best = ',num2str(10^matPropVar_best(1)),'; mu_best = ', num2str(10^matPropVar_best(2))]);
            
             
            toc
            
            matPropVarList{expt} = matPropVar_best;
            
            
            %%
            % CheckConstantModulusModel = A2*CaTilde + B2*ReTilde + E;
            % figure; plot(t2_num,CheckConstantModulusModel);
            % file_name = ['simnum',num2str(simNo),'_NHKV_IMRSolver2_G',num2str(G),'_mu',num2str(k),'_a',num2str(j),'_ln',num2str(i)];


            %%
%             Rdiff = diff(R2_num);
%             inc_idx = find(Rdiff>0);
%             diffzeros_idx = find(Rdiff(1:(end-1)).*Rdiff(2:(end))<0)+1;
% 
%             localmax_idx = [1; diffzeros_idx(2:2:end)];
%             localmin_idx = diffzeros_idx(1:2:end);
%             Rratios = R2_num(localmax_idx)/R2_num(1);
%             tcs = t2_num(localmin_idx);
%             tcs_star = tcs*Uc/R0;
%             tmaxs_star = t2_num(localmax_idx)*Uc/R0;
              %tpeaksall_star = [0; reshape([tcs_star,tmaxs_star(2:end)],size(diffzeros_idx))];

%             soln_mx{i,j,k,l}.tcs_star = tcs_star;
%             soln_mx{i,j,k,l}.Rratios = Rratios;
%             soln_mx{i,j,k,l}.tmaxs_star = tmaxs_star;
%             soln_mx{i,j,k,l}.t2_num = t2_num;
%             soln_mx{i,j,k,l}.R2_num = R2_num;
% 
%             % soln_mx{i,j,k}.U = U;
%             soln_mx{i,j,k,l}.P = P;
%             % soln_mx{i,j,k}.J = J;
%             soln_mx{i,j,k,l}.T = T;
%             soln_mx{i,j,k,l}.C = C;
%             soln_mx{i,j,k,l}.Cdel = Cdel;
%             soln_mx{i,j,k,l}.tdel = tdel;
%             soln_mx{i,j,k,l}.Tdel = Tdel;
% 
%             soln_mx{i,j,k,l}.LSQErr = LSQErr;
%             soln_mx{i,j,k,l}.LSQErrNo = LSQErrNo;

            % toc
%        end
        %If you want to save intermediate files during runs, uncomment this line
        % save([savename 'G_' num2str(G_ooms(i)) '.mat'],'soln_mx','G_ooms','mu_ooms','Rnew','t');
%        end
%    end
   % toc
%    end
    
    % cd(['/Users/yangjin/Documents/MATLAB/Franck/IMR-master/data/numsim/',num2str(simNo),'/']);
    % cd('/Users/yangjin/Documents/MATLAB/Franck/IMR-master/');
    %cd('E:\Jin\Franck\IMR-master\');
    
    % save([file_name '.mat'],'t2_num','R2_num','C','T','P');
    % save([savename '.mat'],'soln_mx','G_ooms','mu_ooms','Rnew','t');
    % cd(fp)
    
    %%
    
%     if strcmp(model,'fung')==1
%         G_ooms = alpha_ooms;
%     end
%     
%     LSQErrMatrix = zeros(length(G_ooms),length(mu_ooms)); LSQErrNoMatrix = LSQErrMatrix;
%     for tempi = 1:length(G_ooms)
%         for tempj = 1:length(mu_ooms)
%             LSQErrMatrix(tempi,tempj) = soln_mx{tempi,tempj,1,1}.LSQErr;
%             LSQErrNoMatrix(tempi,tempj) = soln_mx{tempi,tempj,1,1}.LSQErrNo;
%         end
%     end

    
    
    % 
    % LSQErrMatrix = zeros(length(G_ooms),1);
    % for tempi = 1:length(G_ooms)
    %     LSQErrMatrix(tempi) = soln_mx{tempi,1,5,1}.LSQErr;
    % end
    % figure; plot(G_ooms,LSQErrMatrix); set(gca,'yscale','log');
    % 
    % LSQErrMatrix = zeros(length(mu_ooms),1);
    % for tempi = 1:length(mu_ooms)
    %     LSQErrMatrix(tempi) = soln_mx{1,tempi,5,1}.LSQErr;
    % end
    % figure; plot(mu_ooms,LSQErrMatrix); set(gca,'yscale','log');

    %%
 
 
%     [mu_oomsGrid,G_oomsGrid] = meshgrid(mu_ooms,G_ooms);
%     figure; surf(G_oomsGrid,mu_oomsGrid,log10(LSQErrMatrix));  
%     axis([ (G_ooms(1)), (G_ooms(end)),(mu_ooms(1)), (mu_ooms(end)),-9, -7]); caxis([-9, -7]);
%     % axis auto; caxis auto;
%     colormap gray; set(gca,'fontsize',20); xlabel('log(G)');ylabel('log(\mu)');zlabel('LSQ error');
% 
% 
%     [mu_Peak,G_Peak,LSQErrMin] = findpeak((-log10(LSQErrMatrix)),1); 
%     G_oomsfit_pp = interp1(1:length(G_ooms),G_ooms,'pchip','pp');
%     G_Peak = ppval(G_oomsfit_pp,G_Peak); 10^G_Peak
%     mu_oomsfit_pp = interp1(1:length(mu_ooms),mu_ooms,'pchip','pp');
%     mu_Peak = ppval(mu_oomsfit_pp,mu_Peak); 10^mu_Peak
%     LSQErrMin
% 
%     hold on; p1=plot3(G_Peak,mu_Peak,-LSQErrMin+0.01,'ro','MarkerFaceColor','r');
%     LSQErrMinExact  = interp2(mu_oomsGrid,G_oomsGrid,-log10(LSQErrMatrix),log10(0.01), log10(2.97e3) );
%     % hold on; p2=plot3(log10(2.97e3),log10(0.01),-LSQErrMinExact+0.02,'bo','MarkerFaceColor','b');
%     % lgd = legend([p1,p2],'Numerical','Exact');
% 
%     axis tight; view([-1, -2.5, 3.5])
%     axis([ (G_ooms(1)), (G_ooms(end)),(mu_ooms(1)), (mu_ooms(end)),-LSQErrMin-0.5, -LSQErrMin+2]); caxis([-LSQErrMin-0.5, -LSQErrMin+2]);
%     % axis([3.4,  3.82,(mu_ooms(1)), (mu_ooms(end)),-LSQErrMin-0.5, -LSQErrMin+2]); caxis([-LSQErrMin-0.5, -LSQErrMin+2]);
%     
%        
%     % save(['LSQErrMatrix_stiff_expt',num2str(expt),'_PeakAll3_',num2str(PickNo),'.mat'],'LSQErrMatrix','LSQErrNoMatrix','G_ooms','mu_ooms');
%     % save(['LSQErrMatrix_numsim',num2str(simNo),'_PeakAll.mat'],'LSQErrMatrix');
% 
%     % ****** Generate G_accept and mu_accept LSQMatrix mask ******
% %      % LSQErrMatrix 
% %     [ G_accept, mu_accept ] = find( -log10(LSQErrMatrix)-LSQErrMin > - 0.5);
% %       
% %     LSQErr_acceptMask = zeros(length(G_ooms),length(mu_ooms));
% %     for tempi = 1:length(mu_accept)
% %         LSQErr_acceptMask(G_accept(tempi),mu_accept(tempi)) = 1;
% %     end
% %     figure; surf(G_oomsGrid,mu_oomsGrid,LSQErr_acceptMask );
% %     set(gca,'fontsize',20); xlabel('log(G)');ylabel('log(\mu)');
% %     axis tight;    view(2); colormap gray; %  view([-1,-2.5, 3.5]); 
% %     
% %     LSQErr_acceptMaskGlobal = LSQErr_acceptMaskGlobal+LSQErr_acceptMask.*log10(LSQErrMatrix);
% %     
% %     end
% % 
% %     figure; surf(G_oomsGrid,mu_oomsGrid,LSQErr_acceptMaskGlobal);
% %     set(gca,'fontsize',20); xlabel('log(G)');ylabel('log(\mu)');
% %     axis tight;    view(2); colormap gray; %  view([-1,-2.5, 3.5]); 

  
    %% 
%     save(['LSQErrMatrix_numsim',num2str(simNo),'_NH_expt',num2str(expt),'_PeakAll_',num2str(PickNo),'_100.mat'],'LSQErrMatrix','LSQErrNoMatrix','alpha_ooms','mu_ooms');

    %% LSQErrMatrix Hessian matrix analysis
    
%     if mu_Peak<-2.9999
%         mu_Peak=mu_Peak+0.1;
%     end
%     Glin = G_ooms(1):G_ooms_step:G_ooms(end); mulin=mu_ooms(1):mu_ooms_step:mu_ooms(end);
%     [muNodes,GNodes]=meshgrid(mulin,Glin);
%     [LSQErrMatrix2] = griddata(G_oomsGrid(:) ,mu_oomsGrid(:) ,log10(LSQErrMatrix(:)) ,GNodes,muNodes,'cubic');
%     [LSQErrMatrix2] = gridfit(G_oomsGrid(:) ,mu_oomsGrid(:) ,log10(LSQErrMatrix(:)),Glin,mulin,'smoothness',10); LSQErrMatrix2=LSQErrMatrix2';
%     figure; surf(GNodes ,muNodes ,LSQErrMatrix2); 
%     hold on; p1=plot3(G_Peak,mu_Peak,-LSQErrMin+0.0001,'ro','MarkerFaceColor','r');
%      LSQErrMinExact  = interp2(mu_oomsGrid,G_oomsGrid,-log10(LSQErrMatrix),log10(0.01), log10(2.97e3) );
%     axis tight; view([-1,-2.5, 3.5])
%     axis([ (G_ooms(1)), (G_ooms(end)),(mu_ooms(1)), (mu_ooms(end)),-LSQErrMin-0.01, -LSQErrMin+1]); caxis([-LSQErrMin-0.2, -LSQErrMin+1]);
%     
%     temp = ((GNodes-G_Peak).^2+(muNodes-mu_Peak).^2) ;
%     [temprow,tempcol]=find(temp==min(min( temp )));
%     tempDist=6;
%     DfAxis = [temprow-tempDist,temprow+tempDist,tempcol-tempDist,tempcol+tempDist];
%     [DLSQErrDx,DLSQErrDy,DDLSQErrDxx,DDLSQErrDyy,DDLSQErrDxy] = funGradient(DfAxis,LSQErrMatrix2);
%  
% 	Hessian11=median(reshape(DDLSQErrDxx(tempDist+1,tempDist+1),1,1))/(0.02^2);
%     Hessian22=median(reshape(DDLSQErrDyy(tempDist+1,tempDist+1),1,1))/(0.02^2);
%     Hessian12=median(reshape(DDLSQErrDxy(tempDist+1,tempDist+1),1,1))/(0.02^2);
%  
%     invHessian = inv([Hessian11,Hessian12;Hessian12,Hessian22]);
%     stdHessianG = sqrt(invHessian(1,1));
%     stdHessianmu = sqrt(invHessian(2,2));
%     
%     disp([num2str(10^G_Peak),'; ',num2str(10^mu_Peak),'; ',num2str(G_Peak),'; ',num2str(stdHessianG),'; ',num2str(mu_Peak),'; ',...
%         num2str(stdHessianmu),'; ']);

%%


end
end
% cd('/Users/yangjin/Documents/MATLAB/Franck/IMR-master');
%cd('E:\Jin\Franck\IMR-master\');

 

%%
% sz = size(soln_mx);
% sG1 = sz(1);
% try
%     sM = sz(2);
% catch
%     sM = 1;
% end
% try
%     sG = sz(3);
% catch
%     sG = 1;
% end
%%
%
% LSQ = cell(sG1,sM,sG);
% LSQminidx = [inf 0 0 0];
% try
%     for i=1:sG1
%        %figure(i)
%         for j=1:sM
%             for k=1:sG
%
%             %figure(10000*i+100*j+k)
%             %figure(10000*i+100+k)
%              %figure(505)
%
%                 [R2max idx] = max(soln_mx{i,j,k}.R2_num);
%
%                 hold on;
%
%                 t2_num = soln_mx{i,j,k}.t2_num-soln_mx{i,j,k}.t2_num(idx);
%                 R2_num = soln_mx{i,j,k}.R2_num;
%                 %plot(t2_num,R2_num,'Color',[(i-1)/sG1(1) (k-1)/sG(1) (j-1)/sM(1)]);
%
%                 [R0,t0] =  calcR0(Rnew(expt,:)*1E-6,t);
%                 t2_numexp=t-t0;
%                 R2exp=Rnew(expt,:)*1E-6;
%                 %plot(t2_numexp,R2exp, ' *');
%
%                 distRmat = bsxfun(@minus,R2_num,R2exp);
%                 disttmat = bsxfun(@minus,t2_num,t2_numexp);
%                 distmat = sqrt(disttmat.^2+distRmat.^2);
%                 weight = ones(1,101);
%                 %weight([1:10 18:end]) = 0; %for fitting the first peak
%                 %weight([1:10 19:26]) = 0; %for ignoring the second peak
%                 weight([1:9 70:end]) = 0;
%
%
%                 LSQ{i,j,k} = nansum(weight.*min(distmat,[],1));
%
%                 if LSQ{i,j,k}<LSQminidx(1)
%                     LSQminidx = [LSQ{i,j,k} i j k];
%                 end
%             end
% %             for n=k:-1:1
% %                 [maxRNHZ(n) idx(n)] = max(soln_mx{i,j,n}.R2_num);
% %             end
%
%         end
%
%     end
% catch
% end
%
% i=LSQminidx(2); j=LSQminidx(3); k=LSQminidx(4);
% log10([soln_mx{i,j,k}.G soln_mx{i,j,k}.mu soln_mx{i,j,k}.G1])
% [k j i]
% %%
% %k=7; j=8; i=16;
% %figure(9999); hold on;
% R2_num = soln_mx{i,j,k}.R2_num;
% [R2_nummax idx] = max(soln_mx{i,j,k}.R2_num);
% t2_num = soln_mx{i,j,k}.t2_num-soln_mx{i,j,k}.t2_num(idx);
% %plot(t2_num,R2_num);%,'Color','red');
% %plot(t2_num/max(R2_num)*10,R2_num/max(R2_num)./(1-((1-0.1301).*(exp(-t2_num/max(R2_num)*10/1.15))+0.1301)));%,'Color','red');
% %hold on; scatter(t2_numexp,R2_numexp,16,[0 0 0],'filled');
%
% [R2_nummax idx] = max(soln_mx{i,j,k}.R2_num);
%
% %save('dataprocoutputs.mat','t2_num' , 'R2_num','U','P','Z', 'T','C', 'Tm','tdel','Tdel','Cdel');
% end