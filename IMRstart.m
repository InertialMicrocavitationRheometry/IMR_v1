%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% INERTIAL MICROCAVITATION RHEOMETRY CODE %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors:
% Jon Estrada, jonathan_estrada@alumni.brown.edu
% Brown Solid Mechanics, PhD '17
% Carlos Barajas, carlobar@umich.edu
% U-M Mechanical Engineering BS '16

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% PRELIMINARIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 0. GO TO PREPROCESSING TO GENERATE RofTdata.mat FILE FROM 
% EXPERIMENTAL DATA AND/OR MAKE SURE TO HAVE YOUR RofTdata.mat FILE 
% AVAILABLE IN DIRECTORY SPECIFIED BELOW

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% RUN PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
warning('off','all');

%%%%%%%%%%%%%%%%%% FILE DIRECTORY SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% STEP 1. ADD GENERAL FOLDER PATH FOR MATLAB (EXAMPLE PROVIDED)
addpath(genpath('/mnt/hdd/imr/'));

% STEP 2. ADD THE PATH TO DATA LOCATION (EXAMPLE PROVIDED)
fp = '/mnt/hdd/imr/testdata/160420/11kPa_PA/';

% STEP 3. FILE LOADING AND SAVING. RofTdata.mat, contains vars Rnew and t
load([fp 'RofTdata.mat']);  % Rnew has size [num_expts num_video_frames]
savename = '170821_sweep';  % File name of saved data
allRmax = max(Rnew,[],2);   % Creating variable for looping

% STEP 4. CREATE A PARALLEL POOL (OPTIONAL)
parallel = 1;               % 0: SERIAL CODE, 1: PARALLEL CODE
if parallel == 1
    curCluster = parcluster('local');
    curCluster.NumWorkers = 4;
    saveProfile(curCluster);
    pool = parpool(4);
end

% STEP 5. SET THE INDEX CHOICE OF EXPERIMENTS (EXAMPLE PROVIDED)
%expts = [12 14:19]; %water
%expts = [2,3,4,5,8,10,14,15,16,18,20,23,24]; %collagen
expts = 1;

%%%%%%%%%%%%%%%%%%%%%%% NUMERICAL SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% STEP 6. SELECT THE MODELS (OPTIONS PROVIDED)
%model = 'linkv';  (linear Kelvin-Voigt)
%model = 'neoHook'; (neo-Hookean Kelvin-Voigt)
%model = 'sls'; (Standard Linear Solid)
%model = 'nhzen'; (neo-Hookean Standard Solid)
model = 'neoHook';
Pext_type = 'IC';     %'IC' for Flynn, 'ga' for gaussian bubble growth
Tgrad = 1;            % 1: Temperature transfer between bubble and material
Cgrad = 1;            % 1: Vapor-non-condensible gas diffusion
Tmgrad = 0;           % 0: Off means cold liquid assumption
disptime = 0;         % 1 = Displays time to complete simulation
Dim = 1;              % 1 = displays results in dimensional form
comp = 1;             % 0 uses Rayleigh-Plesset, 1 uses Keller-Miksis
% STEP 7. SET NUMERICAL RESOLUTION
NT = 500;             % Mesh points in bubble, resolution should be >=500
NTM = 10;             % Mesh points in the material, should be 10

% STEP 8. SET THE IN SIMULATION TIME
tspan = 1.3E-4; 
% tspan = 2.25E-4;

%%%%%%%%%%%%%%%%%%%%%%% MATERIAL SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% STEP 9. SET THE RANGE OF G and MU (most conveniently done as powers of 10)
% G_ooms = 1:0.2:5;   
G_ooms = 1:0.5:5;   
%G_ooms = 3.0:0.1:4.0;      %soft PA 
%G_ooms = 3.65:0.05:4.15    %stiff PA
%mu_ooms = -1.4:0.05:-0.9;
mu_ooms = -1.4:0.1:-0.9;
%mu_ooms = [-inf -4:0.25:-1.25 -1.05];
%mu_ooms = [-1.65:0.05:-0.9];
%mu_ooms = [-2.25:0.25:-0.5];
%mu_ooms = [-inf -1.65:0.05:-0.9];
G1_ooms = inf;
%Note, mu/G1 should go to 0 in the limit, so G1 for K-V models should be infinite

% STEP 8. SETTING MATERIAL CONSTANTS
P_inf = 101325;             % (Pa) Atmospheric Pressure
rho = 998.2;                % (kg/m^3) Material Density
Uc = sqrt(P_inf/rho);       % (m/s) Characteristic velocity
               
% STEP 9. CHECK WORKSPACE IS CORRECT AND DATA IS AVAILABLE IN PATH ABOVE
addpath imr_code
IMRmain

%FOLLOWING PART OF CODE NEEDS TO BE CLEANED UP
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
%                 [R2max idx] = max(soln_mx{i,j,k}.R2);
%
%                 hold on;
%
%                 t2 = soln_mx{i,j,k}.t2-soln_mx{i,j,k}.t2(idx);
%                 R2 = soln_mx{i,j,k}.R2;
%                 %plot(t2,R2,'Color',[(i-1)/sG1(1) (k-1)/sG(1) (j-1)/sM(1)]);
%
%                 [R0,t0] =  calcR0(Rnew(expt,:)*1E-6,t);
%                 t2exp=t-t0;
%                 R2exp=Rnew(expt,:)*1E-6;
%                 %plot(t2exp,R2exp, ' *');
%
%                 distRmat = bsxfun(@minus,R2,R2exp);
%                 disttmat = bsxfun(@minus,t2,t2exp);
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
% %                 [maxRNHZ(n) idx(n)] = max(soln_mx{i,j,n}.R2);
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
% R2 = soln_mx{i,j,k}.R2;
% [R2max idx] = max(soln_mx{i,j,k}.R2);
% t2 = soln_mx{i,j,k}.t2-soln_mx{i,j,k}.t2(idx);
% %plot(t2,R2);%,'Color','red');
% %plot(t2/max(R2)*10,R2/max(R2)./(1-((1-0.1301).*(exp(-t2/max(R2)*10/1.15))+0.1301)));%,'Color','red');
% %hold on; scatter(t2exp,R2exp,16,[0 0 0],'filled');
%
% [R2max idx] = max(soln_mx{i,j,k}.R2);
%
% %save('dataprocoutputs.mat','t2' , 'R2','U','P','Z', 'T','C', 'Tm','tdel','Tdel','Cdel');
% end