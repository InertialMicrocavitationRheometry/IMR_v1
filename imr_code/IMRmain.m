%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% RUN IMR CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for expt = expts                    % looping through the number of experiments
    mkdir([fp num2str(expt)]);      % making folders for each of the experiments
    cd([fp num2str(expt)]);         % going into each of the folders
    soln_mx = cell(length(G_ooms),length(mu_ooms),length(G1_ooms));  % initializing the solution matrix
    %%%%%%%%%%%%%
    for k = 1:length(G1_ooms)
        for j= 1:length(mu_ooms)
            parfor i=1:length(G_ooms)
%             for i=1:length(G_ooms) 
                
             
                % CREATING SOLUTION STRUCTURE
                soln_mx{i,j,k} = struct('G',10^G_ooms(i),'mu',10^mu_ooms(j),'G1',10^G1_ooms(k),...
                    'tcs_star',[],'Rratios',[],'tmaxs_star',[],'t2',[],'R2',[],'U',[],'P',[],'J',[],'T',[],'C',[],...
                    'Cdel',[],'tdel',[],'Tdel',[]);            

                % OBTAINING PROPERTIES TO DETERMINE EQUILIBRIUM RADIUS
                G = soln_mx{i,j,k}.G;
                mu = soln_mx{i,j,k}.mu;
                G1 = soln_mx{i,j,k}.G1; 
                [R0,t0] =  calcR0(Rnew(expt,:)*1E-6,t); % Inital radius from a fit (R)
                eqR = median(Rnew(expt,62:end))*10^-6;
                R_eq = eqR/R0;
                % Calculating limits of the equilibrium radius by the
                % bisection method to get initial partial pressure of 
                % the non-condensible gas. Initial guess (purely empirical) 
                % and range (e.g. ~226 for 11kPa, 92 for 1.3kPa, 20.6 for water)
                P_guess = 226;
                P = [8, 400];
                R_eq_lims = [IMRCalc_Req(R0,Tgrad,Cgrad,P(1),G,G1,mu),...
                                IMRCalc_Req(R0,Tgrad,Cgrad,P(2),G,G1,mu)];
                R_eqf  = IMRCalc_Req(R0,Tgrad,Cgrad,P_guess,G,G1,mu);
                
                % WHILE LOOPING TO DETERMINE EQUILIBRIUM RADIUS
                error = 1;
                while abs(error) > 0.00001
                    if R_eq > R_eqf
                        P(1) = P_guess;
                        R_eq_lims(1) = IMRCalc_Req(R0,Tgrad,Cgrad,P_guess,G,G1,mu);
                    else
                        P(2) = P_guess;
                        R_eq_lims(2) = IMRCalc_Req(R0,Tgrad,Cgrad,P_guess,G,G1,mu);
                    end
                    P_guess = mean(P);
                    R_eqf  = IMRCalc_Req(R0,Tgrad,Cgrad,P_guess,G,G1,mu);
                    error = abs(R_eq-R_eqf);
                end

                % WHILE LOOPING TO DETERMINE EQUILIBRIUM RADIUS
                if strcmp(Pext_type,'IC')
                    Pext_Amp_Freq = [P_guess; 0];
                    %Pext_Amp_Freq = [226; 0]; % Tune first number to match equi radii
                elseif strcmp(Pext_type,'ga')
                    Pext_Amp_Freq = [P_guess; dt; tw;];
                    Rmax = R0;
                    R0 = eqR;
                end
      
                tic
                %Variables:
                % t2 = simulation time
                % R2 = simulation radius
                % R2dot = velocity of bubble wall
                % P = bubble pressure
                % S = stress integral
                % T = temperature inside bubble
                % C = relative concentration of vapor
                % Tm = temperature of wall
                [t2, R2, R2dot, P, S, T, C, Tm, tdel, Tdel, Cdel] = IMRsolver...
                    (model, G, G1, mu, tspan, R0, NT, NTM, ...
                    Pext_type, Pext_Amp_Freq , disptime, Tgrad, Tmgrad, Cgrad, Dim, comp);
                toc
                
                % FINDING THE INDICES OF MAXIMUM AND MINIMUM BUBBLE SIZE
                Rdiff = diff(R2);
                inc_idx = find(Rdiff>0);
                diffzeros_idx = find(Rdiff(1:(end-1)).*Rdiff(2:(end))<0)+1;
                localmax_idx = [1; diffzeros_idx(2:2:end)];
                localmin_idx = diffzeros_idx(1:2:end);
                % CALCULATING THE TIMES OF MAXIMUM AND MINIMUM BUBBLE SIZE
                Rratios = R2(localmax_idx)/R2(1);
                tcs = t2(localmin_idx);
                tcs_star = tcs*Uc/R0;
                tmaxs_star = t2(localmax_idx)*Uc/R0;
                %tpeaksall_star = [0; reshape([tcs_star,tmaxs_star(2:end)],size(diffzeros_idx))];
                
                % ASSIGNING THE VALUES TO THE SOLUTION MATRIX
                soln_mx{i,j,k}.tcs_star = tcs_star;
                soln_mx{i,j,k}.Rratios = Rratios;
                soln_mx{i,j,k}.tmaxs_star = tmaxs_star;
                soln_mx{i,j,k}.t2 = t2;
                soln_mx{i,j,k}.R2 = R2;
                %soln_mx{i,j,k}.U = U;
                %soln_mx{i,j,k}.P = P;
                %soln_mx{i,j,k}.J = J;
                %soln_mx{i,j,k}.T = T;
                %soln_mx{i,j,k}.C = C;
                soln_mx{i,j,k}.tdel = tdel;
                soln_mx{i,j,k}.Tdel = Tdel;
                soln_mx{i,j,k}.Cdel = Cdel;
                
                [i j k];
            end
            %If you want to save intermediate files during runs, uncomment this line
            %save([savename 'G_' num2str(G_ooms(i)) '.mat'],'soln_mx','G_ooms','mu_ooms','Rnew','t');
        end
    end
    % SAVING DATA INTO THE SAVE FILE 
    save([savename '.mat'],'soln_mx','G_ooms','mu_ooms','Rnew','t');
    cd(fp)      % needed for the next experiment
end