%   Naringerin Metabolic pathway, Anthitetic controller, metabolic
%   extended biosensor and QdoR biosensor model.
%   Main execution script
%   Updated 04/30/2020 by Yadira Boada, Alejandro Vignoni
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% General parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ODmax = 12;
p = parameters(ODmax);

%Closed loop gain.  
p.phc = 6.5096e-04; %RBS of the constitutive promoter
p.ph = 15.6230;     %RBS of the inducible promoter (Antithetic controller)
p.phc_cI = 6.5096e-04;%RBS of the constitutive promoter
p.ph_cI = 54.89;    %RBS of the inducible promoter (Direct controller)

% Open loop (uncomment if you want to run the simulation of the Open-loop
%p.phc =  6.5493;      %Open loop gains. Comparative plots
%p.ph = 0;
%p.phc_cI = 6.5493;
%p.ph_cI = 0;


%Enzymes (molecules). Max values from each enzyme range.
p.TAL = 20*1.6e5;  
p.CL4 = 15*4.32e5;  
p.CHI = 10*3.54e5;
p.F3H = 2.81;
p.FLS = 5.84;
MAL_PERCENT = 1;


%System size
NumberStates = 16; 

%Initial OD
ODinitial = 0.001;
Cellinitial = ODinitial*p.Vext*p.OD_to_cells; 

                     
%% 0 Null initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
tfin = 60*8;     %simulation time
step = 0.1;
tspan = 0:step:tfin-step;
options = odeset('AbsTol',1e-8,'RelTol',1e-6);      % for ode function 

p.Mal3 = 0;             %Input: 3 Malonyl-CoA
p.Mal30 = p.Mal3;     %Input: 3 Malonyl-CoA
Initial = [zeros(1,NumberStates-2) Cellinitial 0];  %ini conditions[species, cells, ahle]
p.Size = length(Initial)-1;
[t0,x0] = ode23t(@(t,x) model(t,x,p),tspan, Initial, options);
Initial(14) = 00;
Initial(2) = 1000;
[t0,y0] = ode23t(@(t,x) model_cI_repressor(t,x,p),tspan, Initial, options);

%% 1 Adding Malonyl  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Initial = x0(end,1:end); %Initial conditions
Initial_y0 = y0(end,1:end); %Initial conditions
tfin = 60*5;             %Tiempo de simulacion (min)
tspan = 0:step:tfin-step;

p.Mal3 = 1.17e3;     %Mean amount Malonyl-CoA=3.54e-5 (M) from table
                     %Maximum Malonyl-CoA=3.09e-3  - Minimum 4.05e-7 (M)from table
p.Mal31 = p.Mal3;     %
[t1,x1] = ode23t(@(t,x) model(t,x,p),tspan, Initial, options);
[t1,y1] = ode23t(@(t,x) model_cI_repressor(t,x,p),tspan, Initial_y0, options);

%% 2 Adding ahle. Closing the loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ahle_nM = 50;     % Induction in the lab [nM]
p.Mal32 = p.Mal3;     %Input: 3 Malonyl-CoA

nM = 1e-9;       %nM in Molarity
to_molecules = p.Vext*p.nA*nM;
ahle0 = ahle_nM * to_molecules;

Initial = [x1(end,1:end-1) ahle0]; %Initial conditions
Initial_y0 = [y1(end,1:end-1) ahle0]; %Initial conditions
tfin = 60*15; % Tiempo de simulacion (min)
tspan = 0:step:tfin-step;
[t2,x2] = ode23t(@(t,x) model(t,x,p),tspan, Initial, options);
[t2,y2] = ode23t(@(t,x) model_cI_repressor(t,x,p),tspan, Initial_y0, options);

%% 3 Continue closed loop  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.Mal33 = p.Mal3;     %Input: 3 Malonyl-CoA
Initial = x2(end,1:end); %Initial conditions
Initial_y0 = y2(end,1:end); %Initial conditions

tfin = 60*20; % Tiempo de simulacion (min)
tspan = 0:step:tfin-step;
[t3,x3] = ode23t(@(t,x) model(t,x,p),tspan, Initial, options);
[t3,y3] = ode23t(@(t,x) model_cI_repressor(t,x,p),tspan, Initial_y0, options);

%% 4 Malonyl Perturbation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.Mal3 = MAL_PERCENT*p.Mal3;     %Input: 3 Malonyl-CoA
p.Mal34 = p.Mal3;     %Input: 3 Malonyl-CoA
Initial = x3(end,1:end); %Initial conditions
Initial_y0 = y3(end,1:end); %Initial conditions
tfin = 60*36;
tspan = 0:step:tfin-step;
[t4,x4] = ode23t(@(t,x) model(t,x,p),tspan, Initial, options);
[t4,y4] = ode23t(@(t,x) model_cI_repressor(t,x,p),tspan, Initial_y0, options);

%% Species from the Antithetic Controller simulation
sigma = [x1(:,1); x2(:,1); x3(:,1); x4(:,1)]; 
asigma = [x1(:,2); x2(:,2); x3(:,2); x4(:,2)];  
sa_complex = [x1(:,3); x2(:,3); x3(:,3); x4(:,3)]; 
luxR = [x1(:,4); x2(:,4); x3(:,4); x4(:,4)]; 
ahl = [x1(:,5); x2(:,5); x3(:,5); x4(:,5)]; 
chs = [ x1(:,6); x2(:,6); x3(:,6); x4(:,6)];  
qdoR = [x1(:,7); x2(:,7); x3(:,7); x4(:,7)]; 
Ltyrosine = [x1(:,8); x2(:,8); x3(:,8); x4(:,8)]; 
pC_acid = [x1(:,9); x2(:,9); x3(:,9); x4(:,9)]; 
p_CoA = [x1(:,10); x2(:,10); x3(:,10); x4(:,10)]; 
nchalcone = [ x1(:,11); x2(:,11); x3(:,11); x4(:,11)]; 
naringenin = [x1(:,12); x2(:,12); x3(:,12); x4(:,12)]; 
dykaempferol = [x1(:,13); x2(:,13); x3(:,13); x4(:,13)]; 
kaempferol = [x1(:,14); x2(:,14); x3(:,14); x4(:,14)]; 
ahle = [x1(:,16); x2(:,16); x3(:,16); x4(:,16)];
od = [x1(:,15); x2(:,15); x3(:,15); x4(:,15)]./(p.OD_to_cells*p.Vext);
malonyl = [p.Mal31.*ones(length(t1),1); p.Mal32.*ones(length(t2),1);...
           p.Mal33.*ones(length(t3),1); p.Mal34.*ones(length(t4),1)];
% Species from the Direct Controller simulation
cI = [y1(:,2); y2(:,2); y3(:,2); y4(:,2)];  
LUXR = [y1(:,4); y2(:,4); y3(:,4); y4(:,4)]; 
AHL = [y1(:,5); y2(:,5); y3(:,5); y4(:,5)]; 
CHS = [ y1(:,6); y2(:,6); y3(:,6); y4(:,6)];  
QDOR = [y1(:,7); y2(:,7); y3(:,7); y4(:,7)]; 
LTYrosine = [y1(:,8); y2(:,8); y3(:,8); y4(:,8)]; 
PC_Acid = [y1(:,9); y2(:,9); y3(:,9); y4(:,9)]; 
P_COA = [y1(:,10); y2(:,10); y3(:,10); y4(:,10)]; 
NChalcone = [ y1(:,11); y2(:,11); y3(:,11); y4(:,11)]; 
NAringenin = [y1(:,12); y2(:,12); y3(:,12); y4(:,12)]; 
DYkaempferol = [y1(:,13); y2(:,13); y3(:,13); y4(:,13)]; 
KAempferol = [y1(:,14); y2(:,14); y3(:,14); y4(:,14)]; 
AHLe = [y1(:,16); y2(:,16); y3(:,16); y4(:,16)];
OD_cI = [y1(:,15); y2(:,15); y3(:,15); y4(:,15)]./(p.OD_to_cells*p.Vext);
time = [t1; t1(end)+t2; t1(end)+t2(end)+t3; t1(end)+t2(end)+t3(end)+t4];

%% Metabolites production mg/L  
%Save data
total_product = [];  
total_mal = [];

%metabolites
%molecular_weight=[L-ty,p-Co acid,p-CoA,N chal,Naringenin,Dikaem,Kaemp, Malonyl];   %g/mol
molecular_weight=[181.19, 164.047, 913.67, 272.25, 272.25, 288.25, 286.23, 853.6]';   %g/mol

production_gL = (([x4(end,8:14), malonyl(end)]'.*molecular_weight.*p.OD_to_cells.*ODmax)./p.nA).*1e3;
production_gL_cI = (([y4(end,8:14), malonyl(end)]'.*molecular_weight.*p.OD_to_cells.*ODmax)./p.nA).*1e3;

