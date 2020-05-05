%   Naringerin Metabolic pathway, Anthitetic controller, metabolic
%   extended biosensor and QdoR biosensor model.
%   Parameters: structure contains all rates and constants
%   Updated 04/30/2020 by Yadira Boada, Alejandro Vignoni
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [p] = parameters(ODmax)
%%%%%%%%%%%%%%%%%%%%%%%%  General parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    p.doubling = 82;                            % doubling time [min]
    p.mu = log(2)/p.doubling;                   % growth rate [1/min]
    p.nA = 6.023e23;                            % Avogadro's number: # particles/mol
    p.Vcell =  1.1e-15;                         % typical volume of E. coli (liters). Source: Bionumbers
    p.Vext = 4e-3;                              %culture medium volume [l] in a plate reader.  microfluidic device = 1e-9 liters
                                                % From Solvej Siedler, Novel biosensors based on flavonoid
    p.OD_to_cells = 8e11;                       % cells/liter for OD=1, Agilent, E. coli Cell Culture
    p.cellmax = ODmax*p.Vext*p.OD_to_cells;     % maximum number of cells
    p.nM_molecules = 1e15/(p.Vcell*p.nA);       % Conversion factor from number of particles to concentration (1 nMolar=nanomols/liter)

    %Copy number 
    p.CN = 10;   % plasmid number  pACYC184 (10 copies/cell)
    %Copy number antisigma
    p.CNa = 10;   % plasmid number  pACYC184 (10 copies/cell)
    %Copy number CHS
    p.CNh = 10;
    
    %Sigma
    p.dms = log(2)/3;            	% degradation rate mRNA [1/min]
    p.ds = 0.0003;                  % degradation rate  [1/min]. Rapid degradation: RR Burgess, doi: 10.1006/rwgn.2001.1192
    p.ps =  3; 			   			% translation rate  [1/min] [2.9801 - 5.9603] from our rates calculator
    p.ks = 2*0.99; 		            % transcription rate [1/min] [0.99338 - 9.9338] from our rates calculator
    
    p.kd20 = 1000;        % dissociation cte to promoter [molecules], Annunziata 2017
    
    %Anti-sigma
    p.alpha = 0.01;             % basal expresion pLux
    p.dma = log(2)/3;           % mRNA degradation rate [1/min]
    p.da = 0.0003;               % protein degradation rate [1/min]
    p.pa =  0.8*3.96;               % translation rate  [1/min] [3.9648 - 7.9295]
    p.ka = 1.5*1.32;              % transcription rate [1/min] [1.3216 - 13.2159] from our rates calculator
    
    %cI
    p.kdcI = 30;
    p.kd_lamcI = 1000;                                % dissociation cte to promoter [molecules]. Strong promoter
    p.phc_cI = 6.5096e-04;
    p.ph_cI = 0.001*1.54;
    
    %SigmaComplex
    p.kdc = 0.01;                  % dissociation constant [molecules] Annunziata 2017 an orthogonal multi-input
    p.k_c = 10*1.8e-3;                % [1/min] Annunziata 2017 an orthogonal multi-input
    p.kc = p.k_c/p.kdc;            % binding rate sigma to anti-sigma [min^-1 molecules^-1]
    p.dc = 0.001;                   % degradation rate [1/min] Annunziata 2017 an orthogonal multi-input
    
    %LuxR
   
    p.beta1 = 0.01;                % basal expresion p20
    p.dmR = log(2)/3;              % mRNA degradation rate  [1/min]
    p.dR = 0.02;                   % degradation rate [1/min]
    p.pR = 2.34;                  % translation rate  [1/min] [1.1749 - 2.3499] from our rates calculator
    p.kR = 2*0.39;                   % transcription rate [1/min] [0.39164 - 3.9164] from our rates calculator
    
    % Monomer LuxR.AHL
    p.kd1 = 100;                   % dissociation constant of R to A [nM], Urbanowski etal. 2004
    p.k_1 = 10;                   % unbinding rate LuxR to AHL [1/min]
    p.k1 = p.k_1/p.kd1;            % binding rate LuxR to AHL [1/min]
    p.dRA = log(2)/5;              % degradation rate of (LuxR.A) [1/min]. Buchler et al. 2004 Monomer half-life is just few minutes.
    
    % Dimer (R.A)2
    p.kd2 = 20;                   %dissociation cte (LuxR.A) to (LuxR.A) [nM], Buchler et al. 2003
    %Koren, R. & Hammes, G. G. (1976) Biochemistry 15, 1165�1171.
    %Northrup, S. H. & Erickson, H. P. (1992) Proc. Natl. Acad. Sci. USA 89, 3338�3342.
    p.k_2 = 1;                     % dissociation rate [1/min]
    p.k2 = p.k_2/p.kd2;               % binding rate LuxR to AHL [1/min]
    p.kdlux = 600;           % dissociation cte (LuxR.A)2 to promoter [nM], Bucler et al [1 1000]nM
    
    %Sigma dimer
    p.kds = 1000;
    
    %CHS (389 amino acids)
    p.beta = 1.5;  
    p.dmh = log(2)/3;              % degradation rate mRNA [1/min]
    p.dh = 0.0003;                 % degradation rate  [1/min]. 
    
    %p.ph =   1.557*2*1.54;        %Open loop gain. Only comparative plots
    p.phc =  2*1.54;               % CONSTITUTIVE translation rate   
    p.ph = 60*1.54;                %translation rate [1/min] [1.5424 - 3.0848] from our rates calculator  %RBS of the constitutive promoter
    p.kh =  48*7.5*0.51;           %   transcription rate [1/min] [0.51414 - 5.1414] from our rates calculator
        
    %QdoR (842 bp)
    p.dmq = log(2)/3;             % degradation rate mRNA [1/min]
    p.dq = 0.0003;                % degradation rate  [1/min]. 
    p.pq =  1.2*2.13; %b*dms;         % translation rate  [1/min] [2.1378 - 4.2755] from our rates calculator
    p.kq = 0.71; %200/b;          % transcription rate [1/min] [0.71259 - 7.1259] from our rates calculator
    p.kdq = 150;                     % kd = 1 nM dissociation constant complex (Q.Kae)2 from promoter PqdoI.  'Novel biosensors based on 
                                                   %flavonoid-responsive transcriptional regulators introduced into E. coli'
    
    %Dissociation constant Kaempferol from QdoR. kdq>kdk
    p.kdk = 75;              %kd^n = 5 uM.  Ref.: Novel biosensors based on 
                                                   %flavonoid-responsive transcriptional regulators introduced into E.
    
    %AHL and AHLe
    p.D = 2;        %kinetic rate of AHL external transport [1/min] across the cell membrane, calculated
    
    p.dA = 0.0004;    %[0.05 0.03 min^-1]Degradation from Bionumbers online
    p.dAe = 0.0000481;          % Horswill et al., 2007  %0.0164, Degradation rate for AHL. From Kaufmann etal. 2005. Similar to You etal. Nature, 2004
    %[0.05 0.03 min^-1]Degradation from Bionumbers online
    %0.000282 Degradation rate for external AHL. Fitted using half-life of 180 minutes, from Englmann etal. 2007
    % In Kauffmann & Sartorio, 2005 they use 0.0018
    
    %TAL,L-tyrosine
    p.KLt = 2e6;   %L-tyrosine production rate (this is 400mg/L with OD=10 )
    p.KmLt = 1.9e-5*1e9;  %[1.9E-5M, 1.6E-4]M 
    p.catp = 0.02*60;    % 0.02- 4.32 (1/sec)
    
    %4CL, p-Coumeric acid
    p.catA = 8.20e-3*60;           %[8.20e-3, 8.87e1] 1/sec
    p.KmP =  1.4e-5*1e9;            %[1.40E-5, 0.000432]M
       
    %CHS,p-CoA
    p.catNc =0.007*60;    %[0.007,0.042]	 s^(-1)
    p.Ka =  1/(1.0e-6*1e9);    %[1.00E-6, 3.50E-2]M
    p.Kb =  1/(1.0e-6*1e9);    %malonyl-CoA  https://www.brenda-enzymes.org/enzyme.php?ecno=2.3.1.74
    
    %CHI, Naringenin chalcone
    p.catN = 0.07*60;    %[0.07,23]	 s^(-1)
    p.KmNc = 2.8e-5*1e9; %[2.40E-06, 3.54E-04]M   https://www.brenda-enzymes.org/enzyme.php?ecno=5.5.1.6
    
   %F3H, Naringenin
    p.catD= 0.5*1*60;     %[1.00E+00	3.90E+00]s^(-1)
    p.KmN = 100000*5e-6*1e9;   %[5.00E-06, 2.18E-04]M https://www.brenda-enzymes.org/enzyme.php?ecno=1.14.11.11
                                                % this value is 2293 times less effective
   %FLS, Dihydrokaempferol
    p.catK= 0.1*1*60;       %6.6 (1/sec) https://www.brenda-enzymes.org/enzyme.php?ecno=1.14.20.6
    p.KmD = 10e-6*1e9;    %[1.00E-06,5.84E-05]M    https://www.brenda-enzymes.org/enzyme.php?ecno=1.14.20.6

