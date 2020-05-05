%   Naringerin Metabolic pathway, Anthitetic controller, metabolic
%   extended biosensor and QdoR biosensor model.
%   ODE model with Direct controller (cI promoter)
%   Updated 04/30/2020 by Yadira Boada, Alejandro Vignoni
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dxdt] = model(t,x,p)

Ncell = 1;
Size = p.Size; 

for k = 1:Ncell
    m = Size*(k-1)+1;      

%% Genetic model
%x1 = NON   
c1 = 0*(p.ps*p.CN*p.ks./( p.dms+p.mu ));
dxdt(m,1)= c1*( p.alpha +(1-p.alpha)* x(m+4)^2./( p.kdlux*(p.kd2*p.CN./x(m+3))^2 + x(m+4)^2) )-...
          p.k_c./p.kdc*x(m)*x(m+1) + p.k_c*x(m+2) - (p.ds+p.mu)*x(m);
      
%% 
%x2 = cI
c2 =p.kdq*p.CNa;
c3 =p.kdk+ x(m+13); %Kae=x(m+13)
p.pa = 0.05*p.pa;
dxdt(m+1,1)= p.pa*p.CNa*p.ka./(p.dma+p.mu)*( p.alpha +...
             (1-p.alpha)*(c2^2*c3^2)./( c2^2*c3^2 + (p.kdk*x(m+6))^2))-...
             (p.da+p.mu)*x(m+1);
%%
%x3 = NON
dxdt(m+2,1) = 0*(p.k_c/p.kdc*x(m)*x(m+1) - p.k_c*x(m+2)-(p.dc+p.mu)*x(m+2));


%x4 = LuxR
dxdt(m+3,1)= p.pR*p.CN*p.kR./(p.dmR+p.mu) - (p.dR+p.mu)*x(m+3); 

%x5 = AHLint
dxdt(m+4,1) = p.D*p.Vcell/p.Vext*x(Size*Ncell+1)-p.D*x(m+4)- (p.dA+ p.mu)*x(m+4);
           
%x6 = CHS
c6 = p.CN*p.kh./( p.dmh+p.mu );
dxdt(m+5,1) = c6*p.beta*p.phc_cI+ c6*p.ph_cI*( p.alpha + ...
              (1-p.alpha)*p.kd_lamcI*(p.kdcI*p.CN)./( p.kd_lamcI*(p.kdcI*p.CN) + x(m+1)^2 ) )-...
              (p.dh+p.mu)*x(m+5);

%x7 = AHLext
%x8 = QdoR
dxdt(m+6,1) =  p.pq*p.CN*p.kq/(p.dmq+p.mu)-(p.dq+p.mu)*x(m+6); 

%% Metabolic model
%x9 = L-tyrosine
dxdt(m+7,1) =  p.KLt - p.catp*p.TAL*x(m+7)./(p.KmLt+x(m+7))-p.mu*x(m+7); 

%x10 = p-Coumeric acid
dxdt(m+8,1) = p.catp*p.TAL*x(m+7)./(p.KmLt+x(m+7)) -...
              p.catA*p.CL4*x(m+8)./(p.KmP+x(m+8))-p.mu*x(m+8); 

%x11 = p-CoA, CHS=x(m+5)
dxdt(m+9,1) = p.catA*p.CL4*x(m+8)./(p.KmP+x(m+8)) -...
              p.catNc*x(m+5)*(x(m+9)*p.Mal3./( 1/(p.Ka*p.Kb)+x(m+9)/p.Kb+p.Mal3/p.Ka+x(m+9)*p.Mal3))-...
              p.mu*x(m+9); 
                    
%x12 = Naringenin chalcone
dxdt(m+10,1) = p.catNc*x(m+5)*(x(m+9)*p.Mal3./( 1/(p.Ka*p.Kb)+x(m+9)/p.Kb+p.Mal3/p.Ka+x(m+9)*p.Mal3))-...
               p.catN*p.CHI*x(m+10)./(p.KmNc+x(m+10))- p.mu*x(m+10);
           
%x13 = Naringenin
dxdt(m+11,1) = p.catN*p.CHI*x(m+10)./(p.KmNc+x(m+10))-...
               p.catD*p.F3H*x(m+11)./(p.KmN+x(m+11)) - p.mu*x(m+11);

%x14 = Dihydrokaempferol
dxdt(m+12,1) = p.catD*p.F3H*x(m+11)./(p.KmN+x(m+11))-...
               p.catK*p.FLS*x(m+12)./(p.KmD+x(m+12)) - p.mu*x(m+12);
                      
%x15 = Kaempferol
dxdt(m+13,1) = p.catK*p.FLS*x(m+12)./(p.KmD+x(m+12)) - p.mu*x(m+13);

%x16 = Number of cells
dxdt(m+14,1) = p.mu*x(m+14)*(1-x(m+14)/p.cellmax);

end
%AHLext
    m=1;  %x(m+4)=AHL of Cell1
    dxdt((Size*Ncell+1),1) = -p.D*x(m+14)*p.Vcell/p.Vext*x(Size*Ncell+1) +...
                              p.D*x(m+14)*sum(x(m+4:Size:Size*Ncell)) - p.dAe*x(Size*Ncell+1);

end