function   [RING_opt, pathlength,S,O, I] = scan_pathlength(RING)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script by Lina HOUMMI - last edited - Nov. 7th 2020                                                                     %
% INPUT : -RING, lattice structure                                                                                        %
%OUTPUT : - RINGopt, lattice of minimised pathlength                                                                      %
%         - pathlength, vector gathering the distance of all generated rings with the linear formula for the pathlength   %
%         - S, matrix of all generated sextupole sets, to restore all generated lattices                                  %
%         - O, matrix of all generated octupoles sets, to restore all generated lattices                                  %
%         - I, index of the minimised ring, for easier location within S and O                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RING1 = RING;% for linear optics calculation
RING2 = RING; % for derivation of the chromaticity matrix

n_oct = 0;% number of octupoles in the lattice
n_chrom = 0 ; % number of sextupoles in the lattice
cr_dim = 2;% only the linear chromaticities are fixed
index = []; % position of the sextupoles in the lattice
index_oct = [];% position of the octupoles in the lattice

%RING2 : all sextupoles and octupoles off
for i=1:length(RING2)
    if (strcmp(RING2{i}.Class,'Sextupole')==1),
        RING2{i}.PolynomB(1,3)=0;
        n_chrom = n_chrom + 1;
        index = [index, i]; 
    elseif (strcmp(RING2{i}.Class,'Marker')==1),
        n_oct = n_oct+1;
        index_oct = [index_oct, i];
    end
end
n_chrom=n_chrom/2; %for respect of the symmetry

A=rand(2,n_chrom); %Initialisaiton of the chromaticity matrix
n_omega = n_chrom - cr_dim;

%Scan over the sextupoles to extract their contribution to the
%chromaticity.
fprintf('Computing the chromaticity matrix A')
for k=1:n_chrom
    [lindata,tune,chrom]=atlinopt(RING2,0,1:length(RING2)+1); 
    nat_x = chrom(1);
    nat_y = chrom(2);
    if (strcmp(RING2{index(k)}.Class,'Sextupole')==1),
        RING2{index(k)}.PolynomB(1,3)=RING1{index(k)}.PolynomB(1,3);
        [lindata,tune,chrom]=atlinopt(RING2,0,1:length(RING2)+1);
        c_mags(k).k2 = RING1{index(k)}.PolynomB(1,3);
        A(1,k) = (chrom(1)-nat_x/2)/RING2{index(k)}.PolynomB(1,3);
        A(2,k) = (chrom(2)-nat_y/2)/RING2{index(k)}.PolynomB(1,3);
        RING2{index(k)}.PolynomB(1,3)=0;
    end
end

fprintf('Done')
Ap=pinv(A); %Moore-Penrose pseudoinverse
Ident=eye(n_chrom);
B=Ident-mtimes(Ap,A);
% Making matrix Q1
Q1=B(:,1:n_omega);
% Making Q1 triangular  
for i=1:(n_chrom-2)
    Q1(i,i+1:n_omega)=0;
end
%Making Q1t matrix
Q1t=transpose(Q1);% link between the sextupole space and the variable space of constant chromaticities. cf MOGA-Bmad

% Nonlinear magnets parameters
marging = 0.1; % relative variation imposed on the sextupole strengths: default = 10%
marging_oct = 500; % relative variation imposed on the octupoles strengths: default = \pm 500 m^-3
fprintf('Defining the omega bounds...')

for k=1:n_chrom
    c_mags(k).uir = c_mags(k).k2*(1+marging);
    c_mags(k).lir = c_mags(k).k2*(1-marging);
end

%Constraints on omega
for j=1:n_omega
    omega_bound_lir = 0;
    omega_bound_uir = 0;
    for k=1:n_chrom
        if (Q1t(j,k)<0)  
             omega_bound_lir = omega_bound_lir + Q1t(j,k)*c_mags(k).uir;
             omega_bound_uir = omega_bound_uir + Q1t(j,k)*c_mags(k).lir;
        else
            omega_bound_lir = omega_bound_lir + Q1t(j,k)*c_mags(k).lir;
            omega_bound_uir = omega_bound_uir + Q1t(j,k)*c_mags(k).uir; 
        end
    end
end
fprintf('Done')

N=100;% number of generated rings
[lindata,tune,chrom]=atlinopt(RING,0,1:length(RING)+1);
beta = cat(1,lindata.beta);
alpha = cat(1,lindata.alpha);
s = cat(1,lindata.SPos);

%Parameters for pathlength tracking and plots
NJx =101;
Jxmax = 5E-7;
Jx = linspace(0,Jxmax,NJx);
xamp = sqrt(2*Jx*beta(1,1))*10^3; %amplitude
xampred = sqrt(2*Jx)*10^3; % normalised amplitude
xmax = sqrt(2*Jxmax)*10^3; %maximum amplitude

fprintf('Start scanning...')
% [ pathlength, S] = pathfit(N,RING,omega_bound_lir,omega_bound_uir,n_chrom,Q1,index,Jxmax) %scan with the sextupoles only
 
 %Scan with sextupoles and octupoles
[pathlength, S, O] = pathfitoct(N,RING,omega_bound_lir,omega_bound_uir,n_chrom,n_oct, Q1,index, index_oct, marging_oct, Jxmax) 


%Selection of the ring with the pathlength the closest to the linear formula
[min_norme,I]=min(pathlength-1)
%Creation of the optimised ring
RING_opt = RING;
        Nt = (n_chrom);
        for i=1:Nt, %respect de la symetrie
            RING_opt{index(i)}.PolynomB(1,3)= S(i,I)+ RING{index(i)}.PolynomB(1,3);
            RING_opt{index(n_chrom*2+1-i)}.PolynomB(1,3)=S(i,I)+ RING{index(i)}.PolynomB(1,3);
        end
        for i=1:n_oct, 
            RING_opt{index_oct(i)}.PolynomB(1,4)= O(i,I)+ RING{index_oct(i)}.PolynomB(1,4);
        end
        RING_opt= chromaticity(RING_opt); %Setting the chromaticity at (1,1)

fprintf('Scan over, thank you for waiting')
fprintf('Best solution :')

%% Comparaisons
figure(41)
set(gcf,'color','w')
set(gca,'fontsize',20');
hold on
[lindata,tune,chrom]=atlinopt(RING_opt,0,1:length(RING_opt)+1);
beta = cat(1,lindata.beta);
alpha = cat(1,lindata.alpha);
chromX = chrom(1)
chromY = chrom(2)
DC_lin = -2*pi*chromX*Jx; %linear formula of the pathlenth with the chromaticity
plot(xampred, DC_lin,'r','Linewidth',2)
xlabel('Normalised Amplitude (m/sqrt(beta))')
ylabel('Averaged DC (m)')

%% Parameters
[lindata,tune,chrom]=atlinopt(RING,0,1:length(RING)+1);
beta = cat(1,lindata.beta);
alpha = cat(1,lindata.alpha);
variabsxt=atVariableBuilder(RING,{'SXD2E','SXF1E'},{{'PolynomB',{1,3}}});
ConstrChrom=[...
    atlinconstraint(1,{{'chromaticity',{1}}},chromX,chromX,1)...
    atlinconstraint(1,{{'chromaticity',{2}}},chromY,chromY,1)];
% 
tol=1e-15;
RING_corr=atmatch(RING,variabsxt,ConstrChrom,tol,500,4);

%Tracking of the averaged pathlength of the nominal ring, with corrected chromaticities for comparison with the optimised ring
Nturns = 100;
C=zeros(1,NJx);
C_opt=zeros(1,NJx);
Cmoy = zeros(1,NJx);
w =ones(1,Nturns);
Nphi = 101;
phase = linspace(-pi,pi,Nphi);
% phase = 0;
    for j =1:Nphi,
        for i=1:NJx,
%%         Tracking
         phi = phase(j);
         X0 = [sqrt(2*Jx(i)*beta(1,1))*cos(phi);-sqrt(2*Jx(i)/beta(1,1))*(sin(phi)+alpha(1,1)*cos(phi));0;0;0;0];
         ROUT = ringpass(RING_corr,X0,Nturns);
         ROUT_opt = ringpass(RING_opt,X0,Nturns);
         z=ROUT(6,Nturns);
         zopt=ROUT_opt(6,Nturns);%DeltaC directement
         C(i) = C(i) + z/Nturns;
         C_opt(i) = C_opt(i) + zopt/Nturns;
        end
    end
     C = C/Nphi;
     C_opt = C_opt/Nphi;
figure(41)
plot(xampred,C,'g','Linewidth',3)
hold on
plot(xampred,C_opt,'b.','Linewidth',3)

%Function driving the scan with sextupoles
function [pathlength, S] = pathfit(N,RING,omega_bound_lir,omega_bound_uir,n_chrom,Q1,index,Jxmax)
     pathlength=zeros(1,N); %initialisation
     a = omega_bound_lir;
     b = omega_bound_uir;
     S=[]; %matrix gathering all sextupole strengths
     for j=1:N,
        str = ['\n Etape j = ', num2str(j)];
        fprintf(str)
        RING3=RING;
        r = (b-a).*rand(n_chrom-2,1) + a;% random generation of omega
		s = Q1*r;% from omega to sextupoles
		S=[S,s];%translated strengths of the sextupoles
        Nt = (n_chrom);
 
        for i=1:Nt, %symmetry respected for the non-interleaved principle // applied to a hybrid lattice
            RING3{index(i)}.PolynomB(1,3)= s(i)+ RING{index(i)}.PolynomB(1,3);
            RING3{index(n_chrom*2+1-i)}.PolynomB(1,3)= s(i)+ RING{index(i)}.PolynomB(1,3);
        end
        RING3_corr = chromaticity(RING3); % correction of the chromaticity
               
              [lindata,tune,chrom]=atlinopt(RING3_corr,0,1:length(RING3)+1);
              beta = cat(1,lindata.beta);
              alpha = cat(1,lindata.alpha);
              chromX3 = chrom(1)
              chromY3 = chrom(2)
			  
			  %Comparison of the pathlength of all generated rings, with the nominal lattice and the linear formula. Updated during the scan
              figure(10)
              hold on
              da(RING3_corr,100) % comparison of the on-momentum dynamic apertures
              NJx =101;
              Jx = linspace(0,Jxmax,NJx);        
              Nturns = 100;
              Nphi = 101;
              C=zeros(1,NJx);
              DC_lin = -2*pi*chromX3*Jx;
              xamp = sqrt(2*Jx*beta(1,1))*10^3;
              xampred = sqrt(2*Jx)*10^3;
              xmax = sqrt(2*Jxmax)*10^3;
              phase = linspace(-pi,pi,Nphi);
                             for k =1:Nphi,
                                           for i=1:NJx,
                                           phi = phase(k);
                                           X0 = [sqrt(2*Jx(i)*beta(1,1))*cos(phi);-sqrt(2*Jx(i)/beta(1,1))*(sin(phi)+alpha(1,1)*cos(phi));0;0;0;0];
                                           ROUT = ringpass(RING3_corr,X0,Nturns);
                                           z=ROUT(6,Nturns);
                                           C(i) = C(i) + z/Nturns; %averaged over the turns
                                           end
                             end
         C = C/Nphi; %averaged over the input phases
         pathlength(j) = abs(C(NJx)/DC_lin(NJx)); %distance from the linear formula
         figure(42)
         hold on
         plot(xampred,C,'Linewidth',2)
     end    
     

%Function driving the scan with sextupoles ad octupoles
function [pathlength, S, O] = pathfitoct(N,RING,omega_bound_lir,omega_bound_uir,n_chrom,n_oct, Q1,index, index_oct, marging_oct, Jxmax)
     pathlength=zeros(1,N);
     a = omega_bound_lir;
     b = omega_bound_uir;
     
     S=[];
     O=[];
     for j=1:N,
        str = ['\n Etape j = ', num2str(j)];
        fprintf(str)
        RING3=RING;
        r = (b-a).*rand(n_chrom-2,1) + a;
        o = marging_oct.*rand(n_oct,1);
        s = Q1*r;
        S=[S,s];%random generation of a set of sextupole strengths
        O = [O,o]; %random generation of a set of octupoles strengths
        Nt = (n_chrom);
       
        for i=1:Nt, 
            RING3{index(i)}.PolynomB(1,3)= s(i)+ RING{index(i)}.PolynomB(1,3);
            RING3{index(n_chrom*2+1-i)}.PolynomB(1,3)= s(i)+ RING{index(i)}.PolynomB(1,3);
        end
        
        for i=1:n_oct, %no symmetry required for octupoles
              RING3{index_oct(i)}.PolynomB(1,4)= o(i)+ RING{index_oct(i)}.PolynomB(1,4);
        end
        RING3_corr = chromaticity(RING3);
               
              [lindata,tune,chrom]=atlinopt(RING3_corr,0,1:length(RING3)+1);
              beta = cat(1,lindata.beta);
              alpha = cat(1,lindata.alpha);
              chromX3 = chrom(1)
              chromY3 = chrom(2)

              figure(10)
              hold on
              NJx =101;
              Jx = linspace(0,Jxmax,NJx);        
              Nturns = 100;
              Nphi = 101;
              C=zeros(1,NJx);
              DC_lin = -2*pi*chromX3*Jx;
              xamp = sqrt(2*Jx*beta(1,1))*10^3;
              xampred = sqrt(2*Jx)*10^3;
              xmax = sqrt(2*Jxmax)*10^3;
              phase = linspace(-pi,pi,Nphi);
                             for k =1:Nphi,
                                           for i=1:NJx,
                                           phi = phase(k);
                                           X0 = [sqrt(2*Jx(i)*beta(1,1))*cos(phi);-sqrt(2*Jx(i)/beta(1,1))*(sin(phi)+alpha(1,1)*cos(phi));0;0;0;0];
                                           ROUT = ringpass(RING3_corr,X0,Nturns);
                                           z=ROUT(6,Nturns);
                                           C(i) = C(i) + z/Nturns;
                                           end
                             end
         C = C/Nphi;
         pathlength(j) = abs(C(NJx)/DC_lin(NJx));
         
         figure(42)
         hold on
         plot(xampred,C,'Linewidth',2)
     end    
     
     %Correction of the chromaticities, for comparison
function RING_corr = chromaticity(RING)
variabsxt=atVariableBuilder(RING,{'SXD2E','SXD1E','SXF1E','SXF2E','SXF3E'},{{'PolynomB',{1,3}}});
chromX = 1;
chromY = 1;
ConstrChrom=[...
    atlinconstraint(1,{{'chromaticity',{1}}},chromX,chromX,1)...
    atlinconstraint(1,{{'chromaticity',{2}}},chromY,chromY,1)];
tol=1e-15;
RING_corr=atmatch(RING,variabsxt,ConstrChrom,tol,500,4);
