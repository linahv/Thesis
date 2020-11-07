function  [RINGopt, alpha1all, S] =alpha1opt(RING)
%MAIN FUNCTION
ring_fun = RING;
%Slicing of some elements of the ring, for increases precision, in the sommation of the integral to calculate alpha1
for k=1:length(RING),
    if (strcmp(RING{length(RING)-k+1}.Class,'Bend')==1),
        nslices = 10;
        elem=splitelem(RING{length(RING)-k+1},nslices);
        ring_fun=cat(1,ring_fun(1:length(RING)-k),elem,ring_fun(length(RING)-k+2:length(ring_fun)));
    elseif (strcmp(RING{length(RING)-k+1}.Class,'Sextupole')==1),
        nslices = 1;
        elem=splitelem(RING{length(RING)-k+1},nslices);
        ring_fun=cat(1,ring_fun(1:length(RING)-k),elem,ring_fun(length(RING)-k+2:length(ring_fun)));
    else
        nslices = 5;
        elem=splitelem(RING{length(RING)-k+1},nslices);
        ring_fun=cat(1,ring_fun(1:length(RING)-k),elem,ring_fun(length(RING)-k+2:length(ring_fun)));
    end
end

%Definition of different rings
RING1=ring_fun'; % for linear optics calculation
RING2 = ring_fun'; % for derivation of the chromaticity matrix

n_chrom = 0 ; % number of sextupoles in the lattice
cr_dim = 2;% only the linear chromaticities are fixed
index = []; % position of the sextupoles in the lattice

%RING2 : all sextupoles off
for i=1:length(RING2)
    if (strcmp(RING2{i}.Class,'Sextupole')==1),
        RING2{i}.PolynomB(1,3)=0;
        n_chrom = n_chrom + 1;
        index = [index, i];
    end
end

A=rand(2,n_chrom); %Initialisaiton of the chromaticity matrix
n_omega = n_chrom - cr_dim;

%Scan over the sextupoles to extract their contribution to the
%chromaticity.
j = waitbar(0,'Computing the chromaticity matrix A');
for k=1:n_chrom
    waitbar(k/n_chrom)
    [lindata,tune,chrom]=atlinopt(RING2,0,1:length(RING2)+1);
    nat_x = chrom(1);
    nat_y = chrom(2);
    if (strcmp(RING2{index(k)}.Class,'Sextupole')==1),
        RING2{index(k)}.PolynomB(1,3)=RING1{index(k)}.PolynomB(1,3);
        [lindata,tune,chrom]=atlinopt(RING2,0,1:length(RING2)+1);
        c_mags(k).k2 = RING1{index(k)}.PolynomB(1,3);
        A(1,k) = (chrom(1)-nat_x)/RING2{index(k)}.PolynomB(1,3);
        A(2,k) = (chrom(2)-nat_y)/RING2{index(k)}.PolynomB(1,3);
        RING2{index(k)}.PolynomB(1,3)=0;
    end
end
close(j)
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
Q1t=transpose(Q1); % link between the sextupole space and the variable space of constant chromaticities. cf MOGA-Bmad

% Parameters
marging = 0.1; % relative variation imposed on the sextupole strengths: default = 10%

fprintf('Defining the omega bounds...') % omega, defined in MOGA-Bmad, corresponds to the variables in the space of constant chromaticities
for k=1:n_chrom
    c_mags(k).uir = c_mags(k).k2*(1+marging);
    c_mags(k).lir = c_mags(k).k2*(1-marging);
end

%Constraints on omega
l = waitbar(0,'Constraints on omega under calculations');
for j=1:n_omega
    waitbar(j/n_omega)
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
close(l)
fprintf('Done')

N=1000; % number of generated rings
y=1:1:N; % for analysis
fprintf('Start scanning...')
 [alpha1all, S, d0, dd0, d1] = scan_disp1_list(RING1, N, omega_bound_lir, omega_bound_uir, n_chrom, index,Q1); % scan of disp1 and alpha1

% Printing the evolution of alpha1 for all generated rings
figure(10)
set(gcf,'color','w')
set(gca,'fontsize',16');
plot(y, alpha1all, 'b.');
xlabel('Iteration number')
ylabel('alpha1')
title('Variation of alpha1 during the scan of sextupoles')

%Selection of the ring of minimum alpha1
[min_norm,I] = min(alpha1all(1,:))

% Relative variation of the sextupole strengths for the selected ring
figure(11)
set(gcf,'color','w')
set(gca,'fontsize',16');
y=1:1:n_chrom;
var = zeros(1,n_chrom)
for i=1:n_chrom
    var(i) = (S(i,I))/RING1{index(i)}.PolynomB(1,3)*100;
end
figure(11)
plot(y, var, 'k.');
xlabel('Sextupole index')
ylabel('Relative variation of the optimised sextupoles')

fprintf('Scan over, thank you for waiting')

%Creation of the optimised ring
RINGopt = RING1;
for i=1:n_chrom,
        RINGopt{index(i)}.PolynomB(1,3)= S(i,I)+ RING1{index(i)}.PolynomB(1,3);
end

%Comparison of the first-order disperion of the nominal ring (red) and the optimised ring (blue)
figure(12)
title('First-order dispersion of the optimised ring')
N1 = 1000;
[d0, dd0]= disp0(RINGopt, N1);
[d1, dd1]=disp1(RINGopt,d0,dd0,N1);
[lindata,tune,chrom]=atlinopt(RINGopt,0,1:length(RING1)+1);
SPos=cat(1,lindata.SPos);
sPos=SPos(:,1);

figure(12)
set(gcf,'color','w')
set(gca,'fontsize',16');
plot(sPos, d1, 'b-');
% atplotsyn(gca,RING1)
xlabel('s (m)')
ylabel('disp1 (m)')

[d0, dd0]= disp0(RING1, N1);
[d1, dd1]=disp1(RING1,d0,dd0,N1);
figure(12)
hold on
plot(sPos, d1, 'r')


%Function manipulating the scanning
function [alpha1all, S, d0, dd0, d1] = scan_disp1_list(RING1, N, omega_bound_lir, omega_bound_uir, n_chrom, index,Q1)
N1 = 1000; % for calculation of the dispersion - corresponds to the number of harmonics in the Fourier sum
a = omega_bound_lir;
b = omega_bound_uir;
S=[]; % matrix gathering all generated sextupole strengths

%Evolution o fhte first-order dispersion for all generated rings - updated during the scan
figure(30)
set(gcf,'color','w')
set(gca,'fontsize',16');
atplotsyn(gca,RING1')
xlabel('s (m)')
ylabel('disp1 (m)')
hold on
[lindata,tune,chrom]=atlinopt(RING1,0,1:length(RING1)+1);
SPos=cat(1,lindata.SPos);
sPos=SPos(:,1);

alpha1all = zeros(1,N); % vector of all alpha1
for j=1:N,
    str = ['\n Etape j = ', num2str(j)]; % mind the steps!
    fprintf(str)
    RING3=RING1;
    r = (b-a).*rand(n_chrom-2,1) + a;% random generation of omega
    s = Q1*r; % from omega to sextupoles
    S=[S,s]; %translated strengths of the sextupoles
    
    for i=1:n_chrom,
        RING3{index(i)}.PolynomB(1,3)= s(i)+ RING1{index(i)}.PolynomB(1,3);
    end
    [d0, dd0]= disp0(RING3, N1); % calculation of disp0 : suggest to drop this analytical calculation for a fast scan
    [d1, dd1]=disp1(RING3,d0,dd0,N1); % analytical calculation of the first-order dispersion
    alpha1all(j) = alf1(RING3, d0, dd0, d1); % calculation of alpha1
    figure(30)
    plot(sPos, d1, 'k'); % update of figure(30), comparing all first-order dispersion on the lattice
end

%From HOalpha.m, calculation of the first-order momentum compaction factor by sommation of the zeroth- and first-order dispersions and derivative.
 function alpha1 = alf1(RING, d0, dd0, d1)
        
        %Calcul des parametres lineaires de l'anneau - beta, phase/element
        alpha1=0;
        circumference=findspos(RING,length(RING)+1);
        dist=zeros(1,length(RING));
        for k=1:length(RING),
            dist(k)=findspos(RING,k+1)-findspos(RING,k);
        end
        %Sommation sur ring
        for j=1:length(RING),
            %Contribution des dipoles et dipoles combines
            if (strcmp(RING{j}.Class,'Bend')==1);
                rho = RING{j}.Length/RING{j}.BendingAngle;     
                alpha1=alpha1+(1/2*dd0(j)^2)*dist(j)+1/rho*d1(j)*dist(j);
                
                %Contribution des autres elements
            else
                alpha1=alpha1+1/2*dd0(j)^2*dist(j);
            end
        end
        alpha1=alpha1/circumference;
 

%function splitting the elements of the lattice for increased precision of the integrals 
function newelems=splitelem(elem,nslices)
%             nslices= 32;
newelems=atdivelem(elem,ones(1,nslices)./nslices,'KeepAxis');
