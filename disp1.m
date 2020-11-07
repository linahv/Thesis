function [disp1x, disp1p] = disp1(RING,disp0,disp0p,N)
% Wiedemann's method of calculation of the higher order dispersion

RING1=RING;

%Fourier harmonics
F1 = F1n(RING1,disp0,N);

%Twiss parameters and phase
[lindata,tune,chrom]=atlinopt(RING1,0,1:length(RING1)+1);
beta=cat(1,lindata.beta);
betax=beta(:,1);
alpha=cat(1,lindata.alpha);
alphax=alpha(:,1);
muxy=cat(1,lindata.mu);
nux=muxy(length(RING1)+1,1)/2/pi;%normalisation by 2pi required, A.T. convention
phi=muxy(:,1)/nux; %normalisation by nux required

disp1x = zeros(1, length(RING1)+1);
disp1p = zeros(1, length(RING1)+1);
x = zeros(1,N);
z = zeros(1,N);
for ix=1:length(RING1)+1,
    for k=1:N,
    x(k) = cos((k-1).*phi(ix))/(nux^2-(k-1)^2);
    z(k) = -(k-1)*sin((k-1).*phi(ix))/(nux^2-(k-1)^2);
    end
    y =F1;
    disp1x(ix)  = (x*y')*sqrt(betax(ix));
    disp1p(ix) = (z*y')/nux./sqrt(betax(ix));
end
disp1x = disp1x - disp0;
disp1p = disp1p - disp0p -alphax'./betax'.*(disp0+disp1x);

%% Comparaison with the polynomial fit
DDP=1e-3;
N=20;
LRING=length(RING1)+1;
delta = linspace(-DDP, +DDP,N);
delta4eval = linspace(-DDP, +DDP, 2*N+1);
Deta = zeros(N,LRING);
eta = zeros(N,LRING);
dDP = 1e-6;
%Variation of the orbit
for k =1:length(delta)
    for i=1:LRING
    [orbP,o1P]=findorbit4(RING1,delta(k)+0.5*dDP,i);
    [orbM,o1M]=findorbit4(RING1,delta(k)-0.5*dDP,i);
    dispersion=(orbP-orbM)/dDP;
    eta(k,i) =dispersion(1,:);
    Deta(k,i)=dispersion(2,:);
    end
end

% fit curve up to the  3rd order
porder =3;
Deta1=zeros(LRING,1);
eta1=zeros(LRING,1);
for i=1:LRING,
    pvalue_eta = polyfit(delta, eta(:,i)', porder);
    pvalue_Deta = polyfit(delta, Deta(:,i)', porder);
    peta = polyval(pvalue_eta, delta4eval);
    pDeta = polyval(pvalue_Deta, delta4eval);
    eta1(i) = pvalue_eta(end-1)/2; 
    Deta1(i) = pvalue_Deta(end-1)/2;    
end

SPos=cat(1,lindata.SPos);
s=SPos(:,1);
figure(3)
set(gcf,'color','w')
set(gca,'fontsize',16');
plot(s,disp1x,'r.')
hold on
plot(s, eta1, 'k')
atplotsyn(gca,RING1)
ylabel('Second-order dispersion (m)')
xlabel('s (m)')

figure(4)
set(gcf,'color','w')
set(gca,'fontsize',16');
plot(s,disp1p,'r.')
hold on
plot(s, Deta1, 'k')
atplotsyn(gca,RING1)
ylabel('Second-order dispersion derivative (rad)')
xlabel('s (m)')
Dispp = zeros(1, length(RING1)+1);
for ix=2:length(RING1),
    ds = findspos(RING1,ix+1)-findspos(RING1,ix-1);
    Dispp(ix) = (disp1x(ix+1)-disp1x(ix-1))/ds;
end

SPos=cat(1,lindata.SPos);
s=SPos(:,1);
figure(24)
plot(s, Dispp,'k')
hold on
plot(s, disp1p,'r.')
