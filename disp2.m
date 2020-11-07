function disp2 = disp2(RING,disp0,disp1,N)
% Wiedemann's method of calculation of the higher order dispersion
RING1=RING;
%Fourier harmonics
F2 = F2n(RING1,disp0,disp1,N);
%Twiss parameters and phase
[lindata,tune,chrom]=atlinopt(RING1,0,1:length(RING1)+1);
beta=cat(1,lindata.beta);
betax=beta(:,1);
muxy=cat(1,lindata.mu);
nux=muxy(length(RING1)+1,1)/2/pi;
phi=muxy(:,1)/nux;
%% Initialisation
disp2 = zeros(1, length(RING1)+1);
x = zeros(1,N);
for ix=1:length(RING1)+1,
    for k=1:N,
    x(k) = cos((k-1).*phi(ix))/(nux^2-(k-1)^2);
    end
    y =F2;
    disp2(ix) = (x*y')*sqrt(betax(ix));
end
disp2 = disp2 - disp1;
%% Comparaison with the polynomial fit
DisplayFlag = 1;
DDP=1e-3;
N=20;
LRING=length(RING1)+1;
delta = linspace(-DDP, +DDP,N);
delta4eval = linspace(-DDP, +DDP, 2*N+1);
eta = zeros(N,LRING);
Deta = zeros(N,LRING);

REFPTS=ones(length(RING)+1,1);
dDP = 1e-6;

for k =1:length(delta)
    for i=1:LRING
    [orbP,o1P]=findorbit4(RING1,delta(k)+0.5*dDP,i);
    [orbM,o1M]=findorbit4(RING1,delta(k)-0.5*dDP,i);
    dispersion=(orbP-orbM)/dDP;
    eta(k,i)=dispersion(1,:);
    Deta(k,i)=dispersion(2,:);
    end
end

% fit curve upto 3rd order
porder =3;
eta2=zeros(LRING,1);
for i=1:LRING,
    pvalue = polyfit(delta, eta(:,i)', porder);
    pvalue_Deta = polyfit(delta, Deta(:,i)', porder);
    peta = polyval(pvalue, delta4eval);
    pDeta = polyval(pvalue_Deta, delta4eval);
    eta2(i) = pvalue(end-2)/3;    
end
% 
SPos=cat(1,lindata.SPos);
s=SPos(:,1);
figure(5)
set(gcf,'color','w')
set(gca,'fontsize',16');
plot(s, eta2, 'k')
hold on
plot(s, disp2,'r.')
atplotsyn(gca,RING1)
xlabel('s (m)')
ylabel('Third-order dispersion (m)')