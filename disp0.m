function [disp0,disp0p] = disp0(RING,N)
% Wiedemann's method of calculation of the higher order dispersion
RING1=RING;

%Fourier harmonics
F0 = F0n(RING1,N);

%Twiss parameters and phase
[lindata,tune,chrom]=atlinopt(RING1,0,1:length(RING1)+1);
beta=cat(1,lindata.beta);
betax=beta(:,1);
alpha=cat(1,lindata.alpha);
alphax=alpha(:,1);
muxy=cat(1,lindata.mu);
nux=muxy(length(RING1)+1,1)/2/pi; %normalisation by 2pi required, A.T. convention
phi=muxy(:,1)/nux; %normalisation by nux required
disp= cat(2,lindata.Dispersion);
dispx=disp(1,:)'; %for validation
Ddispx=disp(2,:)';%for validation

disp0  = zeros(1, length(RING1)+1);
disp0p = zeros(1, length(RING1)+1);


x = zeros(1,N);
z = zeros(1,N);
for ix=1:length(RING1)+1,
    for k=1:N,
    x(k) = cos((k-1).*phi(ix))/(nux^2-(k-1)^2);
    z(k) = -(k-1)*sin((k-1).*phi(ix))/(nux^2-(k-1)^2)/nux./betax(ix);
    end
    y =F0;
    disp0(ix)  = (x*y')*sqrt(betax(ix));
    disp0p(ix) = (z*y')*sqrt(betax(ix));
end
disp0p = disp0p -alphax'./betax'.*disp0;

%Comparison with the first-order dispersion from Twiss, and plots
SPos=cat(1,lindata.SPos);
s=SPos(:,1);
figure(1)
set(gcf,'color','w')
set(gca,'fontsize',16');
plot(s,disp0,'r.')
hold on
plot(s,dispx,'k')
atplotsyn(gca,RING1)
xlabel('s (m)')
ylabel('First-order dispersion (m)')

figure(2)
set(gcf,'color','w')
set(gca,'fontsize',16');
plot(s,disp0p,'r.')
hold on
plot(s,Ddispx,'k')
atplotsyn(gca,RING1)
xlabel('s (m)')
ylabel('First-order dispersion derivative (rad)')