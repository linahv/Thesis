function F2 = F2n(RING,disp0,disp1,N)
RING1=RING;
F2 = zeros(1,N);

[lindata,tune,chrom]=atlinopt(RING1,0,1:length(RING1)+1);
beta=cat(1,lindata.beta);
betax=beta(:,1);
muxy=cat(1,lindata.mu);
nux=muxy(length(RING1)+1,1)/2/pi;
phi=muxy(:,1)/nux;
ds = zeros(1,length(RING1)+1);
 for ix=1:length(RING1),
         if (strcmp(RING1{ix}.Class,'Bend')==1),
             rho = RING1{ix}.Length/RING1{ix}.BendingAngle;
             K = RING1{ix}.PolynomB(1,2);
             L = RING1{ix}.Length;
             ds(ix)=(K-4*K/rho*disp0(ix))*(findspos(RING1,ix+1)-findspos(RING1,ix));
         elseif (strcmp(RING1{ix}.Class,'Quadrupole')==1),
             K1 = RING1{ix}.PolynomB(1,2);
             L = RING1{ix}.Length;
             ds(ix)=+K1*(findspos(RING1,ix+1)-findspos(RING1,ix));
         elseif (strcmp(RING1{ix}.Class,'Sextupole')==1),
             K2 = RING1{ix}.PolynomB(1,3);
             L = RING1{ix}.Length;
             ds(ix)= -2*K2*disp0(ix)*(findspos(RING1,ix+1)-findspos(RING1,ix));
         end
 end

for k = 1:N,
    x = cos((k-1).*phi).*sqrt(betax).*disp1';
    y = ds;
   F2(k) = y*x*nux/(pi);
end
F2(1)=F2(1)/2;