function F0 = F0n(RING,N)
RING1=RING;
F0 = zeros(1,N);

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
            ds(ix)=(findspos(RING1,ix+1)-findspos(RING1,ix))/rho;
        end
 end

for k = 1:N,
    x = cos((k-1).*phi).*sqrt(betax);
    y = ds;
   F0(k) = y*x*nux/(pi);
end
F0(1)=F0(1)/2;