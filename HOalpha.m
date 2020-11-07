function [alpha0,alpha1,alpha2]=HOalpha(RING, N)
% N : number of Fourier harmonics for the expansion. N>1000 for a precise enough value.

%Slicing of the ring for the integral calculation
ring_fun = RING;
 for k=1:length(RING),
                elem=splitelem(RING{length(RING)-k+1});
                ring_fun=cat(1,ring_fun(1:length(RING)-k),elem,ring_fun(length(RING)-k+2:length(ring_fun)));
 end
ring=ring_fun;
[d0, dd0]= disp0(ring, N);
[d1, dd1]=disp1(ring,d0,dd0,N);
d2(:,1)= disp2(ring,d0,d1,N);
alpha0=alf0(ring)
alpha1=alf1(ring)
alpha2=alf2(ring)

    function alpha0 = alf0(RING)
        %Calculation of the Twiss parameters
        alpha0=0;
        circumference=findspos(RING,length(RING)+1);
        dist=zeros(1,length(RING));
        for k=1:length(RING),
            dist(k)=findspos(RING,k+1)-findspos(RING,k);
        end
        %Integral over the ring
        for j=1:length(RING),
            %Contribution of dipoles and combined-function magnets
            if (strcmp(RING{j}.Class,'Bend')==1);
                rho = RING{j}.Length/RING{j}.BendingAngle;
                alpha0=alpha0+1/rho*d0(j)*dist(j);
            end
        end
        alpha0=alpha0/circumference;
    end
	
    function alpha1 = alf1(RING)
        %Calculation of the Twiss parameters
        alpha1=0;
        circumference=findspos(RING,length(RING)+1);
        dist=zeros(1,length(RING));
        for k=1:length(RING),
            dist(k)=findspos(RING,k+1)-findspos(RING,k);
        end
        %Integral over the ring
        for j=1:length(RING),
            %Contribution of dipoles and combined-function magnets
            if (strcmp(RING{j}.Class,'Bend')==1);
                rho = RING{j}.Length/RING{j}.BendingAngle;     
                alpha1=alpha1+(1/2*dd0(j)^2)*dist(j)+1/rho*d1(j)*dist(j);
                
            %Contribution of other elements
            else
                alpha1=alpha1+1/2*dd0(j)^2*dist(j);
            end
        end
        alpha1=alpha1/circumference;
    end

    function alpha2 = alf2(RING)
       %Calculation of the Twiss parameters
        alpha2=0;
        circumference=findspos(RING,length(RING)+1);
        dist=zeros(1,length(RING));
        for k=1:length(RING),
            dist(k)=findspos(RING,k+1)-findspos(RING,k);
        end
         %Integral over the ring
        for j=1:length(RING),
            %Contribution of dipoles and combined-function magnets
            if (strcmp(RING{j}.Class,'Bend')==1);
                rho = RING{j}.Length/RING{j}.BendingAngle;
                
                alpha2=alpha2+dd0(j)*dd1(j)*dist(j)+(-1/rho*d0(j)*dd0(j)^2/2+1/rho*d2(j)*dist(j));
            %Contribution of other elements
            else
                alpha2=alpha2+dd0(j)*dd1(j)*dist(j);
            end
        end
        alpha2=alpha2/circumference;
    end
end
%Split function
function newelems=splitelem(elem)  
            nslices= 30;
            newelems=atdivelem(elem,ones(1,nslices)./nslices,'KeepAxis');
end

