function [acc,ave_ene,virial_coeff] = compute_force(pos,epsilon,BoxSize,DIM,N)

    %compute forces on positions using LJ-potential
    Sij = zeros(1,DIM); %Box scaled units
    Rij = zeros(1,DIM); %Real space units
    
    %set all variables to zero
    ene_pot = zeros(1,N);
    acc = zeros(N,DIM);
    virial = 0.0;
    Rcutoff = 2.5;
    phicutoff = 4.0/(Rcutoff^12)-4.0/(Rcutoff^6);
    
    %loop over all pairs of particles
    for i = 1:N-1
        for j = i+1:N
            Sij(:) = pos(i,:) - pos(j,:);
            for l = 1:DIM   %periodic interaction
                if abs(Sij(l)) > 0.5
                    Sij(1,l) = Sij(1,l) - sign(Sij(1,l))*1.0;
                end
            end
            Rij = BoxSize*Sij; %scale the box to the real units in this case reduced LJ units
            Rsqij = sum(Rij.*Rij);
            if Rsqij < Rcutoff^2.0
                %calculate LJ potential inside cutoff
                %we calculate parts of the LJ potential at a time to
                %improve the efficiency of the computation
                   rm2 = 1.0/Rsqij;            %1/r^2
                   rm6 = rm2^3.0;           %1/r^6
                   rm12 = rm6^2.0;            %1/r^12
                   phi =  epsilon*(4.0*(rm12-rm6)-phicutoff);          %4[1/r^12-r^6]-phi(Rc)  using shift LJ potential
                % the following is dphi = -(1/r)(dV/dr)
                   dphi =  epsilon*24.0*rm2*(2.0*rm12-rm6);          %24[2/r^14-1/r^8]
                   ene_pot(i) =  ene_pot(i) + 0.5*phi;           %accumulate energy
                   ene_pot(j) =  ene_pot(j) + 0.5*phi;           %accumulate energy
                
                   virial = virial - dphi*Rsqij;            %virial is needed to calculate the pressure
                   acc(i,:) = acc(i,:)+dphi*Sij;           %accumulate force;
                   acc(j,:) = acc(j,:)-dphi*Sij;            %Fji = -Fij
                
            end
        end
    end
    
    ave_ene = sum(ene_pot)/N;
    virial_coeff = -virial/DIM;
    

end