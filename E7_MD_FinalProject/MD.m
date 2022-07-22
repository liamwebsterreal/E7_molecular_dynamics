function [ene_kin_aver,ene_pot_aver,temperature,pressure,pos] = ...
        MD(pos,NSteps,deltat,TRequested,DumpFreq,epsilon,BoxSize,DIM)
    
    global density
    global volume
    
    % vectors to store parameters values at each step
    N = size(pos,1);
    ene_kin_aver = zeros(1,NSteps);
    ene_pot_aver = zeros(1,NSteps);
    temperature = zeros(1,NSteps);
    virial = zeros(1,NSteps);
    pressure = zeros(1,NSteps);
%     ene_pot = zeros(1,N);
    
    vel = randn(N,DIM)-0.5;
    acc = randn(N,DIM)-0.5;
    
    %open file which we will save the outputs to
    fid = fopen('traj.xyz','w');
    
    for k = 1:NSteps
        % refold positions according to periodic boundary conditions
        for i = 1:DIM
            period = find(pos(:,i)>0.5);
            pos(period,i) = pos(period,i)-1.0;
            period = find(pos(:,i)<-0.5);
            pos(period,i) = pos(period,i)+1.0;
        end
        %r(t+dt) modify positions accoding to velocity and acceleration
        pos = pos + deltat*vel + 0.5*(deltat^2.0)*acc;  %step 1
        
        %calculate temperature
        [ene_kin_aver(k),temperature(k)] = Calculate_Temperature(vel,BoxSize,DIM,N);
        
        %rescale velocities and take half step
        chi = sqrt(TRequested/temperature(k));
%         chi = 1.0;
        vel = chi*vel + 0.5*deltat*acc;%step 2
        
        %compute forces a(t+dt), ene_pot, virial
        [acc,ene_pot_aver(k),virial(k)] = compute_force(pos,epsilon,BoxSize,DIM,N);%step 3
        
        %complete the velocity step
        vel = vel + 0.5*deltat*acc; %step 4
        
        %calculate temperature
        [ene_kin_aver(k),temperature(k)] = Calculate_Temperature(vel,BoxSize,DIM,N);
        
        %calculate pressure
        pressure(k) = density*temperature(k) + virial(k)/volume;
        
        %print output to file every DumpFreq number of steps
        if (mod(k,DumpFreq)==0)
            fprintf(fid,'%4i \n',N); %write the number of particles to file
            %write all of quantites to this step to the file
            fprintf(fid,'total energy = %12E  Temperature = %12E \n',ene_kin_aver(k)+ene_pot_aver(k),temperature(k));
            for n=1:N %write positions to file
                fprintf(fid,'X ');
                for l=1:DIM
                    fprintf(fid,'%12E ',pos(n,l)*BoxSize);
                end
                fprintf(fid,'\n');
            end
            if (DIM==2)
                figure
                for i=1:N
                    plot(pos(i,1)*BoxSize,pos(i,2)*BoxSize,'o','MarkerFaceColor','b');
                    hold on
                end
                xlim([-0.5*BoxSize,0.5*BoxSize]);
                ylim([-0.5*BoxSize,0.5*BoxSize]);
                print(2,'-dbmp',sprintf('images/%d',k));
                close;
            end
        end
        
    end
    fclose(fid);
end