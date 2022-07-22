function [ene_kin_aver,temperature] = Calculate_Temperature(vel,BoxSize,DIM,N)

    ene_kin = 0.0;
    
    for i = 1:N
        real_vel = BoxSize*vel(i,:);
        ene_kin = ene_kin + 0.5*sum(real_vel.*real_vel);
    end
    
    ene_kin_aver = 1.0*ene_kin/N;
    temperature = 2.0*ene_kin_aver/DIM;
    
end