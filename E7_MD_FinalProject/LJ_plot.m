clear;clc
close all

r = linspace(0.01,3.0,500);
epsilon = 1.0;
sigma=1.0;
for i=1:length(r)
    E_LJ(i) = 4*epsilon*((sigma/r(i))^12-(sigma/r(i))^6);
end

Rcutoff = 2.0;
phicutoff = 4.0/(Rcutoff^12)-4.0/(Rcutoff^6);
E_LJ_shift = E_LJ - phicutoff;

plot(r,E_LJ,'r-','LineWidth',1);
hold on
plot(r,E_LJ_shift,'b-','LineWidth',1);
xlim([0,3.0]);
ylim([-1.5,1.5])
legend('LJ potnetial','shift LJ potential')
xlabel('r/ \sigma')
ylabel('E_{LJ}/\epsilon')
