clear;clc
close all

%% initilize particles

DIM = 2;
N = 32;

BoxSize = 10.0; %6.35
global volume
volume = BoxSize^DIM;
global density 
density = N/volume;

load position.dat
pos = position;
pos = pos/BoxSize;

MassCenter = sum(pos)/N;
for i=1:DIM
    pos(:,i) = pos(:,i)-MassCenter(i);
end

figure
for i=1:N
    plot(pos(i,1),pos(i,2),'o','MarkerSize',10,'MarkerFaceColor','b');
    text(pos(i,1),pos(i,2),num2str(i));
    hold on
end
axis off


%% setting up the simulation
NSteps = 10000;
deltat = 0.0032; %time step in reduced time units
TRequested = 0.5; %reduced temperature
DumpFreq = 100; % save the position to file every DumpFreq steps
epsilon = 1.0;

% main MD 
[ene_kin_aver,ene_pot_aver,temperature,pressure,pos]= MD(pos,NSteps,deltat,TRequested,DumpFreq,epsilon,BoxSize,DIM);

animation(NSteps/DumpFreq);
%% post-process

figure
subplot(4,1,1);
plot(ene_kin_aver,'k-');
ylim([0.485,0.514]);
ylabel('E_{K}');
subplot(4,1,2);
plot(ene_pot_aver,'k-');
ylim([-1.4,0]);
ylabel('E_{P}');
subplot(4,1,3);
plot(temperature,'k-');
ylim([0.485, 0.515]);
ylabel('T')
subplot(4,1,4)
plot(pressure,'k-');
ylim([-0.6,1.0])
ylabel('P')






