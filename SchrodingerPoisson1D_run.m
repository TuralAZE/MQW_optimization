%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% layers structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first column is the conduction band offset in eV
% second column is the length of the layer in nm
% third column is the n doping volumique of that layer in 1e18cm-3 

% You have to put a resonable amount of doping! Otherwise, it will diverge 

clear all
close all
clc

%%


Epsi   = 4.32;    %dielectric constant of quantum well
Epsi0=8.854e-12;
e=1.602e-19;

parmeter_STO =[
    0       0       0.97       0
    1.56    2.55e20 0.999   0.103
    3.70    6.04e20 1.071   0.112
    4.17    6.80e20 1.093   0.119
    6.25    1.02e21 1.162   0.149
    12.5    2.04e21 1.254   0.222
    25.00   4.07e21 1.329   0.334
];                       % La % in STO, charge_density, effective_mass, chemical_potential

La = 0:0.01:25; 
charge_density_STO = interp1(parmeter_STO(:,1),parmeter_STO(:,2),La );
effective_mass_STO = interp1(parmeter_STO(:,1),parmeter_STO(:,3),La );
chmical_potential_STO = interp1(parmeter_STO(:,1),parmeter_STO(:,4),La );

La_STO = 3;               % La % in STO

STO=chmical_potential_STO(1,find(La==La_STO));       %Quantum well (eV)
LAO=2.4;    %Quantum barrier band offset (eV)

meff_LAO   = 1.4;       %Effective mass
meff_STO   = effective_mass_STO(1,find(La==La_STO));

Lattice_LAO = 0.379;    %Lattice constant
Lattice_STO = 0.3905;


Nd=charge_density_STO(1,find(La==La_STO));  %cm-3


M=[
LAO     Lattice_LAO* 5   0
STO     Lattice_STO* 2   Nd*1e-18
LAO     Lattice_LAO* 1   0
STO     Lattice_STO* 1   0
LAO     Lattice_LAO* 5   0
];              %unit set of MQWs composition


% Lw = Lattice_LAO*5;
% Lb = Lattice_STO*2;
% 
% DF   = 10;        % Electrical field discontinuity [MV/cm]
% Ls = 0.2;         % doping spike thickness [nm] (in order to get the E-field)
% Fb = +DF*(Lw+2*Ls)/(Lw+Lb+4*Ls)*1e6*1e2;   %[V/m]
% Fw = -DF*(Lb+2*Ls)/(Lw+Lb+4*Ls)*1e6*1e2;   %[V/m]
% 
% dopS = DF*1e6*1e2*Epsi*Epsi0/e;   % charge/m2 MUST BE added on the interface for GaN/AlN for Wurtzite
% dopV = dopS/(Ls*1e-9);            % charge/m3
% dopV = dopV*1e-6*1e-18;           % charge 1e18cm-3
% 
% M=[
% LAO     Ls   0
% LAO     Lw - 2*Ls   0
% LAO     Ls   0
% 
% STO     Ls   -dopV
% STO     Lb - 2*Ls   163
% STO     Ls   +dopV
% 
% LAO     Ls   0
% LAO     Lw - 2*Ls   0
% LAO     Ls   0
% 
% STO     Ls   -dopV
% STO     Lb - 2*Ls   163
% STO     Ls   +dopV
% 
% LAO     Ls   0
% LAO     Lw - 2*Ls   0
% LAO     Ls   0
% ];              %unit set of MQWs composition


%%
Nloops = 50;                  % number of loops
n      = 3;                   % number of solution asked per model
dz     = 1e-11;               % resolution of the grid [m]

    meff = meff_STO;
%     meff = 0.07;


 [z,Ec,psic,V0,Vtot,Ef]=SchrodingerPoisson1D_f(M,Epsi,meff,STO,LAO,n,dz,Nloops);


    Wavefunction = zeros(length(V0),n);
    Energy = zeros(n,1);

    Wavefunction=psic;
    Energy=Ec;

% for m=1:n+1
% 
% if m==1
%     meff = meff_STO;
% %     meff = 0.07;
% else
%     meff = round(sum((V0(1,:)./LAO.*meff_LAO+meff_STO)'.*abs(psic(:,m-1).*psic(:,m-1).*dz)),4);
% end
% 
%  [z,Ec,psic,V0,Vtot,Ef]=SchrodingerPoisson1D_f(M,Epsi,meff,STO,LAO,n,dz,Nloops);
% 
%  if m==1
%     Wavefunction = zeros(length(V0),n);
%     Energy = zeros(n,1);
%     else
%     Wavefunction(:,m-1)=psic(:,m-1);
%     Energy(m-1)=Ec(m-1);
%  end
%  
% end

%% dipole_element

Z = zeros(n,n);
E = zeros(n,n);
psicxpsic = zeros(n,n);

for i=1:n
    for j=1:n
        Z(i,j)=sum((Wavefunction(:,i).*Wavefunction(:,j).*z(1,:)'.*dz),'all');
        
        E(i,j)=Energy(i)-Energy(j);
        
        psicxpsic(i,j)=sum((Wavefunction(:,i).*Wavefunction(:,j).*dz),'all');
    end
end

%%

figure('position',[100 100 1000 700],'color','w');
subplot(1,1,1,'fontsize',15)
hold on;grid on;
col=colormap(jet);

ScF    = 0.1;                 % scaling factor to plot the wave function [Without Dimension]
F0     = -0e5;%-2e7;%-6e7;       % Electric field [Volt/meter]
T      = 300;                 % Temperature [Kelvin], react on the Fermi function only

for i=1:n
    PSIc(:,i)=abs(Wavefunction(:,i)).^2/max(abs(Wavefunction(:,i)).^2)*ScF + Energy(i); % normalisation for the plotting
end


    grid off
    set(gca,'color',col(1,:))
    shading flat
    hcb=colorbar;
    title(hcb,'\fontsize{8}cm-3')
    
    plot(z*1e9,V0,  'w--','linewidth',1)
    plot(z*1e9,Vtot,'w-' ,'linewidth',1)
    
    plot([z(1) z(end)]*1e9,[1 1]*Ef,'g','linewidth',1)
    text(z(end)*1e9*0.95,Ef+0.01,'\color{green}Fermi')

for i=1:n
    plot(z*1e9,PSIc(:,i),'color','r','linewidth',1)
end

xlabel('z (nm)')
ylabel('Energy (eV)')

title(strcat('\fontsize{12}T=',num2str(T),'K; meff=',num2str(meff),'; Epsilon=',num2str(Epsi),'; dz=',num2str(dz*1e9),'nm;'))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Perturbation

V=-4:0.08:4;

[E_cal_V,z_cal_V] = perturbation(Z,E,n,V);

figure(2);
plot(V,squeeze(z_cal_V(1,2,:)),V,squeeze(z_cal_V(2,3,:)),V,squeeze(z_cal_V(1,3,:))); xlabel('voltage'); ylabel('dipole element[m]');

figure(3);
plot(V,squeeze(z_cal_V(1,2,:)).*squeeze(z_cal_V(2,3,:)).*squeeze(z_cal_V(1,3,:))); xlabel('voltage'); ylabel('coupled dipole element[m^3]');

figure(4);
plot(V,squeeze(E_cal_V(1,2,:)),V,squeeze(E_cal_V(2,3,:)),V,squeeze(E_cal_V(1,3,:))); xlabel('voltage'); ylabel('Energy[eV]');

