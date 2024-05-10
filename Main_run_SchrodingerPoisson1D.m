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



%%  InGaAs/AlInAs MQWs Input parameter
mode=2;

InGaAs=0;
AlInAs=0.52;
AlInAs2=0.2;

meff_AlInAs   = 0.076;
meff_InGaAs   = 0.0426;
meff = meff_InGaAs; 
Epsi   = 11.6;
bw = 0.5869;    % basic width of the layer

M=[
AlInAs         bw* 7     0           meff_AlInAs
InGaAs         bw* 10    1.5         meff_InGaAs
AlInAs         bw* 2     0           meff_AlInAs
InGaAs         bw* 7     0           meff_InGaAs
AlInAs         bw* 2     0           meff_AlInAs
InGaAs         bw* 6     0           meff_InGaAs
AlInAs         bw* 3     0           meff_AlInAs
InGaAs         bw* 4     0           meff_InGaAs
AlInAs         bw* 7     0           meff_AlInAs
];              %unit set of MQWs composition
adjustment_factor = 1.00; % Increase widths by 10%
M(:, 2) = M(:, 2) * adjustment_factor;

%% BZO/BSO MQWs Input parameter
%mode=4;

%Epsi   = 4.4;    % dielectric constant of quantum well / Physical properties of transparent perovskite oxides (Ba,La)SnO3 with high electrical mobility at room temperature 
% Epsi0=8.854e-12;
% e=1.602e-19;
% 
% parmeter_BSO =[
%     0       0       0.19    0
%     1       1.6e19  0.21    0.1
%     3       8.4e19  0.25    0.4
%     7       4.4e20  0.30    0.8
% ];                       % La % in BSO, charge_density, effective_mass, (chemical_potential)X from Supplemental Material: Electron effective mass and mobility limits in degenerate perovskite stannate BaSnO3
% 
% La = 0:0.01:25; 
% charge_density_BSO = interp1(parmeter_BSO(:,1),parmeter_BSO(:,2),La );
% effective_mass_BSO = interp1(parmeter_BSO(:,1),parmeter_BSO(:,3),La );
% chemical_potential_BSO = interp1(parmeter_BSO(:,1),parmeter_BSO(:,4),La );
% 
% La_BSO_1 = 0.72;               % La % in BSO
% 
% La_BSO_2 = 0.72;               % La % in BSO
% 
% BSO_1=chemical_potential_BSO(1,find(La==La_BSO_1));       %Quantum well (eV)
% BSO_2=chemical_potential_BSO(1,find(La==La_BSO_2));       %Quantum well (eV)
% BZO=4.524-2.799;    %Quantum barrier band offset (eV)
% 
% % Effects of Internal Relaxation under Inplane Strain on the Structural, Electronic and Optical Properties of Perovskite BaZrO3
% % Electronic band structure and optical phonons of BaSnO3 and Ba0.97La0.03SnO3 single crystals: Theory and experiment
% 
% meff_BZO   = 0.99;       %Effective mass / Electronic band-offsets across Cu2O/BaZrO3 heterojunction and its stable photo-electro-chemical response: First-principles theoretical analysis and experimental optimization
% meff_BSO_1   = effective_mass_BSO(1,find(La==La_BSO_1));
% meff_BSO_2   = effective_mass_BSO(1,find(La==La_BSO_2));
% meff = meff_BSO_1;
% 
% Lattice_BZO = 0.4194;    %Lattice constant
% Lattice_BSO = 0.4116;
% 
% 
% Nd_1=charge_density_BSO(1,find(La==La_BSO_1));  %cm-3
% Nd_2=charge_density_BSO(1,find(La==La_BSO_2));  %cm-3
% 
% 
% M=[
% BZO         Lattice_BZO* 7     0           meff_BZO
% BSO_1       Lattice_BSO* 8      Nd_1*1e-18  meff_BSO_1
% BZO         Lattice_BZO* 1      0           meff_BZO
% BSO_2       Lattice_BSO* 3      Nd_2*1e-18  meff_BSO_2
% BZO         Lattice_BZO* 7     0           meff_BZO
% BSO_2       Lattice_BSO* 3      Nd_2*1e-18  meff_BSO_2
% BZO         Lattice_BZO* 7     0           meff_BZO
% ];              %unit set of MQWs composition



%% resolution setting

Nloops = 100;                  % number of loops
n      = 5;                   % number of solution asked per model
dz     = 1e-11;               % resolution of the grid [m]
count=1;                      % 파라미터 스윕 위한 횟수 설정

%% Schrodinger_Poisson calculation    
if mode==2         % for the precise result compared to Wang program
newWavefunc={};
newEnergy={};
newEc={};
newz={};
for ind= 1:count    
%바꾸고싶은 파라미터를 지정하기.
%예시
 [z,Ec,psic,V0,Vtot,Ef,Ntot,me]=SchrodingerPoisson1D_f(M,Epsi,meff,InGaAs,AlInAs,n,dz,Nloops);
 %[z,Ec,psic,V0,Vtot,Ef]=SchrodingerPoisson1D_f(M,Epsi,meff,BSO_1,BZO,n,dz,Nloops);

    Wavefunction = zeros(length(V0),n);
    Energy = zeros(n,1);
    Wavefunction=psic;
    Energy=Ec;
    newWavefunc{end+1}=Wavefunction;
    newEnergy{end+1}=Energy;
    newEc{end+1}=Ec;
    newz{end+1}=z;
end
end

%% dipole_element
Z = zeros(n,n,count);
E = zeros(count,n,n);
psicxpsic = zeros(n,n,count);
for ind=1:count

for i=1:n
    for j=1:n
        disp(newWavefunc);
        tempwf=cell2mat(newWavefunc(ind));
        tempz=cell2mat(newz(ind));
        nz=squeeze(tempz);
        B=squeeze(tempwf);
        Z(i,j,ind)=sum((B(:,i).*B(:,j).*nz(1,:)'.*dz),'all');
        %ind는 parameter의 인덱스 값으로, 각 index안의 값이 스윕햇을때 dipole element값.
        tempe=cell2mat(newEnergy(ind));
        E(ind,i,j)=tempe(i)-tempe(j);
        %psicxpsic는 파동함수의 내적들로, 비대각성분은 0이어야만함!
        psicxpsic(i,j,ind)=sum((B(:,i).*B(:,j).*dz),'all');
    end
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

%z=z';
%Vtot=Vtot';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Perturbation
tempEc=cell2mat(newEnergy);
E12 = tempEc(2,:) - tempEc(1,:);
E13 = tempEc(3,:) - tempEc(1,:);
E14 = tempEc(4,:) - tempEc(1,:);
%얻고싶은 에너지 차이 값들을 위와 같이 지정해서 구하세요.
% 400nm에 건 전압을 계산하는 것이다.(2nd order까지 계산. 너무 크다면 3차항 이상도 필요해지니 조심)
voltstep=0.08;

volt=4;%max 값을 넘기세요.
V=-volt:voltstep:volt;%E level 이 서로 Cross한 뒤부터는 더이상 Perturbation Theory가 성립하지 않으니 조심해야함. FWM의 경우에는 E level 이 너무 붙어있으므로 Cross하기 쉽다.
[E_cal_V,z_cal_V,first_wave_0V,second_wave_0V,E_basis1,E_basis2] = perturbation(Z,E,n,V,Wavefunction,V0,volt);

figure(2);
plot(V,squeeze(z_cal_V(1,2,:)),V,squeeze(z_cal_V(2,3,:)),V,squeeze(z_cal_V(3,4,:))); xlabel('voltage'); ylabel('dipole element[m]');
legend('Z12','Z23','Z13');

figure(3);
plot(V,squeeze(z_cal_V(1,2,:)).*squeeze(z_cal_V(2,3,:)).*squeeze(z_cal_V(3,4,:))); xlabel('voltage'); ylabel('coupled dipole element[m^3]');
figure(4);
plot(V,squeeze(E_cal_V(1,2,:)),V,squeeze(E_cal_V(2,3,:)),V,squeeze(E_cal_V(3,4,:))); xlabel('voltage'); ylabel('Energy[eV]');
legend('E12','E23','E13');
want = 0;
while want~=-99 %Energy level cross 일어나면, 절대로 Perturbation Theory가 맞지않음.상전이가 일어난 상태라서 해밀토니안을 재정의해야함.
    want= input('plot하고자 하는 Voltage값을 입력하세요.,그만하고자 하면 -99를 입력하십시오.')
    if want==-99
        break;
    end
    k= find(abs(V-want) < (voltstep)/2.1); %에너지를 잘못치더라도, 가까운 에너지 격자 값을 찾기 위함임.
    A = (V(k)/V(end)) % Perturbation을 구할때 Vmax로 구하므로, 우리가 원하는 Volt에 맞는 값들은 람다^order가 곱해져야함.
    E_basis=A*E_basis1+A^2*E_basis2+Energy(1); %새로운 바닥 에너지 상태는, 원래 바닥에너지+바닥에너지의 섭동량.
    new_2ndE=E_cal_V(1,2,k)+E_basis;
    WaveCorr= Wavefunction;
    WaveCorr = WaveCorr+ A*first_wave_0V + (A^2)*second_wave_0V; %Wave Function의 2nd order perturbation Theory.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Here, I re-define the energy grid in order optimize the meshing
    dE1=1e-4; dE2=1e-2;

    new_E1 = E_basis:dE1:E_basis+0.1 ;
    new_E2 = new_E1(end):dE2:max(Vtot);
    new_En=sort([new_E1 new_E2]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fermime=dot(WaveCorr(:,1).*WaveCorr(:,1).*dz,me);%새로 고친 파동함수로 질량의 기대값을 다시 고친다.
    new_Ec=zeros(n);
    for i=1:n
        new_Ec(i)=E_basis+squeeze(E_cal_V(1,i,k));%새로 고쳐진 에너지 레벨을 다시구한다.
    end
    %Fermi level 을 구할때에 Effective mass 는 Ground state를 이용하여 mass의 기대값을
    %구하는 형태로 진행함.
    h    = 6.62606896E-34;              %% Planck constant [J.s]
    hbar = h/(2*pi);
    e    = 1.602176487E-19;             %% electron charge [C]
    m0   = 9.10938188E-31;              %% electron mass [kg]
    Epsi0= 8.854187817620E-12;          %% Vaccum dielectric constant [F/m]
    kB   = 1.3806488E-23;               %% Boltzmann's constant [J/K]

    ro=[];
    for i=1:length(new_Ec)
        ro( new_En>new_Ec(i),i) = e*fermime*m0/(pi*hbar^2);
        ro( new_En<new_Ec(i),i) = 0;
    end

    [Ef_new,NN_new,roEf_new]=find_Ef_f(new_Ec,new_En,ro,Ntot,T);
    New_V0=-((V(k)*2.5e6)*z)+V0;%포텐셜 그래프에 전기장을 도입하기 위해 400nm로 나누고 Voltage를 곱한다.(eV단위) 전자는 음전하를 띄므로 마이너스부호. 가장 왼쪽이 V=0의 기준이다.
    New_Vtot=-((V(k)*2.5e6)*z)+Vtot;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('position',[100 100 1000 700],'color','w');
subplot(1,1,1,'fontsize',15)
hold on;grid on;
col=colormap(jet);
    for i=1:n
    newPSIc(:,i)=abs(WaveCorr(:,i)).^2/max(abs(WaveCorr(:,i)).^2)*ScF + new_Ec(i); % normalisation for the plotting
end


    grid off
    set(gca,'color',col(1,:))
    shading flat
    hcb=colorbar;
    title(hcb,'\fontsize{8}cm-3')
    
    plot(z*1e9,New_V0,  'w--','linewidth',1)
    plot(z*1e9,New_Vtot,'w-' ,'linewidth',1)
    
    plot([z(1) z(end)]*1e9,[1 1]*Ef_new,'g','linewidth',1)
    text(z(end)*1e9*0.95,Ef_new+0.01,'\color{green}Fermi')

for i=1:n
    plot(z*1e9,newPSIc(:,i),'color','r','linewidth',1)
end

xlabel('z (nm)')
ylabel('Energy (eV)')

title(strcat('\fontsize{12}T=',num2str(T),'K; meff=',num2str(meff),'; Epsilon=',num2str(Epsi),'; dz=',num2str(dz*1e9),'nm;'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end



%Susceptibility plotting

% Constants
h    = 6.62606896E-34;              %% Planck constant [J.s]
hbar = h/(2*pi);
e    = 1.602176487E-19;             %% electron charge [C]
e0 = 8.854187817e-12;               % Vacuum permittivity in F/m
Ne = 5e24;                          % Electron density in m^-3, adjust accordingly

%Dipole and Energy Values
% Assuming Z and E matrices are filled out correctly from the code above
Z12 = Z(1,2,1); % Dipole element between state 1 and 2
Z23 = Z(2,3,1); % Dipole element between state 2 and 3
Z34 = Z(3,4,1); % Dipole element between state 3 and 4
Z41 = Z(4,1,1); % Dipole element between state 4 and 1
E120 = E12*e
E130 = E13*e
E140 = E14*e
omega12 = E120/hbar; % Transition frequency between state 1 and 2
omega13 = E130/hbar; % Transition frequency between state 1 and 3
omega14 = E140/hbar; % Transition frequency between state 1 and 4
% Damping Factors
gamma12=0.01039547*e/hbar; % Damping factor for transition 1-2
gamma13=0.0156*e/hbar; % Damping factor for transition 1-3
gamma14=0.0148*e/hbar; % Damping factor for transition 1-4
% Frequency Range for Plotting
omega = linspace(0.01*e/hbar, 0.15*e/hbar, 1000); % Define appropriately
% Calculate chi3 for Each Frequency
chi3 = zeros(size(omega));
for idx = 1:length(omega)
    chi3(idx) = (Ne * e^4 * Z12 * Z23 * Z34 * Z41) / ...
                (e0 * hbar^3 * ...
                (omega(idx) - omega12 - 1i*gamma12) * ...
                (2*omega(idx) - omega13 - 1i*gamma13) * ...
                (3*omega(idx) - omega14 - 1i*gamma14));
end
% Plot Real and Imaginary Parts of Susceptibility
figure;
plot(omega, real(chi3), 'r', 'DisplayName', 'Real(\chi^{(3)})');
hold on;
plot(omega, imag(chi3), 'b--', 'DisplayName', 'Imag(\chi^{(3)})');
xlabel('Angular Frequency (rad/s)');
ylabel('Third-order Susceptibility (\chi^{(3)})');
legend;
title('Third-order Optical Susceptibility vs Angular Frequency');
figure;
plot(omega, abs(chi3), 'g', 'DisplayName', '|χ^{(3)}|'); % using 'g' for green line
xlabel('Angular Frequency (rad/s)');
ylabel('Magnitude of Third-order Susceptibility (|χ^{(3)})|');
legend;
title('Magnitude of Third-order Optical Susceptibility vs Angular Frequency');

