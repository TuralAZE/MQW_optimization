function [E_cal_V,z_cal_V,first_wave_0V,second_wave_0V,E_basis1,E_basis2]=perturbation(Z,E,n,V,Wavefunction,V0,volt)
%%
e=1.602176487e-19;
width=400e-9;
for i=1:n
    for j=1:n
        z_0V(i,j)=Z(i,j);
        E_0V(i,j)=E(j)-E(i);%Energy diff in eV unit
        H_per(i,j)=-Z(i,j)*e*volt/width; %perturbed hamiltonian 원래 행렬성분[SI unit]
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E_basis1=H_per(1,1)/e;%변경된 절대 기준 넘겨주기(eV단위이다.)
E_basis2=0;
for j=2:n
    E_basis2=E_basis2+((H_per(j,1))^2)/e*E_0V(j,1); %에너지는 Electron Volt단위이다.
end
E_basis2=E_basis2/e; %Electron volt단위로 변환.
%%

first_E_0V = zeros(n,n);
second_E_0V = zeros(n,n);
first_z_0V = zeros(n,n);
second_z_0V = zeros(n,n);
first_wave_0V = zeros(length(V0),n);
second_wave_0V = zeros(length(V0),n);
for j=1:n %% 1st order perturbation of wave function.
    for i=1:n
        if ~(j==i)
            first_wave_0V(:,j) = first_wave_0V(:,j)+Wavefunction(:,i)*(H_per(j,i)/(e*E_0V(i,j)));%공식을 이용했다.
        end
    end
end
for i=1:n %% 2nd order perturbation of wave function.참고로.Z의 2차 order와 Energy의 2차 Correction까지는 파동함수의 1차 Correction까지만 고려되는것임.
    %파동함수의 2차 correction 의 개입이 시작되려면 최소 3차 correction이 필요하다.
    %에너지 간격이 20meV 보다 작은 difference가 일어나면 Z와 E의 3차 이상의 Correction을 고려하는 것이 좋음.. 
    for j=1:n
        if ~(j==i)
            for k=1:n
                if ~(k==i)
                    second_wave_0V = second_wave_0V+Wavefunction(:,j)*(H_per(j,k)*H_per(k,i)/(e^2*E_0V(j,i)*E_0V(k,i)));
                end
            end
            second_wave_0V = second_wave_0V-Wavefunction(:,j)*(H_per(j,i)*H_per(i,i)/((e*E_0V(j,i))^2));
            second_wave_0V = second_wave_0V-0.5*Wavefunction(:,i)*(H_per(i,j)*H_per(j,i)/((e*E_0V(j,i))^2));
        end
    end
end

for i=1:n
    for j=1:n



first_E_0V(i,j)=z_0V(j,j)-z_0V(i,i); %energy "difference" 에 대한 perturbation.우선 m unit으로 구하고 나중에 단위를 통일시킴.

%% 2nd order perturbation of energy "diff"
for k=1:n
    if ~(k==j)
        second_E_0V(i,j) = second_E_0V(i,j)+z_0V(k,j)^2/E_0V(k,j);
    end
end
for l=1:n
    if ~(l==i)
        second_E_0V(i,j) = second_E_0V(i,j)-z_0V(l,i)^2/E_0V(l,i);
    end
end
%%<m|Z|n>의 perturbation은..|n>을  1차 term까지 인경우랑..<m| 1차 term 까지인 경우를 더한것.

for k=1:n
    if ~(k==j)
        first_z_0V(i,j) = first_z_0V(i,j)+z_0V(j,k)*z_0V(k,i)/E_0V(k,j);
    end
end
for l=1:n
    if ~(l==i)
        first_z_0V(i,j) = first_z_0V(i,j)+z_0V(i,l)*z_0V(l,j)/E_0V(l,i);
    end
end


%%
for k=1:n
    for l=1:n
        if ~(k==j)&&~(l==i)
           second_z_0V(i,j) = second_z_0V(i,j)+z_0V(j,k)/E_0V(k,j)*z_0V(i,l)/E_0V(l,i)*z_0V(k,l);
        end
    end
end

for l=1:n
    if ~(l==j)
        second_z_0V(i,j) = second_z_0V(i,j)-z_0V(j,j)/E_0V(l,j)*z_0V(j,l)/E_0V(l,j)*z_0V(l,i);
    end
end

for k=1:n
    for l=1:n
        if ~(k==j)&&~(l==j)
           second_z_0V(i,j) = second_z_0V(i,j)+z_0V(k,l)/E_0V(l,j)*z_0V(j,k)/E_0V(k,j)*z_0V(l,i);
        end
    end
end

for k=1:n
    if ~(k==i)
        second_z_0V(i,j) = second_z_0V(i,j)-z_0V(i,i)/E_0V(k,i)*z_0V(i,k)/E_0V(k,i)*z_0V(j,k);
    end
end

for k=1:n
    for l=1:n
        if ~(k==i)&&~(l==i)
           second_z_0V(i,j) = second_z_0V(i,j)+z_0V(k,l)/E_0V(k,i)*z_0V(i,l)/E_0V(l,i)*z_0V(j,k);
        end
    end
end


%% 400nm로 나누고, 단위를 한번에 맞춘다.
E_cal_V(i,j,:)=E_0V(i,j)+first_E_0V(i,j)*(-2.5e6)*V+second_E_0V(i,j)*(-2.5e6)^2*V.^2;
z_cal_V(i,j,:)=z_0V(i,j)+first_z_0V(i,j)*(-2.5e6).*V+second_z_0V(i,j)*(-2.5e6)^2.*V.^2;

    end
end



end