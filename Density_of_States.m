clear all;
warning off;
c=3e8;
hcut=(6.625e-34)/(2*pi);
lambdaB = input('Enter value of the wavelength(in microns) = ');                      
omega=2*pi*(c/lambdaB);
x1=input('Value of Mole Fraction = ');  %Value of Mole Fraction

%Coefficient Function for AlxGa1-xN
a11=9.827-8.216*x1-31.59*(x1^2);         %for AlxGa1-xN
b11=2.736+0.842*x1-6.293*(x1^2);         %for AlxGa1-xN

%Coefficient Function for GaN
a21=9.84;                       %for GaN
b21=2.74;                       %for GaN
  
%Bandgap Calculation of AlxGa1-xN/GaN 
eg11=6.28*x1+3.42*(1-x1)-1.3*x1*(1-x1);    %Bandgap of AlxGa1-xN
eg21=3.42;                                 %Bandgap of GaN

%Refractive Indices of the Layers
n1 =(a11.*((hcut*omega)/eg11)^(-2))*(2-sqrt(1+((hcut*omega)/eg11))-sqrt(1-((hcut*omega)/eg11)))+b11    %Refractive Index of AlxGa1-xN
n2 =(a21*((hcut.*omega)/eg21)^(-2))*(2-sqrt(1+((hcut*omega)/eg21))-sqrt(1-((hcut*omega)/eg21)))+b21    %Refractive Index of GaN
     
%Layer Thickness of each Material      
d1 = input('Layer Thickness of AlxGa1-xN(in microns) = ');        % Layer Thickness of AlxGa1-xN
d2 = input('Layer Thickness of GaN(in microns) = ');              % Layer Thickness of GaN

% for "S" Polarization
N=1000; % No. of Discretization points
u = linspace(0.15, 0.45, N); % Vector of Normalized Frequency
Ds1 = []; % DOS for "S" Polarization

for j=1:N
    o=u(j);
    
    Fs = @(x) x.*((n1^2*d1./sqrt(n1^2-x.^2)).*(sin(2*pi*d1*o.*sqrt(n1^2-x.^2)).*cos(2*pi*d2*o.*sqrt(n2^2-x.^2))+0.5.*(sqrt(n1^2-x.^2)./sqrt(n2^2-x.^2)+sqrt(n2^2-x.^2)./sqrt(n1^2-x.^2)).*sin(2*pi*d1*o.*sqrt(n1^2-x.^2)).*cos(2*pi*d2*o.*sqrt(n2^2-x.^2)))+...
        + (n2^2*d2./sqrt(n2^2-x.^2)).*(sin(2*pi*d2*o.*sqrt(n2^2-x.^2)).*cos(2*pi*d1*o.*sqrt(n1^2-x.^2))+0.5.*(sqrt(n2^2-x.^2)./sqrt(n1^2-x.^2)+sqrt(n1^2-x.^2)./sqrt(n2^2-x.^2)).*sin(2*pi*d2*o.*sqrt(n2^2-x.^2)).*cos(2*pi*d1*o.*sqrt(n1^2-x.^2)))-...
        - (1/o)*((n1^2-n2^2)^2/(4*pi)).*x.^2.*sin(2*pi*d1*o.*sqrt(n1^2-x.^2)).*sin(2*pi*d2*o.*sqrt(n2^2-x.^2))./(sqrt(n1^2-x.^2).*sqrt(n2^2-x.^2)).^3)./(1-(cos(2*pi*d1*o.*sqrt(n1^2-x.^2)).*cos(2*pi*d2*o.*sqrt(n2^2-x.^2))-0.5.*(sqrt(n1^2-x.^2)./sqrt(n2^2-x.^2)+...
        + sqrt(n2^2-x.^2)./sqrt(n1^2-x.^2)).*sin(2*pi*d1*o.*sqrt(n1^2-x.^2)).*sin(2*pi*d2*o.*sqrt(n2^2-x.^2))).^2);
    
    Ds1(j)=abs(quadl(Fs,0,n1-0.001));
    display(sprintf('Calculation for o[%d] is finished',j));
    
end

% for "P" Polarization
N=1000; % No. of Discretization Points
u = linspace(0.1, 0.35, N); % Vector of Normalized Frequency
Dp1 = [];  % DOS for "P" Polarization

for j=1:N
    o=u(j);    
    
     Fp = @(x) x.*((n1^2*d1./sqrt(n1^2-x.^2)).*(sin(2*pi*d1*o.*sqrt(n1^2-x.^2)).*cos(2*pi*d2*o.*sqrt(n2^2-x.^2))+0.5.*((n2^2/n1^2).*sqrt(n1^2-x.^2)./sqrt(n2^2-x.^2)+(n1^2/n2^2).*sqrt(n2^2-x.^2)./sqrt(n1^2-x.^2)).*sin(2*pi*d1*o.*sqrt(n1^2-x.^2)).*cos(2*pi*d2*o.*sqrt(n2^2-x.^2)))+...
         + (n2^2*d2./sqrt(n2^2-x.^2)).*(sin(2*pi*d2*o.*sqrt(n2^2-x.^2)).*cos(2*pi*d1*o.*sqrt(n1^2-x.^2))+0.5.*((n1^2/n2^2).*sqrt(n2^2-x.^2)./sqrt(n1^2-x.^2)+(n2^2/n1^2).*sqrt(n1^2-x.^2)./sqrt(n2^2-x.^2)).*sin(2*pi*d2*o.*sqrt(n2^2-x.^2)).*cos(2*pi*d1*o.*sqrt(n1^2-x.^2)))-...
         - (1/o)*((n1^2-n2^2)^2*(1/n1^2+1/n2^2-1)/(4*pi)).*x.^2.*sin(2*pi*d1*o.*sqrt(n1^2-x.^2)).*sin(2*pi*d2*o.*sqrt(n2^2-x.^2))./(sqrt(n1^2-x.^2).*sqrt(n2^2-x.^2)).^3)./(1-(cos(2*pi*d1*o.*sqrt(n1^2-x.^2)).*cos(2*pi*d2*o.*sqrt(n2^2-x.^2))-0.5.*((n2^2/n1^2).*sqrt(n1^2-x.^2)./sqrt(n2^2-x.^2)+...
         +(n1^2/n2^2).*sqrt(n2^2-x.^2)./sqrt(n1^2-x.^2)).*sin(2*pi*d1*o.*sqrt(n1^2-x.^2)).*sin(2*pi*d2*o.*sqrt(n2^2-x.^2))).^2);
    
    Dp1(j)=abs(quadl(Fp,0,n1-0.001));   
    display(sprintf('Calculation for o[%d] is finished',j));
end

% for "S" Polarization
subplot(2,1,1);
plot(u,Ds1,'b');
legend('S-Polarization');
xlabel('Normalized Frequency'); 
ylabel('Density of States(a.u.)');

% for "P" Polarization
subplot(2,1,2);
plot(u,Dp1,'r');
legend('P-Polarization');
xlabel('Normalized Frequency'); 
ylabel('Density of States(a.u.)');

% Add title to the Overall Plot
ha = axes ('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text (0.5, 1,'\bf Density of States of One Dimensional Photonic Crystal with Al_{x}Ga_{1-x}N/GaN Material Composition ','HorizontalAlignment','center','VerticalAlignment', 'top')

