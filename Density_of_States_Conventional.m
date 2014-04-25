clear all
warning off
clc

%Refractive Indices of the Layers
n1 = input('Refractive Index of 1st Dielectric Material = ');
n2 = input('Refractive Index of 2nd Dielectric Material = ');

%%% Normalized Thicknesses of the Layers
d1 = n2/(n1+n2);
d2 = 1-d1;

% for "S" Polarization
N=input('No. of Discretization points = '); %No. of Discretization points
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
u = linspace(0.15, 0.45, N); % Vector of Normalized Frequency
Dp1 = [];  % DOS for "P" Polarization

for j=1:N
    o=u(j);
    
   Fp = @(x) x.*((n1^2*d1./sqrt(n1^2-x.^2)).*(sin(2*pi*d1*o.*sqrt(n1^2-x.^2)).*cos(2*pi*d2*o.*sqrt(n2^2-x.^2))+0.5.*((n2^2/n1^2).*sqrt(n1^2-x.^2)./sqrt(n2^2-x.^2)+(n1^2/n2^2).*sqrt(n2^2-x.^2)./sqrt(n1^2-x.^2)).*sin(2*pi*d1*o.*sqrt(n1^2-x.^2)).*cos(2*pi*d2*o.*sqrt(n2^2-x.^2)))+...
        + (n2^2*d2./sqrt(n2^2-x.^2)).*(sin(2*pi*d2*o.*sqrt(n2^2-x.^2)).*cos(2*pi*d1*o.*sqrt(n1^2-x.^2))+0.5.*((n1^2/n2^2).*sqrt(n2^2-x.^2)./sqrt(n1^2-x.^2)+(n2^2/n1^2).*sqrt(n1^2-x.^2)./sqrt(n2^2-x.^2)).*sin(2*pi*d2*o.*sqrt(n2^2-x.^2)).*cos(2*pi*d1*o.*sqrt(n1^2-x.^2)))-...
        - (1/o)*((n1^2-n2^2)^2*(1/n1^2+1/n2^2-1)/(4*pi)).*x.^2.*sin(2*pi*d1*o.*sqrt(n1^2-x.^2)).*sin(2*pi*d2*o.*sqrt(n2^2-x.^2))./(sqrt(n1^2-x.^2).*sqrt(n2^2-x.^2)).^3)./(1-(cos(2*pi*d1*o.*sqrt(n1^2-x.^2)).*cos(2*pi*d2*o.*sqrt(n2^2-x.^2))-0.5.*((n2^2/n1^2).*sqrt(n1^2-x.^2)./sqrt(n2^2-x.^2)+...
        + (n1^2/n2^2).*sqrt(n2^2-x.^2)./sqrt(n1^2-x.^2)).*sin(2*pi*d1*o.*sqrt(n1^2-x.^2)).*sin(2*pi*d2*o.*sqrt(n2^2-x.^2))).^2);
    
    Dp1(j)=abs(quadl(Fp,0,n1-0.0001)); 
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
text (0.5, 1,'\bf Density of States of One Dimensional Photonic Crystal ','HorizontalAlignment','center','VerticalAlignment', 'top')


