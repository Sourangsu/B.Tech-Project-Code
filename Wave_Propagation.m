warning off
clear all
clc
L = 30;                      % grating length in microns
z = [0:0.01:L];              % propagation distance
lambda = 1.55;               % input wavelength in microns
lambdaB = 1.5;               % Bragg wavelength in microns

c=3e8;
omega=2*pi*(c/lambda);
hcut=(6.625e-34)/(2*pi);
x=0.25;                      % mole fraction [Code can be varied here]
 
a1=9.827-8.216*x-31.59*(x^2)         %for AlxGa1-xN
b1=2.736+0.842*x-6.293*(x^2)         %for AlxGa1-xN
 
a2=9.84;                             %for GaN
b2=2.74;                             %for GaN
 
eg1=6.28*x+3.42*(1-x)-1.3*x*(1-x)    %for AlxGa1-xN
eg2=3.42;                            %for GaN
 
n1 =(a1.*((hcut*omega)/eg1)^(-2))*(2-sqrt(1+((hcut*omega)/eg1))-sqrt(1-((hcut*omega)/eg1)))+b1;        %for AlxGa1-xN
n2 =(a2*((hcut.*omega)/eg2)^(-2))*(2-sqrt(1+((hcut*omega)/eg2))-sqrt(1-((hcut*omega)/eg2)))+b2;        %for GaN

neff1 = n2-n1                             % mode effective index
kappa = 0.01;                             % coupling coefficient in 1/microns [Code can be varied here]
dbeta1=(1./lambda-1/lambdaB)*2*pi*neff1   % delta beta
S1 = sqrt(kappa^2-dbeta1.^2)
K1= (1/lambdaB)*2*pi*neff1

% To find the backward (aa) and forward (b) field amplitudes inside grating
f1 = ((1*kappa.*exp(-j.*dbeta1.*z))./(1*dbeta1.*sinh(S1.*L)-j.*S1.*cosh(S1*L))).*sinh(S1.*(z-L));
r1 = ((1*exp(-j.*dbeta1.*z))./(1*dbeta1.*sinh(S1.*L)-j.*S1.*cosh(S1*L))).*(dbeta1.*sinh(S1.*(z-L))+j*S1.*cosh(S1.*(z-L)));

% To plot the backward (aa) and forward (b) field amplitudes inside grating
plot(z,abs(f1).^2, 'r')
xlabel ('Grating Length');
ylabel ('Forward Wave');
% Add title to the Overall Plot
ha = axes ('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text (0.5, 1,'\bf Study of Wave Propagation in One Dimensional Photonic Crystal with Al_{x}Ga_{1-x}N/GaN Material Composition ','HorizontalAlignment','center','VerticalAlignment', 'top')

figure
plot(z,abs(r1).^2)
xlabel ('Grating Length');
ylabel ('Backward Wave');
% Add title to the Overall Plot
ha = axes ('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text (0.5, 1,'\bf Study of Wave Propagation in One Dimensional Photonic Crystal with Al_{x}Ga_{1-x}N/GaN Material Composition ','HorizontalAlignment','center','VerticalAlignment', 'top')

