warning off
clear all
clc;
lambdaB1=linspace(1.4,1.7,100);
L1 = 30;                            % grating length in microns [Code can be changed here]
z = [0:0.1:L1];                     % propagation distance
lambda =1.55;                       % input wavelength in microns
c=3e8;
kappa = 0.01;                       %[Code can be changed here]
omega=2*pi*(c/lambda);
hcut=(6.625e-34)/(2*pi);
x=0.3;                              %mole fraction [Code can be changed here]
 
a1=9.827-8.216*x-31.59*(x^2)        %for AlxGa1-xN
b1=2.736+0.842*x-6.293*(x^2)        %for AlxGa1-xN
 
a2=9.84;                            %for GaNi
b2=2.74;                            %for GaN
 
eg1=6.28*x+3.42*(1-x)-1.3*x*(1-x)   %for AlxGa1-xN
eg2=3.42;                           %for GaN
 
n1 =(a1.*((hcut*omega)/eg1)^(-2))*(2-sqrt(1+((hcut*omega)/eg1))-sqrt(1-((hcut*omega)/eg1)))+b1;        %for AlxGa1-xN
n2 =(a2*((hcut.*omega)/eg2)^(-2))*(2-sqrt(1+((hcut*omega)/eg2))-sqrt(1-((hcut*omega)/eg2)))+b2;        %for GaN

neff1 = n2-n1                        % mode effective index

% coupling coefficient in 1/microns
dbeta1=(1./lambdaB1-1/lambda).*(2*pi*neff1);       % delta beta
S1 = sqrt(kappa^2-dbeta1.^2);
K1= (1/lambda)*(2*pi*neff1);
r1=kappa.*(sinh(L1.*S1))./((dbeta1.*(sinh(L1.*S1)))-(1i*S1.*cosh(L1.*S1)));
R1=r1.*r1;

%To plot the graph
plot(lambdaB1,abs(R1),'k')
xlabel ('Wavelength [\mum]');
ylabel ('Reflectivity');

% Add title to the Overall Plot
ha = axes ('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text (0.5, 1,'\bf Reflectivity in One Dimensional Photonic Crystal with Al_{x}Ga_{1-x}N/GaN Material Composition ','HorizontalAlignment','center','VerticalAlignment', 'top')


