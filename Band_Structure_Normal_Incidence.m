%%% This program is designed for calculating the band structure of one-dimensional photonic
%%% crystal under normal electromagnetic incident wave

warning off
clear all
clc 
epsa=input('Dielectric Constant for Material 1 = '); %Dielectric Constant for Material 1(Atom)
epsb=input('Dielectric Constant for the Material 2 = '); %Dielectric Constant for the Material 2(Background Material)
a=1.0;   % actual period of the PBG material
R=0.2*a; % width of the atom in a period
a1=a;
b1=2*pi/a;
n=input('N = ');  %number of grids in the reciprocal space along the +x direction (the same number is used for –x direction)
NumberofPW=(2*n+1);
count=1;
G=(-n:n).*b1;  %grid vectors in the reciprocal space
r=(-n:n);
N=2*NumberofPW+1;
m=floor(N*R/a)+1; 
eps1=[epsb*ones((N-m)/2,1);  %matrix containing the dielectric function in one period in the real space;
epsa*ones(m,1);
epsb*ones((N-m)/2,1)];
eps20=(fftshift(fft(eps1)./N));  %Fourier transformed matrix of eps1;

for x=1:NumberofPW,
 for y=x:NumberofPW,
 b=r(x)-r(y)+(2*n+1)+1;
 eps2(x,y)=eps20(b);
 eps2(y,x)=eps2(x,y);
 end
end

k0=-pi/a:2*pi/a/30:pi/a;
eps2=inv(eps2);  % matrix (epsilon^-1).(G -G')
counter=1;

for ii=1:length(k0),
 k=k0(ii);
 M=abs((k+G')*(k+G)).*(eps2);
E=sort(abs(eig(M)));
freq(:,counter)=sqrt(abs(E(1:10))).*a./2./pi;
counter=counter+1;
end

%Plot for the Graph
tmpx=1:length(k0);
plot(tmpx,freq)
xlabel('Wave Vector')
ylabel('w_{a}/2\pic')

% Add title to the Overall Plot
ha = axes ('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text (0.5, 1,'\bf Band Structure of One Dimensional Photonic Crystal under Normal Incidence ','HorizontalAlignment','center','VerticalAlignment', 'top')

