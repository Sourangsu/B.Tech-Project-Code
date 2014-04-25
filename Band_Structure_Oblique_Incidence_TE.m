warning off
clear all
clc

%The dielectric constants of the Composition 
epsa=8.8;        %AlGaN
epsb=1.46^2; 
epsc=9.7;        %GaN
a=0.75;          %um 
ra=0.4;          %um 
rb=0.05;         %um 
rc=0.25/2;       %um 
 
a1=a; 
b1=2*pi/a; 
 
n=input('Input N= ');  
NumberofPW=(2*n+1); 
 
count=1; 
G=(-n:n).*b1; 
r=(-n:n); 
N=4*n+1; 
mc=round(N*rc/a); 
mb=round(N*rb/a); 
ma=N-2*mc-2*mb; 

eps1=[epsc*ones(mc,1);epsb*ones(mb,1); 
   epsa*ones(ma,1);epsb*ones(mb,1);epsc*ones(mc,1)]; 
eps20=real(fftshift(fft(eps1)./N)); 
 
 
for x=1:NumberofPW, 
  	for y=x:NumberofPW, 
        b=r(x)-r(y)+2*n+1; 
        eps2(x,y)=eps20(b); 
        eps2(y,x)=eps2(x,y);    
   end 
end 
 
ky=0.5*pi/a; 
k0=(-pi/a:2*pi/a/30:pi/a)+i*ky;  
eps2=inv(eps2); 
counter=1; 
for ii=1:length(k0), 
   k=k0(ii); 

	M=abs(k+G.')*abs(k+G).*eps2; %TM wave
	
	E=sort(abs(eig(M))); 
	freq(:,counter)=sqrt(abs(E(1:10))).*a./2./pi; 
	display(sprintf('calculation of k=%f is finished',k)); 
	counter=counter+1; 
end 
tmpx=1:length(k0); 
plot(tmpx,freq,'linewidth',1) 
xlabel('Wave Vector') 
ylabel('wa/2\pic') 

% Add title to the Overall Plot
ha = axes ('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text (0.5, 1,'\bf Band Structure of One Dimensional Photonic Crystal under Oblique Incidence (TE Mode) ','HorizontalAlignment','center','VerticalAlignment', 'top')

