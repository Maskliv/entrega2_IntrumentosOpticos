function[u2]=propFresnel(u1,L,lambda,z) 
% propagation - impulse response approach 
% uniform sampling 
% u1 - Transmitancia
% L - Longitud de la transmitancia
% lambda - Longitud de onda
% z - Distancia de propagación
% u2 - Campo en el plano de observación
 
[M,N]=size(u1);  
dx=L/M; 
k=2*pi/lambda; %Numero de onda 

x=-L/2:dx:L/2-dx; %Plano coordenado
[X,Y]=meshgrid(x,x); 
 
h=1/(1i*lambda*z)*exp(1i*k/(2*z)*(X.^2+Y.^2)); %FuciónDeImpulso
H=fft2(fftshift(h))*dx^2; 
U1=fft2(fftshift(u1));  
U2=H.*U1; 
u2=ifftshift(ifft2(U2));
end