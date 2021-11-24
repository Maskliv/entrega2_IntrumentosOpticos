function[campoPropagado] = transformadaFresnel(campoEntrada,dx,dy,distPropagacion, waveLength, options)
arguments
    campoEntrada
    dx double
    dy double
    distPropagacion double
    waveLength double
    options.dft logical = false
end


% propagation - impulse response approach 
% uniform sampling 
% campoEntrada - Transmitancia
% L - Longitud de la transmitancia
% lambda - Longitud de onda
% z - Distancia de propagación
% u2 - Campo en el plano de observación
 
% Se define el espacio coordenado
[Ny,Nx] = size(campoEntrada);

x = 1-(Nx/2):Nx/2;
y = 1-(Ny/2):Ny/2;
[X,Y] = meshgrid(x,y);

dfx = (1/(Nx*dx));
dfy = (1/(Ny*dy));

% Damos dimensiones de mundo
fX = X*dfx;
fY = Y*dfy;

% Definimos constantes
k = 2*pi/waveLength;

z = distPropagacion;

impulsoEnFourier = exp(1i*k*z)*exp(-1i*pi*waveLength*z*((fX.^2) + (fY.^2))); % transformada de fourier Fución De respuesta al Impulso

if options.dft
    U0Fourier = (dx*dy)*fftshift(DFT_selfMade(campoEntrada));
    UzFourier = U0Fourier.*impulsoEnFourier;
    campoPropagado = (dfx*dfy)*iDFT_selfMade(UzFourier);
else
    U0Fourier = (dx*dy)*fftshift(fft2(campoEntrada));
    UzFourier = U0Fourier.*impulsoEnFourier;
    campoPropagado = (dfx*dfy)*ifft2(UzFourier);
end


end