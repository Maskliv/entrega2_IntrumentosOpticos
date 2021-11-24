function [campoPropagado] = espectroAngular(campoEntrada,dx,dy,distPropagacion, waveLength, options)
arguments
    campoEntrada
    dx double
    dy double
    distPropagacion double
    waveLength double
    options.dft logical = false
end

% Se define el espacio coordenado
[Ny,Nx] = size(campoEntrada);

x = 1-(Nx/2):Nx/2;
y = 1-(Ny/2):Ny/2;
[X,Y] = meshgrid(x,y);

% Damos dimensiones de mundo
fX = X*(1/(Nx*dx));
fY = Y*(1/(Ny*dy));

% Definimos constantes
k = 2*pi/waveLength;

z = distPropagacion;

freqPropagator = exp(1i*k*z*sqrt(1-(waveLength^2)*((fX.^2)+(fY.^2))));

if options.dft
    % Paso 1
    A_0 = fftshift(DFT_selfMade(campoEntrada));
    % Paso 2
    A_z = A_0.*freqPropagator;
    % Paso 3
    campoPropagado = iDFT_selfMade(A_z);
else
    % Paso 1
    A_0 = fftshift(fft2(campoEntrada));
    % Paso 2
    A_z = A_0.*freqPropagator;
    % Paso 3
    campoPropagado = ifft2(A_z);
end



end