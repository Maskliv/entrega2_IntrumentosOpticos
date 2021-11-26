function[campoPropagado] = transformadaFresnel(campoEntrada,dx,dy,distPropagacion, waveLength, options)
arguments
    campoEntrada
    dx double
    dy double
    distPropagacion double
    waveLength double
    options.dft logical = false
    options.zoom logical = false
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

x0 = 1-(Nx/2):Nx/2;
y0 = 1-(Ny/2):Ny/2;
[X0,Y0] = meshgrid(x0,y0);


z = distPropagacion;

% Damos dimensiones de mundo
X = X0*(waveLength*z/(Nx*dx)); % Condicion de muestreo entre planos
Y = Y0*(waveLength*z/(Ny*dy));

X0 = X0*dx;
Y0 = Y0*dy;




% Definimos constantes
k = 2*pi/waveLength;



frenteParabolicoU0 = exp(1i*k/(2*z)*((X0.^2) + (Y0.^2))); % transformada de fourier Fución De respuesta al Impulso
frenteParabolicoUz = (exp(1i*k*z)/1i*waveLength*z)*exp(1i*k/(2*z)*((X.^2) + (Y.^2)));

if options.dft
    Uprima = campoEntrada.*frenteParabolicoU0;
    UdoblePrima = (dx*dy)*fftshift(DFT_selfMade(Uprima));
    campoPropagado = frenteParabolicoUz.*UdoblePrima;
else
    Uprima = campoEntrada.*frenteParabolicoU0;
    UdoblePrima = (dx*dy)*fftshift(fft2(Uprima));
    campoPropagado = frenteParabolicoUz.*UdoblePrima;
end

if options.zoom
    hMin = Ny/2 - round((Ny*dy/2)/(waveLength*z/(Ny*dy)));
    hMax = Ny/2 + round((Ny*dy/2)/(waveLength*z/(Ny*dy)));
    wMin = Nx/2 - round((Nx*dx/2)/(waveLength*z/(Nx*dx)));
    wMax = Nx/2 + round((Nx*dx/2)/(waveLength*z/(Nx*dx)));
    campoPropagado = campoPropagado(hMin:hMax, wMin:wMax);
end



end