function [fourierSpectrum2D] = iDFT_selfMade(signal2D)
% Definición de tamaño de la matriz y la matriz de salida
[Ny,Nx] = size(signal2D);
fourierSpectrum2D = zeros([Ny Nx]);

% Se pasa la matriz a decimales para que no ocurran problemas de operacion
% con complejos
signal2D = double(signal2D);

% Se recorren cada una de las posiciones de la matriz
for p = 1:Nx
    for q = 1:Ny
      fourierSpectrum2D(q,p) = fourierSum(signal2D,p,q,Nx,Ny);
    end
end

end

function indexSum = fourierSum (signal2D, p, q, Nx, Ny)
% Equivalente a las sumatorias dobles, solo que un poquito más rapido
% ya que se evita el uso de ciclos for
[X,Y] = meshgrid(1:Nx,1:Ny);
matrixExp = (1/(Nx*Ny))*exp(1i*2*pi*((X*p/Nx) + (Y*q/Ny)));
indexSum = sum(signal2D.*matrixExp,'all');


end