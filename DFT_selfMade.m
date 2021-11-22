function [fourierSpectrum2D] = DFT_selfMade(signal2D)
[Ny,Nx] = size(signal2D);

fourierSpectrum2D = zeros([Ny Nx]);

for p = 1:Nx
    for q = 1:Ny
      fourierSpectrum2D(q,p) = fourierSum(signal2D,p,q,Nx,Ny);
    end
end

end

function indexSum = fourierSum (signal2D, p, q, Nx, Ny)
[X,Y] = meshgrid(1:Nx,1:Ny);
matrixExp =exp(-1i*2*pi*((X*p/Nx) + (Y*q/Ny)));
indexSum = sum(double(signal2D).*matrixExp,'all');


end
