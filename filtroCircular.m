function [filtroResultado] = filtroCircular(paresDeCentros,radios,filterSize)

% Funcion para hacer filtros de naturaleza circular
[numCircles, ~] = size(paresDeCentros);

h = filterSize(1,1);
w = filterSize(1,2);


[x,y] = meshgrid(1:w,1:h);

filtro = false(filterSize);
for i = 1:numCircles
    circulo = (x-paresDeCentros(i,1)).^2 + (y-paresDeCentros(i,2)).^2 <= radios(i)^2;
    filtro = circulo|filtro;
end

filtroResultado= double(filtro);

end