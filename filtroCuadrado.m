function [filtroResultado] = filtroCuadrado(paresDeCentros,lados,filterSize)

% Funcion para hacer filtros de naturaleza circular
[numRects, ~] = size(paresDeCentros);

h = filterSize(1,1);
w = filterSize(1,2);


[x,y] = meshgrid(1:w,1:h);

filtro = false(filterSize);
for i = 1:numRects
    rectangulo = abs(x-paresDeCentros(i,1))<=lados(i,1)/2 & abs(y-paresDeCentros(i,2))<=lados(i,2)/2;
    filtro = rectangulo|filtro;
end

filtroResultado = double(filtro);
end
