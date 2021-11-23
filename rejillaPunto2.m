function [rejilla] = rejillaPunto2(Nx,Ny,dx,dy,m, L)

x = 1:Nx;
y = 1:Ny;
[X,Y] = meshgrid(x,y);

X = X*dx;
Y = Y*dy;

rejilla = 0.5 * (1+m*cos(2*pi*X/L));
end