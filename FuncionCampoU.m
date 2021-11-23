function[u]=FuncionCampoU(n,m,z,F0,dx1,lambda)

    u_suma=0;
   

C=size(u1);
P=C(1,2);

    for p=1:P
        for q=1:P
            F=F0(p,q)*exp(-1j*pi*dx1^2.*lambda*z*(p^2 +q^2 ) );
            
            u_suma=u_suma+F*exp(-1j*2*pi/P*(p*n+q*m) );
        end
    end
u=u_suma*dx1.^2;
end