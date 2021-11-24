function[F]=FZeroFunction(p,q,u1)
    Fsuma=0;

[M,N]=size(u1);

    for n=1:N
        for m=1:M
            Fsuma=Fsuma+u1(m,n)*exp(-2j*pi*((p*n/N)+(q*m/M)));
       
        end
    end

 F=Fsuma;
 
end