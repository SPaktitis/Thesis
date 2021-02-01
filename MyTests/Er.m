%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H sunarthsh kanei kai ta duo( er,et ) kai sumferei.

function er = Er( Omega, nt, Dt )

    e=[];
    for ii=1:nt
        tmp1 = exp(-1i*2*pi*(ii-1)*Dt*Omega);
        e = [e;    tmp1];
    end
    er = 1/sqrt(nt) .* e;

end