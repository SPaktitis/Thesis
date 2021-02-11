%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H sunarthsh kanei kai ta duo( er,et ) kai sumferei.

function er = Er( Omega, nt, Dt )

    e=[];
    for ii=1:nt
        e = [e;    exp(-1i*2*pi*(ii-1)*Dt*Omega) ];
    end
    er = 1/sqrt(nt) .* e;

end