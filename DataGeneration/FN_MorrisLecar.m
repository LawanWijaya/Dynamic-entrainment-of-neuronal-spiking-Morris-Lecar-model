% Function for the data generation code
% Morris-Lecar differential equations
function z = FN_MorrisLecar(t, y, item,input, gK, gCa, gL, EK, EL, ECa, C, V1, V2, V3, V4, phi, Amp, L,I0)


% global V1 V2 V3 V4 gL EL gK EK gCa ECa phi Amp L C

V=y(1); n=y(2);

switch item

    case 1

        I = I0;

    case 2

        t0 = input(1);
        if (t>=t0 && t<=t0+L)
            I = I0+Amp;
        else
            I = I0;
        end

    case 3

        t0 = input(1); t1 = input(2);

        if (t>=t0 && t<=t0+L) || (t>=t1 && t<=t1+L)
            I = I0+Amp;
        else
            I = I0;
        end

end

minf=1/2*(1+tanh((V-V1)/V2));
taun=1/cosh((V-V3)/(2*V4));
ninf=1/2*(1+tanh((V-V3)/V4));

z(1)=1/C*(I-gL*(V-EL)-gK*n*(V-EK)-gCa*minf*(V-ECa));
z(2)=phi*(ninf-n)/taun;

z = z';
end