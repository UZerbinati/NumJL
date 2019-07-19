function driftExample()
    f0(x) = if (-1 ≤ x ≤  1) return 0.5 else return 0 end
    D = [-4,4];
    x = D[1]:0.01:D[2];
    T = 0.0;
    while(true)
        PyPlot.clf();
        title(string(L"Drift Equation Evolution $\mu =$",mean(drift.(f0,x,T))))
        plot(x,drift.(f0,x,T));
        T = T+0.1;
        msg=readline();
        if msg=="close"
            return;
        end
    end
end
function transportExample()
    f0(x) = if (-1 ≤ x ≤  1) return 0.5 else return 0 end
    D = [-4,4];
    x = D[1]:0.01:D[2];
    T = 0.0;
    while(true)
        PyPlot.clf();
        plot(x,transport.(f0,1,x,T));
        T = T+0.1;
        msg=readline();
        if msg=="close"
            return;
        end
    end
end
