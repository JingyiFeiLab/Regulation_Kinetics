function mol_fit = nosRNA_fit(x0,time,a_s,b_s,ka,kd,kx_s,b_m,b_ms,b_e,b_p,molecules0)


a_m = x0(1);
kx = x0(2);
%b_m = x0(3);

%kx = x0;

opts = odeset('RelTol',1e-3);
[~,mol_fit] = ode15s(@(t,molecules) sRNA_ODE(time,molecules,0, 0, a_m, b_m, ka, kd, kx, kx_s, b_ms, b_e, b_p),time,molecules0,opts);

end