function mol_fit = kA_fit(time,b_e,a_m,a_s,b_s,b_m,kA,kx,kx_s,b_ms,b_p,molecules0)

opts = odeset('RelTol',1e-3);

[t,mol_fit] = ode15s(@(t,molecules) KA_ODE(time,molecules,a_s, b_s, a_m, b_m, kA, kx, kx_s,b_ms, b_e, b_p,molecules0),time,molecules0,opts);


end