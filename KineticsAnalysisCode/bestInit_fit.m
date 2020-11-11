function mol_fit = bestInit_fit(time,a_s, b_s, k_init, k_elon, k_elon_prime, b_m, k_on, k_off, kx, kx_s, b_ms, b_e, b_p,k_nuc,b_nuc,n_plasmids,k_off_nuc,P,molecules0)

opts = odeset('RelTol',1e-3);

[t,mol_fit] = ode15s(@(t,molecules) Bestinit_ODE(time,molecules,a_s, b_s, k_init, k_elon, k_elon_prime, b_m, k_on, k_off, kx, kx_s, b_ms, b_e, b_p,k_nuc,b_nuc,n_plasmids,k_off_nuc,P),time,molecules0,opts);


end