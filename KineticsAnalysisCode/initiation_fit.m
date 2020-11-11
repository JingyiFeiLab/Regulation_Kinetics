function mol_fit = initiation_fit(x0,time,k_elon,b_m,b_p,n_plasmids,molecules0)


k_i = x0(1);
kx = x0(2);
%b_m = x0(3);

%kx = x0;

opts = odeset('RelTol',1e-3);
[~,mol_fit] = ode15s(@(t,molecules) nosRNA_ODE(time,molecules, k_i, kx, b_m, k_elon, b_p,n_plasmids),time,molecules0,opts);

end