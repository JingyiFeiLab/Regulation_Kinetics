function mol_fit = JustsRNA_fit(x0,time,b_s,molecules0)

a_s = x0(1);

[t,mol_fit] = ode45(@(t,molecules) JustsRNA_ODE(time,molecules,a_s, b_s,molecules0),time,molecules0);

end