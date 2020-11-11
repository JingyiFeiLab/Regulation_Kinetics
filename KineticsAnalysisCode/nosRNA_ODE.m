function dmol_dt = nosRNA_ODE(t,molecules, k_i, kx, b_m, k_elon, b_p, n_plasmids)

%function dydt = vdp1(t,y)
%VDP1  Evaluate the van der Pol ODEs for mu = 1
%
%   See also ODE113, ODE23, ODE45.

%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2014 The MathWorks, Inc.

%dydt = [y(2); (1-y(1)^2)*y(2)-y(1)];

m_i = molecules(1);
m_full = molecules(2);
P = molecules(3);

DNA = n_plasmids;

   
dmi_dt = k_i*DNA - k_elon*m_i;

dm_dt =  k_elon*m_i - b_m*m_full;

dp_dt = kx*m_full - b_p*P;

dmol_dt = [dmi_dt; dm_dt; dp_dt];

end





