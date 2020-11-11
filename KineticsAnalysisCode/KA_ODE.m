function dmol_dt = kA_ODE(t, molecules, a_s, b_s, a_m, b_m, kA, kx, kx_s, b_ms, b_e, b_p,molecules0)

%function dydt = vdp1(t,y)
%VDP1  Evaluate the van der Pol ODEs for mu = 1
%
%   See also ODE113, ODE23, ODE45.

%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2014 The MathWorks, Inc.

%dydt = [y(2); (1-y(1)^2)*y(2)-y(1)];



sRNA1 = molecules(1);
mRNA1 = molecules(2);
MS1 = molecules(3);
P1 = molecules(4);

c = 0.829;

%ds_dt =  a_s - sRNA1*b_s - ka*sRNA1*mRNA1 + kd*MS1 + b_ms*MS1;
ds_dt =  a_s - sRNA1*b_s - kA*sRNA1*mRNA1 + c*(b_ms + b_e)*MS1 ;

dm_dt =  a_m - mRNA1*b_m - kA*sRNA1*mRNA1;

dms_dt = kA*sRNA1*mRNA1 - b_ms*MS1 - b_e*MS1;

dp_dt = kx*mRNA1 + kx_s*MS1 - b_p*P1;

%dmol_dt = [ds_dt + dms_dt; dm_dt + dms_dt; dms_dt; dp_dt];

dmol_dt = [ds_dt; dm_dt; dms_dt; dp_dt];

end





