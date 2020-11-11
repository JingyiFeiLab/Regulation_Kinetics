function [dmol_dt,c_list] = coTranODE(t, molecules, a_s, b_s, a_m, a_m_co, b_m, k_on, k_off, kx, kx_s, b_ms, b_e, b_p,k_on_list,sRNA_list,c_means)%,c_list)

%global c_list
sRNA = molecules(1);
mRNA = molecules(2);
mRNA_co = molecules(3);
MS = molecules(4);
P = molecules(5);


[~,si] = min(abs(sRNA_list - sRNA));
[~,ki] = min(abs(k_on_list - k_on));
c = c_means{ki,si};

%c_list = [c_list; c] 
   
ds_dt =  a_s - sRNA*b_s - k_on*sRNA*mRNA + k_off*(MS + mRNA_co) - c*a_m_co; % sRNA 

dm_dt =  a_m*(1-c)  - b_m*mRNA - k_on*sRNA*mRNA + k_off*(MS+mRNA_co); % full-length mRNA

dmco_dt =  a_m_co*c - k_off*mRNA_co; % full-length mRNA

dms_dt = k_on*sRNA*mRNA - b_ms*MS - b_e*MS - k_off*MS; % sRNA post-transcriptionally bound to mRNA

dp_dt = kx*mRNA + kx_s*(MS+mRNA_co) - b_p*P; % Protein

dmol_dt = [ds_dt; dm_dt; dmco_dt; dms_dt; dp_dt];

end

%k_init = Initiation rate. Figured out by assuming constant elongation rate
%and breaking transcription into initiation and elongation

%k_elon = elongation rate

% To turn this into kA model, just set k_off to 0. 
