function dmol_dt = init_ODE(t, molecules, a_s, b_s, a_m, a_m_co, b_m, k_on, k_off, kx, kx_s, b_ms, b_e, b_p, k_nuc)


sRNA = molecules(1); % Free sRNA
MiS = molecules(2); % co-transcriptionally bound mRNA-sRNA
mRNA = molecules(3); % Free mRNA, post-transcriptional
MS = molecules(4); % post-transcriptionally bound mRNA
P = molecules(5); % Protein

t_elon = 17.7800;
t_ini = 0.25;


mRNA_i = (t_elon/(t_elon+t_ini))*a_m;
   
ds_dt =  a_s - sRNA*b_s - k_on*sRNA*mRNA - k_nuc*sRNA*mRNA_i + k_off*(MS) + a_m_co*MiS; % sRNA 
          % ^   ~= a_m           ^%
          
dmis_dt =  k_nuc*sRNA*mRNA_i - a_m_co*MiS; % sRNA co-transciptionally bound to mRNA

dm_dt =  (a_m -(k_nuc*sRNA*mRNA_i)) + a_m_co*MiS - k_on*(sRNA*mRNA) + k_off* MS - b_m*mRNA; % full-length mRNA

dms_dt = k_on*sRNA*mRNA - b_ms*MS - b_e*MS - k_off*MS; % sRNA bound to fully-transcribed mRNA

dp_dt = kx*mRNA + kx_s*MS - b_p*P; % Protein

dmol_dt = [ds_dt; dmis_dt; dm_dt; dms_dt; dp_dt];

end

%k_init = Initiation rate. Figured out by assuming constant elongation rate
%and breaking transcription into initiation and elongation

%k_elon = elongation rate

% To turn this into kA model, just set k_off to 0. 
