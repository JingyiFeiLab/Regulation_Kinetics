function dmol_dt = Bestinit_ODE(t,molecules,a_s, b_s, k_init, k_elon, k_elon_prime, b_m, k_on, k_off, kx, kx_s, b_ms, b_e, b_p,k_nuc,b_nuc,n_plasmids,k_off_nuc,Pt)


sRNA = molecules(1); % Free sRNA
mRNA_i = molecules(2);
MiS = molecules(3); % co-transcriptionally bound mRNA-sRNA
mRNA = molecules(4); % Free mRNA, post-transcriptional
MS = molecules(5); % post-transcriptionally bound mRNA
P = molecules(6); % Protein

DNA = n_plasmids;


% c = 0.833725;
% e = 0.649125;
% 
% if b_e == 0
%     d = 0.62885; %g
%     f = 0.609425; %h
% else
%     d = 0.9131; %d
%     f = 0.6165; %f
% end

% c = 0.9309;
% e = 0.9309;
% 
% if b_e == 0
%     d = 0.9309; %g
%     f = 0.9309; %h
% else
%     d = 0.9309; %d
%     f = 0.9309; %f
% end

% c = 0.8128;
% e = 0.81228;
% 
% if b_e == 0
%     d = 0.81228; %g
%     f = 0.81228; %h
% else
%     d = 0.81228; %d
%     f = 0.81228; %f
% end

c = 0;
e = 0;

if b_e == 0
    d = 0.9276; %g
    f = 0; %h
else
    d = Pt; %d
    f = 0; %f
end

   
% ds_dt =  a_s - sRNA*b_s - k_on*sRNA*mRNA - k_nuc*sRNA*mRNA_i + k_off*(MS) + k_off_nuc*MiS + (c*b_e + d*b_ms)*(MiS) + (e*b_e + f*b_ms)*(MS); % sRNA 
%           % ^   ~= a_m           ^%
% 
% dmi_dt = k_init*DNA - k_elon*mRNA_i - k_nuc*sRNA*mRNA_i + k_off_nuc*MiS;        
%           
% dmis_dt =  k_nuc*sRNA*mRNA_i - k_elon_prime * MiS  - (b_e + b_ms)*MiS -k_off_nuc*MiS; % sRNA co-transciptionally bound to mRNA
% 
% dm_dt =  k_elon*mRNA_i - k_on*(sRNA*mRNA) + k_off* MS - b_m*mRNA; % full-length mRNA
% 
% dms_dt = k_elon_prime*MiS + k_on*sRNA*mRNA - b_ms*MS - b_e*MS - k_off*MS; % sRNA bound to fully-transcribed mRNA
% 
% dp_dt = kx*mRNA + kx_s*MS - b_p*P; % Protein

ds_dt =  a_s - sRNA*b_s - k_on*sRNA*mRNA - k_nuc*sRNA*mRNA_i + k_off*(MS) + k_off_nuc*MiS; % sRNA 
          % ^   ~= a_m           ^%

dmi_dt = k_init*DNA - k_elon*mRNA_i - k_nuc*sRNA*mRNA_i + d*(k_off_nuc*MiS);        
          
dmis_dt =  k_nuc*sRNA*mRNA_i - k_off_nuc*MiS; % sRNA co-transciptionally bound to mRNA

dm_dt =  k_elon*mRNA_i - k_on*(sRNA*mRNA) + k_off* MS - b_m*mRNA; % full-length mRNA

dms_dt = k_on*sRNA*mRNA - b_ms*MS - b_e*MS - k_off*MS; % sRNA bound to fully-transcribed mRNA

dp_dt = kx*mRNA + kx_s*MS - b_p*P; % Protein
% 
dmol_dt = [ds_dt; dmi_dt; dmis_dt; dm_dt; dms_dt; dp_dt];

end

%k_init = Initiation rate. Figured out by assuming constant elongation rate
%and breaking transcription into initiation and elongation

%k_elon = elongation rate

% To turn this into kA model, just set k_off to 0. 
