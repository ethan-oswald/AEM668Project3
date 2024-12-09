clear;
clc;
close all;

%%Declare vehicle attributes
vehicle.length = 143;
vehicle.cg_loc = 88.42;
vehicle.wing_area = 1950;
vehicle.wing_chord = 15.3;
vehicle.wing_span = 70;
vehicle.wing_taper = 65;
vehicle.generalized_mass = 184;
vehicle.displacement = 0.32;
vehicle.weight = 288000;
vehicle.cockpit_loc = 63.4;
vehicle.inertia_xx = 1.5*10^5;
vehicle.inertia_yy = 7*10^6;
vehicle.inertia_zz = 7.1*10^6;
vehicle.inertia_xz = -5.27*10^4;
vehicle.frequency = 12.6;
vehicle.slope = -0.027;
vehicle.damping = .01;

%%Declare conditions
conditions.height = 5000;
conditions.Mach = 0.6;
conditions.airspeed = 659;
conditions.AOA = 0;
conditions.lift_coef = 0.34;
conditions.drag_coef = 0.028;
conditions.g = 32.2 * (20902000/(20902000+conditions.height))^2;

%%Declare derivatives
deriv.X_u = -.013;
deriv.X_alpha = 19.45;
deriv.X_alpha_dot = 0;
deriv.X_q = -1.913;
deriv.X_sigma_e = 14.83;
deriv.X_sigma_f = .742;
deriv.X_nu = 0;
deriv.X_nu_dot = 0;
deriv.Z_u = -.1001;
deriv.Z_alpha = -283.3;
deriv.Z_alpha_dot = 0;
deriv.Z_q = 16.55;
deriv.Z_sigma_e = -42.22;
deriv.Z_sigma_f = 2.11;
deriv.Z_nu = -2.812;
deriv.Z_nu_dot = -.0968;
deriv.M_u = .003;
deriv.M_alpha = -3.445;
deriv.M_alpha_dot = -.1035;
deriv.M_q = -3.136;
deriv.M_sigma_e = -5.346;
deriv.M_sigma_f = 3.376;
deriv.M_nu = -0.0663;
deriv.M_nu_dot = -.00372;
deriv.Xi_u = 0;
deriv.Xi_alpha = -1075;
deriv.Xi_alpha_dot = 0;
deriv.Xi_q = -79.44;
deriv.Xi_sigma_e = -923;
deriv.Xi_sigma_f = 89.53;
deriv.Xi_nu = 4.219;
deriv.Xi_nu_dot = -.3502;


%Find Longitudinal Matrix and eigenvalues
A_long = zeros(6,6);
A_long(1,1) = deriv.X_u;
A_long(1,2) = deriv.X_alpha;
A_long(1,3) = deriv.X_q;
A_long(1,4) = -conditions.g;
A_long(1,5) = deriv.X_nu;
A_long(1,6) = deriv.X_nu_dot;
A_long(2,1) = deriv.Z_u/conditions.airspeed;
A_long(2,2) = deriv.Z_alpha/conditions.airspeed;
A_long(2,3) = 1 + deriv.Z_q/conditions.airspeed;
A_long(2,5) = deriv.Z_nu/conditions.airspeed;
A_long(2,6) = deriv.Z_nu_dot/conditions.airspeed;
A_long(3,1) = deriv.M_u;
A_long(3,2) = deriv.M_alpha;
A_long(3,3) = deriv.M_q;
A_long(3,5) = deriv.M_nu;
A_long(3,6) = deriv.M_nu_dot;
A_long(4,3) = 1;
A_long(5,6) = 1;
A_long(6,1) = deriv.Xi_u;
A_long(6,2) = deriv.Xi_alpha;
A_long(6,3) = deriv.Xi_q;
A_long(6,5) = deriv.Xi_nu - vehicle.frequency^2;
A_long(6,6) = deriv.Xi_nu_dot - 2*vehicle.damping*vehicle.frequency;

B_long = zeros(6,2);
B_long(1,1) = deriv.X_sigma_e;
B_long(1,2) = deriv.X_sigma_f;
B_long(2,1) = deriv.Z_sigma_e/conditions.airspeed;
B_long(2,2) = deriv.Z_sigma_f/conditions.airspeed;
B_long(3,1) = deriv.M_sigma_e;
B_long(3,2) = deriv.M_sigma_f;
B_long(6,1) = deriv.Xi_sigma_e;
B_long(6,2) = deriv.Xi_sigma_f;

%Calculate eigenvalues
eigen = eig(A_long);

%Modes

%Short-Period Mode
SP.eigenvalue = eigen(3);
a = real(SP.eigenvalue);
b = imag(SP.eigenvalue);
SP.freq = b;
SP.damp_rat = -a/sqrt(a^2+SP.freq^2);


%Phugoid Mode
Phug.eigenvalue = eigen(5);
a = real(Phug.eigenvalue);
b = imag(Phug.eigenvalue);
Phug.freq = b;
Phug.damp_rat = -a/sqrt(a^2+Phug.freq^2);

%Structural Mode
struct.eigenvalue = eigen(1);
a = real(struct.eigenvalue);
b = imag(struct.eigenvalue);
struct.freq = b;
struct.damp_rat = -a/sqrt(a^2+struct.freq^2);

%Output and Comment
fprintf("The eigenvalues for the Short-Period Mode are:\n")
display(eigen(3))
display(eigen(4))
fprintf("The natural frequency for the Short-Period mode is %5.2f rad/s\n", SP.freq);
fprintf("The damping ratio for the Short-Period mode is %5.2f\n", SP.damp_rat);
fprintf("This frequency value signifies a relatively responsive short-period mode,\n" + ...
    "and the damping ratio suggests that the short-period oscillations are highly damped.\n\n");

fprintf("The eigenvalues for the Phugoid Mode are:\n")
display(eigen(5))
display(eigen(6))
fprintf("The natural frequency for the Phugoid mode is %5.2f rad/s\n", Phug.freq);
fprintf("The damping ratio for the Phugoid mode is %5.2f\n", Phug.damp_rat);
fprintf("This frequency value signifies slow, gradual changes without significant oscillations,\n" + ...
    "and the damping ratio suggests that the oscillations are only lightly damped.\n\n");

fprintf("The eigenvalues for the Structural Mode are:\n")
display(eigen(1))
display(eigen(2))
fprintf("The natural frequency for the Structural mode is %5.2f rad/s\n", struct.freq);
fprintf("The damping ratio for the Structural mode is %5.2f\n", struct.damp_rat);
fprintf("This frequency value signifies signifies high-energy oscillations which likely affect the cockpit and wings,\n" + ...
    "and the damping ratio suggests that the oscillations are only lightly damped,\nwhich could result in sustained oscillations when excited.\n\n");



%Create C and D matrices
C = [0 conditions.airspeed -vehicle.cockpit_loc -conditions.airspeed 0 vehicle.displacement];
D = [0 0];

G_ss = ss(A_long, B_long, C, D);
G = tf(G_ss);

%Display Frequency Responses
figure(1)
bode(G(1))
title("Frequency Response for Elevator Deflection")

figure(2)
bode(G(2))
title("Frequency Response for Fin Deflection")

%Calculate magnitude for all modes for elevator
struct.freqresp1 = freqresp(G(1), struct.freq);
struct.mag1 = 20*log10(sqrt(real(struct.freqresp1)^2 + imag(struct.freqresp1)^2));
Phug.freqresp1 = freqresp(G(1), Phug.freq);
Phug.mag1 = 20*log10(sqrt(real(Phug.freqresp1)^2 + imag(Phug.freqresp1)^2));
SP.freqresp1 = freqresp(G(1), SP.freq);
SP.mag1 = 20*log10(sqrt(real(SP.freqresp1)^2 + imag(SP.freqresp1)^2));

%Print values and comment
fprintf(" As can be seen from the Bode Plot for the Elevator Deflection,\n" + ...
    " the elevator heavily excites the phugoid mode,\n" + ...
    " with a magnitude of %5.2f dB at that frequency.\n" + ...
     "It has a moderately exciting effect on the elastic mode, \n" + ...
     "with a magnitude of %5.2f dB at the corresponding frequency.\n" + ...
     "It does not excite the short period mode, with no peak at that frequency, and a\n" + ...
     "magnitude of %5.2f.\n\n", Phug.mag1, struct.mag1, SP.mag1);

%Calculate magnitude for all modes for fin
struct.freqresp2 = freqresp(G(2), struct.freq);
struct.mag2 = 20*log10(sqrt(real(struct.freqresp2)^2 + imag(struct.freqresp2)^2));
Phug.freqresp2 = freqresp(G(2), Phug.freq);
Phug.mag2 = 20*log10(sqrt(real(Phug.freqresp2)^2 + imag(Phug.freqresp2)^2));
SP.freqresp2 = freqresp(G(2), SP.freq);
SP.mag2 = 20*log10(sqrt(real(SP.freqresp2)^2 + imag(SP.freqresp2)^2));

%Print values and comment
fprintf(" As can be seen from the Bode Plot for the Fin Deflection,\n" + ...
    "the fin also heavily excites the phugoid mode,\n" + ...
    "with a magnitude of %5.2f dB at that frequency.\n" + ...
     "It is less effective at exciting the elastic mode than the \n" + ...
     "elevator, with a magnitude of %5.2f dB at the corresponding frequency.\n" + ...
     "It also does not excite the short-period mode,\n" + ...
     "with a magnitude of %5.2f dB at the corresponding frequency.\n\n", Phug.mag2, struct.mag2, SP.mag2);

%Create Actuators
Actuator.bandwidth = 50;
Actuator.tf = tf(Actuator.bandwidth, [1, Actuator.bandwidth]);

%Create Passive SMC
SMC_passive.damp_rat_n = struct.damp_rat;
SMC_passive.omega = struct.freq;
SMC_passive.damp_rat_id = 30 * SMC_passive.damp_rat_n;
SMC_passive.gain = 5;
SMC_passive.tf = tf([1 2*SMC_passive.damp_rat_n*SMC_passive.omega SMC_passive.omega^2], [1 2*SMC_passive.damp_rat_id*SMC_passive.omega SMC_passive.omega^2]);
SMC_passive.fulltf = SMC_passive.gain * SMC_passive.tf;

%Output root locus of Passive SMC and comment
figure(3)
rlocus(SMC_passive.fulltf)
title("Root Locus for the Passive SMC")
fprintf("The Passive SMC is not designed to influence the structural mode directly,\n" + ...
    "but instead relies on inherent system properties rather than active adjustments, as well as\n" + ...
    "controlling the rest of the modes. If you choose too high of a proportional gain,\n" + ...
    "it can excite the structural modes and loss of stability, among other effects.\n\n");

%Create Active SMC
SMC_active.freq_l = struct.freq*2;
SMC_active.freq_h = struct.freq/2;
SMC_active.tf = tf([1 0], [1 SMC_active.freq_h]) * tf([0 SMC_active.freq_l], [1 SMC_active.freq_l]);
SMC_active.gain = 15;
SMC_active.fulltf = SMC_active.gain * SMC_active.tf;

%Output root locus of Active SMC and comment
figure(4)
rlocus(SMC_active.fulltf)
title("Root Locus for the Active SMC")
fprintf("The Active SMC does influence the structural mode directly by damping\n" + ...
    "structural oscillations and adapting to the dynamics loads. If you choose\n" + ...
    "too high of a proportional gain, like the passive SMC, it can excite the\n" + ...
    "structural modes and cause loss of stability, among other effects.\n\n")

%Create feedback loop for elevator and output bode plot
elev.open_loop = Actuator.tf * G(1);
elev.closed_loop = feedback(elev.open_loop,SMC_passive.fulltf);
figure(5)
bode(elev.closed_loop)
title("Frequency Response for the Elevator With Actuator and Passive SMC")

%Create feedback loop for fin and output bode plot
fin.open_loop = Actuator.tf * G(2);
fin.closed_loop = feedback(fin.open_loop, SMC_active.fulltf);
figure(6)
bode(fin.closed_loop)
title("Frequency Response for the Fin With Actuator and Active SMC")

%Collect Eigenvalues for both feedback loops
elev.eigen = eig(elev.closed_loop);
fin.eigen = eig(fin.closed_loop);


%%Display values for Passive SMC
fprintf("After implementation of the Elevator Actuator and Passive SMC,\n" + ...
    "the Phugoid and Short-Period modes become uncoupled, and the actuator becomes\n" + ...
    "coupled with one of the short-period modes, and the phugoid and short-period\n" + ...
    "mode become coupled.\n\n");

fprintf("For the Passive SMC:\n\n")
fprintf("The Proportional Gain for the passive SMC is %5.2f\n", SMC_passive.gain);
fprintf("The Notch Filter's center frequency is %5.2f rad/s\n", SMC_passive.omega);
fprintf("The Notch Filter's zero damping is %5.2f\n", SMC_passive.damp_rat_n);
fprintf("The Notch Filter's damping ratio increase is %5.2f\n\n", SMC_passive.damp_rat_id);

fprintf("The eigenvalues for the Passive SMC are:\n")
display(elev.eigen(3))
display(elev.eigen(4))
a = real(elev.eigen(3));
b = imag(elev.eigen(3));
freq = b;
damp_rat = -a/sqrt(a^2+freq^2);
fprintf("with a natural frequency of %5.2f rad/s and a damping ratio of %5.2f\n\n", freq, damp_rat);

fprintf("The eigenvalues for the Elastic Mode are:\n")
display(elev.eigen(5))
display(elev.eigen(6))
a = real(elev.eigen(5));
b = imag(elev.eigen(5));
freq = b;
damp_rat = -a/sqrt(a^2+freq^2);
fprintf("with a natural frequency of %5.2f rad/s and a damping ratio of %5.2f\n\n", freq, damp_rat);

struct.freqresp3 = freqresp(elev.closed_loop, freq);
struct.mag3 = 20*log10(sqrt(real(struct.freqresp3)^2 + imag(struct.freqresp3)^2));

fprintf("The eigenvalues for the Coupled Actuator and Short-Period Mode are:\n")
display(elev.eigen(1))
display(elev.eigen(2))
a = real(elev.eigen(1));
b = imag(elev.eigen(1));
freq = b;
damp_rat = -a/sqrt(a^2+freq^2);
fprintf("with a natural frequency of %5.2f rad/s and a damping ratio of %5.2f\n\n", freq, damp_rat);

fprintf("The eigenvalues for the Coupled Phugoid and Short-Period Mode are:\n")
display(elev.eigen(7))
display(elev.eigen(8))
a = real(elev.eigen(7));
b = imag(elev.eigen(7));
freq = b;
damp_rat = -a/sqrt(a^2+freq^2);
fprintf("with a natural frequency of %5.2f rad/s and a damping ratio of %5.2f\n\n", freq, damp_rat);

fprintf("The eigenvalue for the Uncoupled phugoid mode is:\n")
display(elev.eigen(9))
freq = abs(elev.eigen(9));
fprintf("with a natural frequency of %5.2f rad/s\n\n", freq);

fprintf("Because of the nature of the Passive Controller, the peak for the\n" + ...
    "the elastic mode is much lower, with a frequency of only %5.2 dB. The\n" + ...
    " coupled phugoid/short-period mode, as well asthe uncoupled phugoid \n" + ...
    "mode become positive, which should make them unstable, but they are \n" + ...
    "overpowered by the other, stable modes.\n\n", struct.mag3);


%%Display values for Active SMC
fprintf("After implementation of the Fin Actuator and Active SMC,\n" + ...
    "the Phugoid mode and the Active SMC becomes uncoupled, and the actuator becomes\n" + ...
    "coupled with the Active SMC\n\n");

fprintf("For the Active SMC:\n\n")
fprintf("The Proportional Gain for the active SMC is %5.2f\n", SMC_active.gain);
fprintf("The high pass filter corner frequency is %5.2f rad/s\n", SMC_active.freq_h);
fprintf("The low pass filter corner frequency is %5.2f rad/s\n\n", SMC_active.freq_l);

fprintf("The eigenvalues for the coupled Active SMC and Actuator are:\n")
display(fin.eigen(1))
display(fin.eigen(2))
a = real(fin.eigen(1));
b = imag(fin.eigen(1));
freq = b;
damp_rat = -a/sqrt(a^2+freq^2);
fprintf("with a natural frequency of %5.2f rad/s and a damping ratio of %5.2f\n\n", freq, damp_rat);

fprintf("The eigenvalues for the Elastic Mode are:\n")
display(fin.eigen(4))
display(fin.eigen(5))
a = real(fin.eigen(4));
b = imag(fin.eigen(4));
freq = b;
damp_rat = -a/sqrt(a^2+freq^2);
fprintf("with a natural frequency of %5.2f rad/s and a damping ratio of %5.2f\n\n", freq, damp_rat);

fprintf("The eigenvalues for the Short-Period Mode are:\n")
display(fin.eigen(6))
display(fin.eigen(7))
a = real(fin.eigen(6));
b = imag(fin.eigen(6));
freq = b;
damp_rat = -a/sqrt(a^2+freq^2);
fprintf("with a natural frequency of %5.2f rad/s and a damping ratio of %5.2f\n\n", freq, damp_rat);

fprintf("The eigenvalue for the Uncoupled active SMC mode is:\n")
display(fin.eigen(3))
freq = abs(fin.eigen(3));
fprintf("with a natural frequency of %5.2f rad/s\n\n", freq);

fprintf("The eigenvalue for the First Uncoupled Phugoid mode is:\n")
display(fin.eigen(8))
freq = abs(fin.eigen(8));
fprintf("with a natural frequency of %5.3f rad/s\n\n", freq);


fprintf("The eigenvalue for the Second Uncoupled Phugoid mode is:\n")
display(fin.eigen(9))
freq = abs(fin.eigen(9));
fprintf("with a natural frequency of %5.3f rad/s\n\n", freq);

fprintf("Because of the nature of the active SMC, all of the peaks have \n" + ...
    "been damped and the bode plot as become stabilized. Similarly to the passive\n" + ...
    "controller, some of the modes are unstable, but overpowered by the stronger,\n" + ...
    "more stable modes.\n")
