clear;
clc;
close all;

Fs = 10000;          % sampling frequency
Fc = 1000;           % cutoff
N  = 100;            % number of taps
order = N-1;
Wn = Fc/(Fs/2);      % normalized cutoff

h = fir1(order, Wn, 'low'); % generating FIR coefficients in float type



h_q = float_to_Q(2, 14, h(:)); % Converting floating-point coefficients into Q(2,14) format. The h_q is a column vector

write_hex16("coef_q214_hex.txt", h_q); % Saving coefficients h_q in 16 bit hex as well as in decimal
write_dec16("coef_q214_dec.txt", h_q);



% Generate 3 sinewaves, quantize, save

freqs = [950 1100 2000];
L = 2000;            % number of samples in each sine
A = 1.0;             % amplitude of sine wave

n = (0:L-1).';
t = n/Fs;
x = cell(1,3);
x_q = cell(1,3);
for i = 1:3
    f = freqs(i);
    x{i} = A*sin(2*pi*f*n/Fs);
    x_q{i} = float_to_Q(2, 14, x{i});

    write_hex16(sprintf("x%d_q214_hex.txt", f), x_q{i});
    write_dec16(sprintf("x%d_q214_dec.txt", f), x_q{i});
end

% MATLAB float outputs (reference)

y_float = cell(1,3);
for i = 1:3
    y_float{i} = filter(h,1,x{i});
end

figure; 
plot(t, y_float{1}); hold on;
plot(t, y_float{2});
plot(t, y_float{3});
grid on;
title('MATLAB floating-point outputs');
xlabel('Time (s)'); ylabel('Amplitude');
legend('950 Hz','1100 Hz','2000 Hz');

% MATLAB fixed-point outputs

y_q = cell(1,3);
y_fx = cell(1,3);

for i = 1:3
    y_q{i}  = fir_fixed_q214(x_q{i}, h_q);
    y_fx{i} = double(y_q{i})/2^14;

    % save reference fixed output too (optional)
    write_dec16(sprintf("y%d_matlab_fixed_dec.txt", freqs(i)), y_q{i});
    write_hex16(sprintf("y%d_matlab_fixed_hex.txt", freqs(i)), y_q{i});
end

figure; 
plot(t, y_fx{1}); hold on;
plot(t, y_fx{2});
plot(t, y_fx{3});
grid on;
title('MATLAB fixed-point outputs (Q2.14 arithmetic)');
xlabel('Time (s)'); ylabel('Amplitude');
legend('950 Hz','1100 Hz','2000 Hz');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ============================================================
% PLOT VERILOG OUTPUTS (3 freqs in SAME figure, exact colors)
% ============================================================

Fs = 10000;
L  = 2000;
freqs = [950 1100 2000];
t = (0:L-1).'/Fs;

% Exact MATLAB default first-3 colors (matches your .fig)
COL = [0    0.447 0.741;
       0.85 0.325 0.098;
       0.929 0.694 0.125];

% ---------- Direct ----------
y950  = read_hex16_q214(sprintf("y%d_direct_hex.txt", freqs(1)), L);
y1100 = read_hex16_q214(sprintf("y%d_direct_hex.txt", freqs(2)), L);
y2000 = read_hex16_q214(sprintf("y%d_direct_hex.txt", freqs(3)), L);

figure;
set(gca,'ColorOrder',COL);   % enforce exact color order
plot(t, y950); hold on;
plot(t, y1100);
plot(t, y2000);
grid on;
title('Verilog output - Direct form (Q2.14 scaled)');
xlabel('Time (s)'); ylabel('Amplitude');
legend('950 Hz','1100 Hz','2000 Hz');

% ---------- Optimized ----------
y950  = read_hex16_q214(sprintf("y%d_opt_hex.txt", freqs(1)), L);
y1100 = read_hex16_q214(sprintf("y%d_opt_hex.txt", freqs(2)), L);
y2000 = read_hex16_q214(sprintf("y%d_opt_hex.txt", freqs(3)), L);

figure;
set(gca,'ColorOrder',COL);
plot(t, y950); hold on;
plot(t, y1100);
plot(t, y2000);
grid on;
title('Verilog output - Optimized symmetric form (Q2.14 scaled)');
xlabel('Time (s)'); ylabel('Amplitude');
legend('950 Hz','1100 Hz','2000 Hz');

% ---------- Genvar ----------
y950  = read_hex16_q214(sprintf("y%d_genvar_hex.txt", freqs(1)), L);
y1100 = read_hex16_q214(sprintf("y%d_genvar_hex.txt", freqs(2)), L);
y2000 = read_hex16_q214(sprintf("y%d_genvar_hex.txt", freqs(3)), L);

figure;
set(gca,'ColorOrder',COL);
plot(t, y950); hold on;
plot(t, y1100);
plot(t, y2000);
grid on;
title('Verilog output - Genvar form (Q2.14 scaled)');
xlabel('Time (s)'); ylabel('Amplitude');
legend('950 Hz','1100 Hz','2000 Hz');

disp('DONE: Plotted 3 freqs per figure with exact MATLAB color order.');

% ============================================================
% Helper: read hex -> int16(two''s comp) -> double/2^14
% (NO readlines; works on old MATLAB)
% ============================================================
function y = read_hex16_q214(fname, L)
    if exist(fname,'file') ~= 2
        error('Missing file: %s', fname);
    end

    fid = fopen(fname,'r');
    if fid < 0
        error('Cannot open file: %s', fname);
    end
    C = textscan(fid, '%s', 'Whitespace', '\n');
    fclose(fid);

    S = C{1};
    if isempty(S)
        error('Empty file: %s', fname);
    end

    % trim + remove empties
    for k = 1:numel(S), S{k} = strtrim(S{k}); end
    S = S(~cellfun('isempty', S));

    M = min(L, numel(S));
    S = S(1:M);

    u = uint16(hex2dec(char(S)));   % hex -> uint16
    s = typecast(u, 'int16');       % reinterpret bits as signed
    y = double(s)/2^14;             % scale Q2.14 -> float

    if M < L
        y = [y; zeros(L-M,1)];
    end
end


% helping functions 

function s = symmetry(h, tol)
if nargin < 2
    tol = 1e-12;
end

h = h(:);
N = numel(h);
s = "Symmetric";

for i = 1:floor(N/2)
    if abs(h(i) - h(N - i + 1)) > tol
        s = "mismatch_at_" + i;
        return
    end
end
end

function q = float_to_Q(a , b, x)
x = double(x);
s1 = pow2(b);
s2 = a+b;
s3 = pow2(s2 - 1);
q = round(s1*x);
q(q > s3 - 1) = s3 - 1;
q(q < - s3) = - s3 ;
if s2 <= 8
    q = int8(q);
elseif s2 <= 16
    q = int16(q);
elseif s2 <= 32
    q = int32(q);
elseif s2 <= 64
    q = int64(q);
else error("Word length %d > 64 not supported with built-in integer types.", s2);
end
end

function write_hex16(fname, x_q)
fid = fopen(fname, "w");
assert(fid > 0, "Can not open %s", fname);
for i = 1:numel(x_q)
    u = typecast(int16(x_q(i)), 'uint16');
    fprintf(fid, "%04X\n", u);    
end
fclose(fid);
end

function write_dec16(fname, x_q)
fid = fopen(fname, "w");
assert(fid > 0, "Can not open %s", fname);
for i = 1: numel(x_q)
    fprintf(fid, "%d\n", x_q(i));    
end
fclose(fid);
end

function y_q = fir_fixed_q214(x_q, h_q)
x_q = x_q(:);
h_q = h_q(:);
N = numel(h_q);
L = numel(x_q);
d = int16(zeros(N,1));
y_q = int16(zeros(L,1));
S = int64(2^14);
for i = 1:L
    d(2:end) = d(1:end - 1);
    d(1) = int16(x_q(i));
    acc = int64(0);
    for k = 1:N 
        acc = acc + int64(d(k))* int64(h_q(k));
    end
    y = acc / S;
    if y > 32767, y = 32767; 
    end
    if y < -32768, y = -32768; 
    end
    y_q(i) = int16(y);
    
end

end

