%filename = 'ProyectoPinkFloyd.wav';
filename = 'test.wav';
[audio,Fs] = audioread(filename);
% Sampling period
T  = 1/Fs;

% Batch size (length)
LBatch = 64;
% [# rows, # columns]
[n, m] = size(audio);
% n is the total number of samples

% Number of batches
N = floor(n/LBatch);        % Integer result of the division

% Frequency Matrix
FreqMatrix = zeros(N,32);

for i = 1:N-1
    % i-th batch
    M = audio(64*i:(64*i+64),1);
    % FFT of the i-th batch
    Y = fft(M, LBatch);
    % Single-sided complex coefficients
    P1 = Y(1:LBatch/2); % Y(1:LBatch/2+1) to have original 33 coefficients
    P1(2:end-1) = 2*P1(2:end-1);
    FreqMatrix(i, :) = P1;
end

% At this point we have the FFT of each batch.
% Now we need to get the Real and imaginary parts of each coefficient,
% and store them in an array of N rows x 64 columns.

RealMatrix = zeros(N, 32);
ImagMatrix = zeros(N, 32);

for i = 1:N
    RealMatrix(i, :) = real(FreqMatrix(i, :));
    ImagMatrix(i, :) = imag(FreqMatrix(i, :));
end

% Now we need to find the maximum value of both Matrices
maxReal = max(RealMatrix(:));
maxImag = max(ImagMatrix(:));
maxTotal = max(maxReal, maxImag);

% Now we need to divide the scale in 16 steps of
% maxTotal/16. Each step has a corresponding 4-bit
% binary code.

step = maxTotal/15;

partition = step:step:maxTotal;
codebook = 0:1:15;

% Real elements are written in the first 32 columns, imaginary 
% elements arw written in the right-most columns.

writeMatrix = zeros(N, 64);
for i = 1:N
    for j = 1:32
       writeR = encoder(RealMatrix(i, j), partition, codebook);
       writeI = encoder(ImagMatrix(i, j), partition, codebook);
       writeMatrix(i, j)    = writeR;
       writeMatrix(i, j+32) = writeI;
    end
end

% CSV File
csvwrite('encodedCSV.csv',writeMatrix);

function code = encoder(value, partition, codebook)
    if value == 0
        code = codebook(1);
    end
    if value ~= 0
        index = 1;
        while partition(index) < value
            index = index + 1;
        end
        code = codebook(index);
    end
end

% % CSV File
% csvwrite('MYCSV.csv',FreqMatrix);
% 
% % Plotting
% Yplot = abs(FreqMatrix(2,:));
% Xplot = Fs*(0:(LBatch/2-1))/LBatch;
% plot(Xplot,Yplot)
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')