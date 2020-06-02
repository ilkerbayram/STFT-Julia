# Utilities for Short Time Fourier Transform
#
# ilker bayram, ibayram@ieee.org, 2020
#
# this file includes functions for computing
# STFT
# ISTFT : inverse ISTFT
# NormalizeW : Normalizing window so the STFT is self-inverting
# Hann_window : produce Hann window
# IndexSTFT : returns time/frequency indices for display

using FFTW

Hann_window = N -> sin.(π * (1:N)/(N+1)).^2;

function NormalizeW(win, hop::Int)
    # applies a normalization to the window so that
    # the STFT is self-inverting

    norm = win.^2
    for i = 1+hop:hop:length(win)
        norm[1:end-i+1] += win[i:end].^2
        norm[i:end] += win[1:end-i+1].^2
    end
    return win ./ sqrt.(norm)
end

function STFT(x, win, hop::Int)
    # computes the STFT of a sequence 'x'
    # using a window 'win'
    # and a hop size 'hop'

    N0 = length(x) # length of x
    W = length(win) # fft size = window length
    # zero pad to simplify code
    K = Int(ceil((N0-W)/hop)) + 1
    N = (K-1) * hop + W
    pad = N - N0
    if pad > 0
        x = [x;zeros(pad)]
    end

    X = Array{Float64,2}(undef, (W,K)) # initialize the STFT of x
    for i = 1:W
        X[i,:] = x[i:hop:i+(K-1)*hop]
    end
    #X = fft(X .* win, [1])
    ## compute windowed fft's
    #for (i,j) = zip(1:hop:N-W+1, 1:K)
    #    X[:,j] = fft(x[i:i+W-1] .* win)
    #end

    return fft(X .* win, [1])
    #return X
end

function ISTFT(X, win, hop::Int)
    # computes the inverse fft of a 2D complex array 'X'
    # given a window 'win'
    # and a hop size 'hop'

    K = size(X,2) # hop count
    W = size(X,1) # window length = fft size
    N = W + (K-1) * hop # length of the signal
    x = zeros(N)im # initialize x
    Z = ifft(X,[1]) .* win
    for (i,j) = zip(1:hop:N-W+1, 1:K)
        #x[i:i+W-1] += ifft(X[:,j]) .* win
        x[i:i+W-1] += Z[:,j]
    end
    return x
end





function IndexSTFT(size_X, hop, fs)
    # returns time and frequency indices given
    # the size of the STFT, 'size_X',
    # hop size, 'hop',
    # sampling frequency, 'fs'
    W = size_X[1] # number of frequency bins = window length
    K = size_X[2] # hop count
    t = range(W/(2*fs), (W/2 + (K-1)*hop)/fs, length = K)
    w = range(0, fs - 1/W, length = W)
    return t, w
end
