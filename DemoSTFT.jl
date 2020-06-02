using WAV
using PyPlot
include("STFT.jl")

# get the input wave file and extract a shorter clip
fname = "tanbur.wav"
y, fs = wavread(fname)
y = y[:];

# choose the window and hop-size
win_dur = 50e-3; # in seconds

winsize = Int(win_dur * fs)
win = Hann_window(winsize)
hop = Int(winsize / 4)

# pad with zeros to ensure zero perfect reconstruction error
y = [y; zeros(winsize)]

# apply normalization to the window so the STFT is self-inverting
win2 = NormalizeW(win,hop)

# compute STFT
Y = STFT(y, win2, hop)

# compute inverse STFT
y2 = ISTFT(Y, win2, hop)

# compute reconstruction error
y2 = y2[1:length(y)]
error = y2 - y

# display results
t = range(0,(length(y)-1)/fs, length = length(y)) # time variable, in seconds

# reconstruction error in time
fig, ax = subplots(1,1,figsize = (8,3))
ax.plot(t, real(y2), label = "reconstruction")
ax.plot(t, abs.(error), label = "error magnitude")
ax.set_xlabel("Time (sec)")
ax.legend()
display(fig)
# time-frequency image, demonstrating the use of
# the function IndexSTFT
t, w = IndexSTFT(size(Y), hop, fs)

wmin = 1000 # minimum frequency to show, in Hz
wmax = 3000 # maximum frequency to show, in Hz

indw = (w .<= wmax) .& (w .>= wmin)

fig, ax = subplots(1,1,figsize = (8,3))
ax.imshow(log.(abs.(Y[indw,:])),extent = [t[1],t[end], wmin * 1e-3,wmax * 1e-3], aspect = 0.5)
ax.set_xlabel("Time (sec)")
ax.set_ylabel("Frequency (KHz) ")
display(fig)
