import numpy as np
import librosa
import matplotlib.pyplot as plt
import scipy.io as sio

# constant function
test_vector = np.ones(200)
stft = librosa.stft(test_vector, hop_length=25, n_fft=100, win_length=100, window='hamming')
noise = np.random.random(stft.shape)
stft = stft + noise
istft = librosa.istft(stft, hop_length=25, win_length=100, window='hamming')
plt.plot(istft)
plt.show()
plt.grid()


# sinusoid
winlen = 100
ovp = int(0.75*winlen)
nfft = winlen

fs = 100
t = np.arange(2.*fs)/fs
x = np.cos(2*np.pi*2*t)

X = librosa.stft(x, hop_length=winlen-ovp, n_fft=nfft, win_length=winlen, window='hamming')
X = X + 1*np.random.random(X.shape)

x_out = librosa.istft(X, hop_length=winlen-ovp, win_length=winlen, window='hamming')

plt.plot(x_out)
plt.show()
plt.grid()

plt.pcolormesh(np.abs(X))
plt.show()

# masking of an audio sample
path = '/Users/fatemeh/Dropbox/MyDocuments/My MATLAB Toolboxes/mcft_toolbox_git/demos/audio_examples/2src/oct2/'

s1_name = path + 'src1_samp1_2src.wav'
s2_name = path + 'src2_samp1_2src.wav'
mix_name = path + 'mix_samp1_2src.wav'

s1 = librosa.load(s1_name, sr=44100)[0]
s2 = librosa.load(s2_name, sr=44100)[0]
mix = librosa.load(mix_name, sr=44100)[0]

#plt.plot(s1)
#plt.show()

ibm_thr = 30

winlen = 4096
ovp = int(0.75*winlen)
nfft = winlen

src1_stft = librosa.stft(s1, hop_length=winlen-ovp, n_fft=nfft, win_length=winlen, window='hamming')
src2_stft = librosa.stft(s2, hop_length=winlen-ovp, n_fft=nfft, win_length=winlen, window='hamming')

mix_stft = librosa.stft(mix, hop_length=winlen-ovp, n_fft=nfft, win_length=winlen, window='hamming')

#plt.pcolormesh(np.abs(src1_stft))
#plt.show()

src1_mag = np.abs(src1_stft)
src2_mag = np.abs(src2_stft)

src_snr = 20*np.log10((src1_mag+1e-16)/(src2_mag+1e-16))

#plt.pcolormesh(np.abs(src_snr))
#plt.show()

ibm = src_snr>ibm_thr

src_masked = ibm * mix_stft

src1_est = librosa.istft(src_masked, hop_length=winlen-ovp, win_length=winlen, window='hamming')

plt.plot(4000+np.arange(1200),src1_est[4000:5200])
plt.show()

py_result = src1_est
sio.savemat('/Users/fatemeh/Desktop/py_result.mat',{'py_src1_est':src1_est})

