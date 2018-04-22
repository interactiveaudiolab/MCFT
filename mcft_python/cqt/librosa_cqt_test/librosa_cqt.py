import sys
import numpy as np
import matplotlib.pyplot as plt

import librosa
from librosa import display as disp

def main():
    '''
    Computes the CQT of a specified audio input file, plots the CQT if specified, then
    reconstructs the signal using the iCQT and computes the reconstruction error. Command
    line arguments in order are: path to the audio file, audio files sample rate, and a 
    boolean for whether or not to graph the CQT.

    Example command line syntax: python librosaCQT.py trombone_D4_tremolo.wav 22050 1
    '''
    path_to_audio = sys.argv[1]
    sr = int(sys.argv[2])
    display = int(sys.argv[3])
    # Initial load in of audio and computing CQT
    audio, sample_rate = librosa.core.load(path_to_audio, sr=sr)
    cqt = librosa.core.cqt(audio, sr=sample_rate)

    # Reconstructing the signal
    audio_re = librosa.core.icqt(cqt, sr=sample_rate)

    # # Plotting the original and reconstructed signals overlaid
    # plt.plot(range(0,len(audio)), audio, 'g-')
    # plt.plot(range(0,len(audio_re)), audio_re, 'b-', alpha=0.3)
    # plt.show()
    if display:
        # Display functionality, taken from example source in librosa docs
        D = librosa.amplitude_to_db(cqt, ref=np.max)
        disp.specshow(D, y_axis='linear')
        plt.colorbar(format='%+2.0f dB')
        plt.tight_layout()
        plt.show()

    # Calculating the reconstruction error 
    # Hardcoded in to normalize the length of the vectors (for trombone_D4_tremolo.wav)
    err_vec = 20*np.log10(np.linalg.norm(audio_re[:82000]-audio[:82000])/np.linalg.norm(audio_re[:82000]))
    err = np.inner(err_vec, err_vec)
    print "Reconstruction error is: ", err

if __name__ == "__main__":
    main()