import sys
import csv
import numpy as np
import matplotlib.pyplot as plt

import librosa
from librosa import display as disp


def load_matlab_cqt(path_to_matlab_cqt):
    # Calls a matlab function in python that returns MATLAB's implementation of the CQT
    with open(path_to_matlab_cqt, 'rb') as f:
        reader = csv.reader(f, delimiter=',')
        rows = []
        for i, row in enumerate(reader):
            row = map(lambda x: x.replace('i','j'), row)
            rows.append(map(complex,row))

    m_cqt = np.zeros((len(rows),len(rows[0])), dtype=complex)
    for i in range(len(rows)):
        m_cqt[i,:] = rows[i]

    return m_cqt


def run_python_cqt(path_to_audio, bins):
    # Calls the python version of the CQT to compare to the MATLAB
    # TODO: Write our own CQT to replace the use of the librosa one
    audio, sr = librosa.core.load(path_to_audio)
    p_cqt = librosa.core.cqt(audio, hop_length=128, sr=sr, n_bins=bins, bins_per_octave=24, fmin=27.5)
    return audio, sr, p_cqt


def reconstruction_error(original, reconstructed):
    # Calculates the difference between two input audio arrays
    length = min(len(original), len(reconstructed))
    return (np.sum(abs(reconstructed[:length]**2 - original[:length]**2)))/(np.sum(abs(original)**2))


def main():
    # Load in command line arguments for audio file and its sample rate
    path_to_audio = sys.argv[1]
    path_to_matlab_cqt = sys.argv[2]
    display = int(sys.argv[3])

    # Compute both implementations of the CQT
    m_cqt = load_matlab_cqt(path_to_matlab_cqt)
    audio, sr, p_cqt = run_python_cqt(path_to_audio, int(m_cqt.shape[0]))

    # # Find the total error between the two CQTs
    # window = min(p_cqt.shape[1], m_cqt.shape[1])
    # diff_matrix = abs(p_cqt[:,:window] - m_cqt[:,:window])
    # cqt_error = np.sum(diff_matrix)
    # # Downsampling is causing CQTs to be different lengths, so comparison doesn't make sense

    if display:
        # Display functionality, taken from example source in librosa docs
        D = librosa.amplitude_to_db(p_cqt, ref=np.max)
        E = librosa.amplitude_to_db(m_cqt, ref=np.max)

        # Plot Python CQT
        plt.subplot(2,1,1)
        plt.title("Python CQT Implementation")
        disp.specshow(D, y_axis='linear')
        plt.colorbar(format='%+2.0f dB')
        plt.tight_layout()
        plt.xticks([0,p_cqt.shape[1]-1], [0,p_cqt.shape[1]-1])

        # Plot MATLab CQT
        plt.subplot(2,1,2)
        plt.title("MATLAB CQT Implementation")
        disp.specshow(E, y_axis='linear')
        plt.colorbar(format='%+2.0f dB')
        plt.tight_layout()
        plt.xticks([0,m_cqt.shape[1]-1], [0,m_cqt.shape[1]-1])
        plt.show()

    # Reconstruct the audio from the python implementation and compare to the original
    re_audio = librosa.core.icqt(p_cqt, hop_length=128, sr=sr, bins_per_octave=24, fmin=27.5)
    librosa.output.write_wav('python_reconstructed_audio.wav',re_audio,sr)
    re_error = reconstruction_error(audio, re_audio)

    print re_error
    return re_error


if __name__ == '__main__':
    main()