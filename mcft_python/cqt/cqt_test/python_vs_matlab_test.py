import sys
import numpy as np 
import librosa
import matlab.engine


def run_matlab_cqt(path_to_audio):
    # Calls a matlab function in python that returns MATLAB's implementation of the CQT
    eng = matlab.engine.start_engine()
    m_cqt = eng.cqt(path_to_audio)
    return m_cqt


def run_python_cqt(path_to_audio, sample_rate):
    # Calls the python version of the CQT to compare to the MATLAB
    # TODO: Write our own CQT to replace the use of the librosa one
    audio, sr = librosa.core.load(path_to_audio, sr=sample_rate)
    p_cqt = librosa.core.cqt(audio, sr=sample_rate)
    return audio, p_cqt


def reconstruction_error(original, reconstructed):
    # Calculates the difference between two input audio arrays
    length = min(len(original), len(reconstructed))
    return np.sum(abs(reconstructed[:length] - original[:length]))


def main():
    # Load in command line arguments for audio file and its sample rate
    path_to_audio = sys.argv[1]
    sample_rate = int(sys.argv[2])

    # Compute both implementations of the CQT
    m_cqt = np.zeros((84, 80000)) # run_matlab_cqt(path_to_audio)
    audio, p_cqt = run_python_cqt(path_to_audio, sample_rate)

    # Find the total error between the two CQTs
    window = min(p_cqt.shape[1], m_cqt.shape[1])
    diff_matrix = abs(p_cqt[:,:window] - m_cqt[:,:window])
    cqt_error = np.sum(diff_matrix)

    # Reconstruct the audio from the python implementation and compare to the original
    re_audio = librosa.core.icqt(p_cqt, sr=sample_rate)
    re_error = reconstruction_error(audio, re_audio)

    return cqt_error, re_error


if __name__ == '__main__':
    main()