import librosa
import numpy as np


audio,sr = librosa.core.load('chirps.wav')
librosa.output.write_wav('librosa_chirps.wav',audio,sr)

librosa_audio,sr = librosa.core.load('chirps.wav')

error = np.sum(abs(audio - librosa_audio))
print(error)