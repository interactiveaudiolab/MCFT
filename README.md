# MCFT
This toolbox provides MATLAB and Python implementations of the Multi-resolution Common Fate Transform. The MCFT of an audio signal outputs a four-dimensional representation capturing the spectro-temporal modulation patterns of the signal. 

MCFT [increases the separability](https://interactiveaudiolab.github.io/assets/papers/pishdadian_pardo_mcft_journal_2018.pdf) of audio mixtures composed of sources with significant time-frequency overlap and different modulation patterns (e.g. two voices singing in unison, each having differernt vibratto). This allows effective source separation for audio scenes where approaches that work on time-frequency representations (e.g. magnitude spectrograms) fail. 

## MATLAB Code Dependencies
The MCFT toolbox uses [the CQT implementation](http://www.cs.tut.fi/sgn/arg/CQT/) by Schörkhuber et al. to compute the time-frequency representation of the input audio signal. 

## Demos
Audio examples, time-frequency plots, and detailed experimental results are provided in [the demo webpage](https://interactiveaudiolab.github.io/MCFT). 

## Citing
If you are using the MCFT for your research, please cite it using one of the following bibtex citation:

```
@article{pishdadian2018multi,
  title={Multi-resolution common fate transform},
  author={Pishdadian, Fatemeh and Pardo, Bryan},
  journal={IEEE/ACM Transactions on Audio, Speech, and Language Processing},
  volume={27},
  number={2},
  pages={342--354},
  year={2018},
  publisher={IEEE}
}

@inproceedings{pishdadian2017multi,
  title={A multi-resolution approach to common fate-based audio separation},
  author={Pishdadian, Fatemeh and Pardo, Bryan and Liutkus, Antoine},
  booktitle={2017 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP)},
  pages={566--570},
  year={2017},
  organization={IEEE}
}
```

## Lincense

The MCFT toolbox is under an [MIT License](https://opensource.org/licenses/MIT)

MIT License

Copyright (c) 2019 Interactive Audio Lab

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.




