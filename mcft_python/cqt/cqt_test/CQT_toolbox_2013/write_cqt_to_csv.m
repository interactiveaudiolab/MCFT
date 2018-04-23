function [] = write_cqt_to_csv(path_to_audio, output_filename)
% Simple wrapper for CQT toolbox to write CQT of given audio out to a CSV
% file.
%
% Input types:
%          path_to_audio: string containing the path to the audio file to
%          be used for computing the CQT
%          
%          output_filename: string containing the name of the file to write
%          the CQT out to
[signal,sample_rate] = audioread(path_to_audio);
min_freq = 27.5;
max_freq = 4186;
bins = 24;
Xcq = cqt(signal,bins,sample_rate,min_freq,max_freq);
csvwrite([output_filename,'.dat'], Xcq.c);