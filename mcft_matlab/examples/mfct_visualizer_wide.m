function varargout = mfct_visualizer_wide(varargin)
% MFCT_VISUALIZER_WIDE MATLAB code for mfct_visualizer_wide.fig
%      MFCT_VISUALIZER_WIDE, by itself, creates a new MFCT_VISUALIZER_WIDE or raises the existing
%      singleton*.
%
%      H = MFCT_VISUALIZER_WIDE returns the handle to a new MFCT_VISUALIZER_WIDE or the handle to
%      the existing singleton*.
%
%      MFCT_VISUALIZER_WIDE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MFCT_VISUALIZER_WIDE.M with the given input arguments.
%
%      MFCT_VISUALIZER_WIDE('Property','Value',...) creates a new MFCT_VISUALIZER_WIDE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mfct_visualizer_wide_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mfct_visualizer_wide_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mfct_visualizer_wide

% Last Modified by GUIDE v2.5 10-Oct-2016 01:19:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mfct_visualizer_wide_OpeningFcn, ...
                   'gui_OutputFcn',  @mfct_visualizer_wide_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end


% --- Executes just before mfct_visualizer_wide is made visible.
function mfct_visualizer_wide_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mfct_visualizer_wide (see VARARGIN)

% Choose default command line output for mfct_visualizer_wide
handles.output = hObject;

% add the path to the toolbox functions
addpath([cd(cd('..')),'/MCFT']); 

% time-domain signal axes
handles.time_signal=handles.axes1;
set(handles.time_signal,'xtick',[],'ytick',[]);
handles.axes_color=get(handles.time_signal,'Color');

% cqt axes
handles.cqt_magft=handles.axes2;
handles.cqt_magsr=handles.axes3;
set(handles.cqt_magft,'xtick',[],'ytick',[]);
set(handles.cqt_magsr,'xtick',[],'ytick',[]);

% filter axes
handles.h_plot=handles.axes4;
handles.H_plot=handles.axes5;
set(handles.h_plot,'xtick',[],'ytick',[]);
set(handles.H_plot,'xtick',[],'ytick',[]);

% filtered signal axes
handles.filtout_ft_plot=handles.axes6;
handles.filtout_sr_plot=handles.axes7;
set(handles.filtout_ft_plot,'xtick',[],'ytick',[]);
set(handles.filtout_sr_plot,'xtick',[],'ytick',[]);

% load button initial image
[open_icon,~,map] = imread(fullfile(matlabroot,'toolbox','matlab','icons','file_open.png'),'PNG');	% Use the hand icon from MATLAB
open_icon = im2double(open_icon);
map = im2double(map);
map(map==0) = NaN;
open_icon = open_icon.*repmat(map,[1,1,3]);
set(handles.load_audio,'CData',open_icon);

% play button initial image
handles.play_icon=nan(11,11,3);
handles.play_icon(:,1:2,:)=0;
handles.play_icon(2:end-1,3:4,:)=0;
handles.play_icon(3:end-2,5:6,:)=0;
handles.play_icon(4:end-3,7:8,:)=0;
handles.play_icon(5:end-4,9:10,:)=0;
handles.play_icon(6:end-5,11,:)=0;
set(handles.play_button,'CData',handles.play_icon);

% initialize audio player
handles.sig_play=[];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mfct_visualizer_wide wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end


% --- Outputs from this function are returned to the command line.
function varargout = mfct_visualizer_wide_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end


% --- Executes on button press in load_audio.
function load_audio_Callback(hObject, eventdata, handles)
% hObject    handle to load_audio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% stop the audioplayer if it's playing
if ~isempty(handles.sig_play)
    if isplaying(handles.sig_play)
        stop(handles.sig_play);
    end
end

% get the audio file name
[file_name, file_path] = uigetfile({'*.wav'});

if isequal(file_name,0)                                                      % Return if 'cancel'
    return
end

audfile_name = fullfile(file_path, file_name);


% read audio data
[x,fs]=audioread(audfile_name);
handles.sig=x;
handles.fs=fs;

% plot the time-domain audio signal
t=(0:length(x)-1)/fs;
axes(handles.time_signal)
plot(t,x,'b-')
axis tight
title('Time-domain Signal','fontsize',12)
xlabel('Time(s)')

% activate play and cqt buttons 
set(handles.play_button,'Enable','on');
set(handles.cqt,'Enable','on');


% reset elements if the load_audio is used more than once
if ~isempty(get(handles.cqt_magft,'Children'))
    
    % reset the cqt plots
    cla(handles.cqt_magft,'reset')
    set(handles.cqt_magft,'box','on','color',handles.axes_color);
    set(handles.cqt_magft,'xtick',[],'ytick',[]);

    cla(handles.cqt_magsr,'reset')
    set(handles.cqt_magsr,'box','on','color',handles.axes_color);
    set(handles.cqt_magsr,'xtick',[],'ytick',[]);
    
    % reset the filter plots
    cla(handles.h_plot,'reset')
    set(handles.h_plot,'box','on','color',handles.axes_color);
    set(handles.h_plot,'xtick',[],'ytick',[]);

    cla(handles.H_plot,'reset')
    set(handles.H_plot,'box','on','color',handles.axes_color);
    set(handles.H_plot,'xtick',[],'ytick',[]);
    
    
    % reset the filter output plots
    cla(handles.filtout_ft_plot,'reset')
    set(handles.filtout_ft_plot,'box','on','color',handles.axes_color);
    set(handles.filtout_ft_plot,'xtick',[],'ytick',[]);
    
    cla(handles.filtout_sr_plot,'reset')
    set(handles.filtout_sr_plot,'box','on','color',handles.axes_color);
    set(handles.filtout_sr_plot,'xtick',[],'ytick',[]);
    
    % reset the input list and 2D filtering button    
    set(handles.filt_text,'Enable','off');
    set(handles.filt_input,'Enable','off');
    set(handles.filt_2d,'Enable','off');
    
    % reset filter parameter buttons and sliders
    set(handles.direction_text,'Enable','off');
    set(handles.direction_up,'Enable','off');
    set(handles.direction_down,'Enable','off');
    set(handles.scale_text,'Enable','off');
    set(handles.s_text,'Enable','off');
    set(handles.rate_text,'Enable','off');
    set(handles.scale_val,'Enable','off');
    set(handles.rate_val,'Enable','off');
    set(handles.r_text,'Enable','off');

end

% update handles structure
guidata(hObject, handles);

end


% --- Executes on button press in play_button.
function play_button_Callback(hObject, eventdata, handles)
% hObject    handle to play_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% play the audio signal
sig_norm=handles.sig;
sig_fs=handles.fs;

sig_norm=sig_norm/max(abs(sig_norm(:)));
handles.sig_play=audioplayer(sig_norm,sig_fs);

state=get(handles.play_button,'Value');

if state==1
  stop_icon=zeros(11,11,3); 
  set(handles.play_button,'CData',stop_icon);
  play(handles.sig_play);
end

if ~isplaying(handles.sig_play)
    disp('not playing')
end

if state==0
  pause(handles.sig_play);
  set(handles.play_button,'CData',handles.play_icon);
end

% update handles structure
guidata(hObject,handles);

end

% --- Executes on button press in cqt.
function cqt_Callback(hObject, eventdata, handles)
% hObject    handle to cqt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% generate plots
fmin = 27.5*2^(0/12); %C2
fmax = 27.5*2^(87/12); %C7
fres = 24; % bins per octave
fs=handles.fs;
dur_x=length(handles.sig)/fs;

Xcq = cqt(handles.sig, fres, fs, fmin, fmax,'gamma',0,'rasterize','full');
X=Xcq.c;
[LF,LT]=size(X);
handles.LF=LF;
handles.LT=LT;

% time and frequency vectors
fvec=Xcq.fbas; 
tvec=linspace(0,dur_x,LT);

% unwrap and frequency normalize the phase
fmat=repmat(fvec,1,LT);
Xph=unwrap(angle(X),[],2)./fmat;
X=abs(X).*exp(1j*Xph);

% adjust X dimensions (to even values)
LF=LF+mod(LF,2);
LT=LT+mod(LT,2);
X=ifft2(fft2(X,LF,LT));

Xmag=abs(X);
Xmag=Xmag/max(Xmag(:)); % normalize for plotting
Xsr=fftshift(fft2(Xmag));
Xmag_sr=abs(Xsr)/max(abs(Xsr(:)));
svec=linspace(-fres/2+fres/LF,fres/2,LF);
FPS=round(LT/dur_x);
rvec=linspace(-FPS/2+FPS/LT,FPS/2,LT);

handles.sig_cqt=X;
handles.SRF=fres;
handles.FPS=FPS;
handles.fvec=fvec;
handles.tvec=tvec;
handles.svec=svec;
handles.rvec=rvec;

axes(handles.cqt_magft)
imagesc(tvec,fvec,Xmag)
set(gca,'ydir','normal','yscale','log');
ylabel('Frequency(Hz)')
xlabel('Time(s)')
colormap(parula)
colorbar;
title('Magnitude CQT in Frequency-Time Domain','fontsize',12)

axes(handles.cqt_magsr)
mesh(rvec,svec,Xmag_sr)
axis tight
set(gca,'ydir','normal');
ylabel('Scale(cyc/oct)')
xlabel('Rate(Hz)')
colormap(parula)
title('Magnitude CQT in Scale-Rate Domain','fontsize',12)

% activate filt_input slider and 2D-filtering button
set(handles.filt_text,'Enable','on');
set(handles.filt_input,'Enable','on');
set(handles.filt_2d,'Enable','on');

% get the initial status of the filter input
handles.filt_input_status=get(handles.filt_input,'Value');

% update handles structure
guidata(hObject,handles);

end


% --- Executes on selection change in filt_input.
function filt_input_Callback(hObject, eventdata, handles)
% hObject    handle to filt_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns filt_input contents as cell array
%        contents{get(hObject,'Value')} returns selected item from filt_input

% update the status of the filter input
handles.filt_input_status=get(handles.filt_input,'Value');

filt_enable=get(handles.filt_2d,'Enable');
if strcmp(filt_enable,'off')
    set(handles.filt_2d,'Enable','on');
end

% reset the elements below the pop-up menue

    % reset the filter plots
    cla(handles.h_plot,'reset')
    set(handles.h_plot,'box','on','color',handles.axes_color);
    set(handles.h_plot,'xtick',[],'ytick',[]);

    cla(handles.H_plot,'reset')
    set(handles.H_plot,'box','on','color',handles.axes_color);
    set(handles.H_plot,'xtick',[],'ytick',[]);

    % reset the filter output plots
    cla(handles.filtout_ft_plot,'reset')
    set(handles.filtout_ft_plot,'box','on','color',handles.axes_color);
    set(handles.filtout_ft_plot,'xtick',[],'ytick',[]);
    
    cla(handles.filtout_sr_plot,'reset')
    set(handles.filtout_sr_plot,'box','on','color',handles.axes_color);
    set(handles.filtout_sr_plot,'xtick',[],'ytick',[]);
        
    % reset filter parameter buttons and sliders
    set(handles.direction_text,'Enable','off');
    set(handles.direction_up,'Enable','off');
    set(handles.direction_down,'Enable','off');
    set(handles.scale_text,'Enable','off');
    set(handles.s_text,'Enable','off');
    set(handles.rate_text,'Enable','off');
    set(handles.scale_val,'Enable','off');
    set(handles.rate_val,'Enable','off');
    set(handles.r_text,'Enable','off');

% update handles structure
guidata(hObject,handles);



end

% --- Executes during object creation, after setting all properties.
function filt_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filt_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in filt_2d.
function filt_2d_Callback(hObject, eventdata, handles)
% hObject    handle to filt_2d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of filt_2d

comp_state=get(handles.filt_2d,'Value');

if comp_state==1
   set(handles.filt_2d,'String','Computing ...');
   pause(1e-6);
   handles=filt_output(handles);
end

comp_state=get(handles.filt_2d,'Value');

if comp_state==0
   set(handles.filt_2d,'String','2D Filtering');
   
  % activate filter parameter buttons and sliders
  set(handles.direction_text,'Enable','on');
  set(handles.direction_up,'Enable','on');
  set(handles.direction_down,'Enable','on');
  set(handles.scale_text,'Enable','on');
  set(handles.s_text,'Enable','on');
  set(handles.rate_text,'Enable','on');
  set(handles.r_text,'Enable','on');
  set(handles.scale_val,'Enable','on');
  set(handles.rate_val,'Enable','on');
  
end
    
% update handles structure
guidata(hObject,handles);   


end

function handles=filt_output(handles)

  % default filter parameters
  LF=handles.LF;
  LT=handles.LT;
  SRF=handles.SRF;
  FPS=handles.FPS;
  beta=1;
  
  [SV,RV]=filt_default_centers(LF,LT,SRF,FPS);
  
  handles.SV=SV;
  handles.RV=RV;

  H_params=struct('samprate_spec',SRF,'samprate_temp',FPS,'time_const',beta);

  % filter input status
  X=handles.sig_cqt;
  [nfft_s,nfft_r]=size(X);
  input_status=handles.filt_input_status;

  % compute the default filter bank
  if input_status==1
      [h,H]=gen_fbank_hsr(SV,RV,nfft_s,nfft_r,H_params); 
      X_input=abs(X);
  elseif input_status==2
      [h,H]=gen_fbank_hsr(SV,RV,nfft_s,nfft_r,H_params); %,X); %don't need to modulate filters
      X_input=X;
  end    
  
  handles.h=h;  
  handles.H=H;

  % compute the filter-bank output
  Z=cqt_to_mcft(X_input,H);
  Z=Z(:,:,1:size(X,1),1:size(X,2));
  handles.Z=Z;
  
  Zsr=zeros(size(Z));
  for i=1:length(SV)
      for j=1:2*length(RV)
          Zsr(i,j,:,:)=fftshift(fft2(squeeze(Z(i,j,:,:))));
      end
  end
  handles.Zsr=Zsr;
          
  % activate the 2D filtering button
  set(handles.filt_2d,'Value',0)
  
  % discretize the sliders based on computed scale and rate values
  LS=length(SV);
  LR=length(RV);
  
  S_step=1/(LS-1); % default rangle of slider values is [0,1]
  R_step=1/(LR-1);
  
  handles.S_step=S_step;
  handles.R_step=R_step;
  
  set(handles.scale_val,'Min',1)
  set(handles.scale_val,'Max',LS)
  set(handles.scale_val,'Value',1)
  set(handles.scale_val, 'SliderStep',[S_step,S_step]);
  set(handles.s_text,'String',num2str(SV(1)));
  
  set(handles.rate_val,'Min',1)
  set(handles.rate_val,'Max',LR)
  set(handles.rate_val,'Value',1)
  set(handles.rate_val, 'SliderStep',[R_step,R_step]);
  set(handles.r_text,'String',num2str(RV(1)));
  
  % set the initial scale and rate indices for plotting
  handles.Sidx_current=1;  
  if get(handles.direction_up,'Value')
   handles.Ridx_current=LR; 
  elseif get(handles.direction_down,'Value')
   handles.Ridx_current=LR+1;
  end
  
  % plot the output for initial parameters: upward,smin,rmin
  handles=plot_filt_output(handles);
   
  
end


function handles=plot_filt_output(handles)

  S_idx=round(handles.Sidx_current); % make sure the index is integer
  R_idx=round(handles.Ridx_current);
    
  % plot the filters in the tf and sr domains  
  
  h=handles.h; 
  H=handles.H;
    
  h_norm=real(h)/max(real(h(:))); % normalize for plotting
  magnorm_H=abs(H)/max(abs(H(:)));
    
  % color bar range
  c_h=[min(h_norm(:)),max(h_norm(:))];
  c_H=[min(magnorm_H(:)),max(magnorm_H(:))];
  
  h_plot=fftshift(squeeze(h_norm(S_idx,R_idx,:,:)),1);
  mag_H=fftshift(squeeze(magnorm_H(S_idx,R_idx,:,:)));
  
  axes(handles.h_plot)
  imagesc(handles.tvec,handles.fvec,h_plot)
  axis tight
  set(gca,'ydir','normal','yscale','log')
  ylabel('Frequency(Hz)')
  xlabel('Time(s)')
  colormap(parula)
  %caxis(c_h)
  colorbar;
  title('Filter Magnitude in F-T','fontsize',12)
  
  axes(handles.H_plot)
  mesh(handles.rvec,handles.svec,mag_H)
  %set(gca,'ydir','normal');
  axis tight
  ylabel('Scale(cyc/oct)')
  xlabel('Rate(Hz)')
  colormap(parula)
  zlim(c_H)
  title('Filter Magnitude in S-R','fontsize',12)
  
  % plot filter outputs in the tf and sr domains  

  Z=handles.Z; 
  Zsr=handles.Zsr;
    
  magnorm_Z=abs(Z)/max(abs(Z(:)));  % normalize for plotting
  magnorm_Zsr=abs(Zsr)/max(abs(Zsr(:)));
    
  % color bar range
  c_ft=[min(magnorm_Z(:)),max(magnorm_Z(:))];
  c_sr=[min(magnorm_Zsr(:)),max(magnorm_Zsr(:))];
    
  magZ_ft=squeeze(magnorm_Z(S_idx,R_idx,:,:));
  magZ_sr=squeeze(magnorm_Zsr(S_idx,R_idx,:,:));
  
  axes(handles.filtout_ft_plot)
  imagesc(handles.tvec,handles.fvec,magZ_ft)
  axis tight
  set(gca,'ydir','normal','yscale','log')
  ylabel('Frequency(Hz)')
  xlabel('Time(s)')
  colormap(parula)
  colorbar;
  caxis(c_ft)
  title('Filter Ouput in Frequency-Time Domain','fontsize',12)
  
  axes(handles.filtout_sr_plot)
  mesh(handles.rvec,handles.svec,magZ_sr)
  %set(gca,'ydir','normal');
  axis tight
  ylabel('Scale(cyc/oct)')
  xlabel('Rate(Hz)')
  colormap(parula)
  zlim(c_sr)
  title('Filter Ouput in Scale-Rate Domain','fontsize',12)
  
end



% --- Executes on slider movement.
function scale_val_Callback(hObject, eventdata, handles)
% hObject    handle to scale_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% read and show the current position of the slider
SV=handles.SV;
S_idx=round(get(handles.scale_val,'Value'));
set(handles.s_text,'String',num2str(SV(S_idx)));

% update the current scale index for plotting
handles.Sidx_current=S_idx;

% update the filter output plots
handles=plot_filt_output(handles);

% update handles structure
guidata(hObject,handles);   


end

% --- Executes during object creation, after setting all properties.
function scale_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scale_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

% --- Executes on slider movement.
function rate_val_Callback(hObject, eventdata, handles)
% hObject    handle to rate_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% read and show the current position of the slider
RV=handles.RV;
%RV
R_idx=round(get(handles.rate_val,'Value'));
%R_ratio=round(RV(R_idx)/RV(1));
%R_string=[num2str(R_ratio),'* Rmin (',num2str(RV(1)),')'];
set(handles.r_text,'String',num2str(RV(R_idx)));

% update the current rate index for plotting
if get(handles.direction_up,'Value')
  handles.Ridx_current=length(RV)-R_idx+1;
elseif get(handles.direction_down,'Value')
  handles.Ridx_current=length(RV)+R_idx;
end

% update the filter output plots
handles=plot_filt_output(handles);

% update handles structure
guidata(hObject,handles);  

end

% --- Executes during object creation, after setting all properties.
function rate_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rate_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end


% --- Executes when selected object is changed in direction_group.
function direction_group_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in direction_group 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get the state of radio buttons 
RV=handles.RV;
if get(handles.direction_up,'Value')
   R_idx=get(handles.rate_val,'Value');
   handles.Ridx_current=length(RV)-R_idx+1; 
elseif get(handles.direction_down,'Value')
   R_idx=get(handles.rate_val,'Value');
   handles.Ridx_current=length(RV)+R_idx;
end

% update the filter output plots
handles=plot_filt_output(handles);

% update handles structure
guidata(hObject,handles);   


end
