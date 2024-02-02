function varargout = MFGRIDSIMULATOR(varargin)

%      MFGRIDSIMULATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MFGRIDSIMULATOR.M with the given input arguments.
%
%      MFGRIDSIMULATOR('Property','Value',...) creates a new MFGRIDSIMULATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MFGRIDSIMULATOR_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MFGRIDSIMULATOR_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MFGRIDSIMULATOR

% Last Modified by GUIDE v2.5 10-Aug-2015 21:27:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MFGRIDSIMULATOR_OpeningFcn, ...
                   'gui_OutputFcn',  @MFGRIDSIMULATOR_OutputFcn, ...
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

function MFGRIDSIMULATOR_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MFGRIDSIMULATOR (see VARARGIN)

% Choose default command line output for MFGRIDSIMULATOR
handles.output = hObject;
NS = str2double(get(handles.Ns,'string'));
DataPAR = cell(NS,19);
DataG = cell(NS,NS);


 

rnames = cell(1,NS);
for i = 1:NS;
rnames{i} = sprintf('Population Number %d',i) ;
DataPAR{i,1} = sprintf('Label %d',i);
end
%get(handles.uitablePAR)
set(handles.uitablePAR,'RowName',rnames);
set(handles.uitablePAR,'ColumnEditable',logical(ones(1,19))); 
set(handles.uitableg,'ColumnEditable',logical(ones(1,NS)));

set(handles.uitablePAR,'Data',DataPAR) ;







guidata(hObject, handles);

function varargout = MFGRIDSIMULATOR_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function Ns_Callback(hObject, eventdata, handles)
function Ns_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function uitablePAR_CreateFcn(hObject, eventdata, handles)
function T_Callback(hObject, eventdata, handles)
function T_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function dt_Callback(hObject, eventdata, handles)
function dt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






function plot_Callback(hObject, eventdata, handles)
disp('Running Simulation')
T = str2double(get(handles.T,'string')); %Total Time in ms
NS = str2double(get(handles.Ns,'string')); %Number of subpopulations
dt = str2double(get(handles.dt,'string')); %Time step size for Euler Integration 

DataPAR = get(handles.uitablePAR,'Data');
DataG = get(handles.uitableg,'Data');


for i = 1:NS
    v0(i) = DataPAR{i,2}; %initial activity(mV)
    theta(i) = DataPAR{i,3}; %threshold potential(mV)
    strength(i) = DataPAR{i,4}; %strength(mA)
    duration(i) = DataPAR{i,5}; %duration(ms)
    frequency(i) = DataPAR{i,6}; %frequency(Hz)
    delay(i) = DataPAR{i,7}; %delayed time(ms)
    variance(i) = DataPAR{i,8}; %noise variance(ms)
    period(i) = (1000/frequency(i))/dt;
end


%weight matrix
for i = 1:NS
    for j = 1:NS 
       G(i,j)= DataG{i,j};
    end
end


%initial activity
vrec = zeros(T/dt+1,NS);
for i = 1:NS
    vrec(1,i)=v0(i);
end

induction1 = zeros(T/dt+1,NS);
induction2 = zeros(T/dt+1,NS);
%induction
for i = 1:NS
    for j = 1:T/(dt)
        if mod(j,period(i))<=duration(i)/dt
            induction1(j,i)=strength(i);
        end
        induction2(j,i)=strength(i)*sin(2*pi*j/period(i));
        
    end
end





 %for i = 1:T/dt %simulating the network on CTLN model
 %    A=vrec(i,:)*(G')+theta;
 %    A(A<0)=0;
 %    vrec(i+1,:) = vrec(i,:) + dt.*(-vrec(i,:)+A);
 %end
 min1=-5;max1=20;
 min2=-5;max2=15;
 mutiple=2;
 fre = zeros((max1-min1)*mutiple,(max2-min2)*mutiple);
 for variable1 =1:(max1-min1)*mutiple
     variable1/((max1-min1)*mutiple)
     for variable2=1:(max2-min2)*mutiple
       theta(1)=(variable1)/mutiple+min1;
       theta(4)=(variable2)/mutiple+min2;
          vuse = zeros(1,NS);
          noise = zeros(1,NS);
          for i = 1:T/dt %simulating the network on Wilson-cowan model
              for j = 1:NS
                  if i<=delay(j)/dt
                      vuse(j)=0;
                  else
                      vuse(j)=vrec(i-delay(j)/dt,j);
                  end
                  noise(j)=variance(j)*randn(1);
              end
             input=vuse*(G')+theta + induction2(i,:)+noise;
             a1=1.5;
             a2=3;
             diff=-vrec(i,:)+1./(1+exp(-a1*(input-a2)))-1./(1+exp(a1*a2));
             vrec(i+1,:) = vrec(i,:) + dt.*diff/20;
          end


        %caculate frequency
%         %caculate frequency
%         count(NS,1) = 0;
%         location=zeros(2,NS);
%         loc=ones(1,NS);
%         for j = 1:NS
%             count(j)=0;
%             for i = T/(2*dt):T/dt
%                 if vrec(i-1,j)<vrec(i,j)
%                     if vrec(i+1,j)<vrec(i,j)
%                         count(j)=count(j)+1;
%                         if loc(j)<3
%                             location(loc(j),j)=i;
%                             loc(j)=loc(j)+1;
%                         end
%                     end
%                 end
%             end
%         end
%         if max(vrec(T/(2*dt):T/dt,3))<1.01*max(vrec(3*T/(4*dt):T/dt,3))
%             if max(vrec(3*T/(4*dt):T/dt,3))-min(vrec(3*T/(4*dt):T/dt,3))>0.05
%                 fre(variable1,variable2)=max(100000/median(location(2,:)-location(1,:)),0);
%             else
%                 fre(variable1,variable2)=0;
%             end
%         else
%             fre(variable1,variable2)=0;
%         end
      

        Nsamps = 10001;
        fsamp = 100;
        Tsamp = 1/fsamp;
        t = (0:Nsamps-1)*Tsamp;
        
        for i= 1:T*(fsamp/1000)
            v_fre1(i)=vrec((1000/fsamp)*i/dt,3);
        end
        
        % Choose FFT size and calculate spectrum
        Nfft = 1000;
         [Pxx1,f1] = pwelch(v_fre1, gausswin(Nfft), Nfft/2, Nfft,fsamp);
%       [Pxx,f] = pwelch(vtry, gausswin(Nfft), Nfft/2, Nfft,fsamp);
        [m1,p1]=max(Pxx1(10:501));
        % Get frequency estimate (spectral peak)
        [~,loc] = max(Pxx1);
        FREQ_ESTIMATE = f1(loc);
        fre(variable1,variable2)=(p1+9)/10;
     end
 end
 
 

        
       
        
        
 
%         figure;
        yaxis=[min1 max1];
        xaxis=[min2 max2];
%         h=imagesc(xaxis,yaxis,pha);
        
%        

        figure
        
        imagesc(xaxis,yaxis,fre)
        set(gca, 'YDir', 'normal')
        ylabel('Input to D2'); %设置纵坐标轴
        xlabel('Input to STN'); %设置横坐标轴
        title('frequency (Input to arky/p proto is 3/1) '); %标题
        colorbar
        
        
%save data

% for a = 1:10
%     for b = 1:10
%         theta(1)=a/10+1;
%         theta(2)=b/10+5;
%         for i = 1:T/dt %simulating the network on Wilson-cowan model
%               for j = 1:NS
%                   if i<=delay(j)/dt
%                       vuse(j)=0;
%                   else
%                       vuse(j)=vrec(i-delay(j)/dt,j);
%                   end 
%               end
%              input=vuse*(G')+theta + induction(i,:);
%              a1=1.5;
%              a2=3;
%              diff=-vrec(i,:)+1./(1+exp(-a1*(input-a2)))-1./(1+exp(a1*a2));
%              vrec(i+1,:) = vrec(i,:) + dt.*diff/10;
%         end
%         %caculate frequency
%         count = zeros(NS,1);
%         for j = 1:NS
%             for i = T/(2*dt):T/dt
%                 if vrec(i-1,j)<vrec(i,j)
%                     if vrec(i+1,j)<vrec(i,j)
%                         count(j)=count(j)+1;
%                     end
%                 end
%             end
%         end
%         fre=count/(T/2000);
%         fre_matrix(a,b)=median(fre)
%     end
% end
% figure(11)
% h = imagesc(fre_matrix);










guidata(hObject, handles);
















function updatetable_Callback(hObject, eventdata, handles)
 NS = str2double(get(handles.Ns,'string'));
 D = get(handles.uitablePAR,'Data');
 F = get(handles.uitableg,'Data');
 m = size(D);
 v = m(2); 
 m = m(1);
 DataPAR = cell(NS,19); 

DataG = cell(NS,NS);


  if m<NS 
 DataPAR(1:m,:) = D;
 DataG(1:m,1:m) = F;
  end
  



rnames = cell(1,NS);
for i = 1:NS;
rnames{i} = sprintf('Population Number %d',i);
DataPAR{i,1} = sprintf('Label %d',i);
end
set(handles.uitablePAR,'Data',DataPAR);
set(handles.uitableg,'Data',DataG);
set(handles.uitablePAR,'RowName',rnames);
set(handles.uitableg,'ColumnEditable',logical(ones(1,NS)));








function save_Callback(hObject, eventdata, handles)
DataPAR = get(handles.uitablePAR,'Data');
DataG = get(handles.uitableg,'Data');

uisave({'DataPAR','DataG'}) 




function Load_Callback(hObject, eventdata, handles)
uiload 

set(handles.uitablePAR,'Data',DataPAR);
set(handles.uitableg,'Data',DataG);
set(handles.uitablep,'Data',DataP);










function par1_Callback(hObject, eventdata, handles)
function par1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to par1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ind1_Callback(hObject, eventdata, handles)
function ind1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ind1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function minpar1_Callback(hObject, eventdata, handles)
function minpar1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minpar1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function maxpar1_Callback(hObject, eventdata, handles)
function maxpar1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxpar1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Npar1_Callback(hObject, eventdata, handles)
function Npar1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Npar1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function par2_Callback(hObject, eventdata, handles)
function par2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to par2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ind2_Callback(hObject, eventdata, handles)
function ind2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ind2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function minpar2_Callback(hObject, eventdata, handles)
function minpar2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minpar2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function maxpar2_Callback(hObject, eventdata, handles)
function maxpar2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxpar2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Npar2_Callback(hObject, eventdata, handles)
function Npar2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Npar2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function rmf_Callback(hObject, eventdata, handles)
function filename_Callback(hObject, eventdata, handles)
function filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function peakfinders_Callback(hObject, eventdata, handles)
function peakfinders_CreateFcn(hObject, eventdata, handles)
% hObject    handle to peakfinders (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function open_Callback(hObject, eventdata, handles)
% hObject    handle to open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiload 
FD = size(DataPAR); 
NS = FD(1); 
set(handles.Ns,'String',NS); 


rnames = cell(1,NS);
for i = 1:NS;
rnames{i} = sprintf('Population Number %d',i);
%DataPAR{i,1} = sprintf('Label %d',i);
end
set(handles.uitablePAR,'Data',DataPAR);
set(handles.uitableg,'Data',DataG);
set(handles.uitablePAR,'RowName',rnames);
set(handles.uitableg,'ColumnEditable',logical(ones(1,NS)));
set(handles.uitableg,'Data',DataG);

 
if exist('par1')==1
set(handles.par1,'string',par1)
end
if exist('par2')==1
set(handles.par2,'string',par2) 
end
if exist('VREC')==1 
    D = size(VREC); 
    set(handles.Npar1,'string',D(1));
    set(handles.Npar2,'string',D(2)); 
    set(handles.minpar1,'string',VREC{1,1}(1,1));
    set(handles.minpar2,'string',VREC{1,1}(1,2));
    set(handles.maxpar1,'string',VREC{D(1),D(2)}(1,1));
    set(handles.maxpar2,'string',VREC{D(1),D(2)}(1,2));
end






function hetnet_Callback(hObject, eventdata, handles)

function numhet_Callback(hObject, eventdata, handles)
function numhet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numhet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in UpdateHET.
function UpdateHET_Callback(hObject, eventdata, handles)
numhet = str2double(get(handles.numhet,'string'))
DataHET = cell(numhet,3); 
for i = 1:numhet
rnames{1,i} = sprintf('Heterogeneous Parameter %d',i) ;
DataHET{i,1} = sprintf('Label %d',i); 
end
set(handles.uitableHET,'ColumnEditable',logical(ones(1,3)));
set(handles.uitableHET,'Data',DataHET); 
set(handles.uitableHET,'Columnname',{'Heterogeneous Parameter','Population Index','Standard Deviation'}); 


% --- Executes during object deletion, before destroying properties.
function uitablePAR_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to uitablePAR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
