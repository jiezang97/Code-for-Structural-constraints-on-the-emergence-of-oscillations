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

inputall1 = zeros(T/dt+1,NS);
inputall2 = zeros(T/dt+1,NS);
induction1 = zeros(T/dt+1,NS);
induction2 = zeros(T/dt+1,NS);

%induction
for i = 1:NS
    for j = 1:T/(2*dt)
        if mod(j,period(i))<=duration(i)/dt
            induction1(j,i)=strength(i);
        end
        induction2(j,i)=strength(i)*sin(2*pi*j/period(i));
        
    end
end




eig(-eye(NS)+G)
 %for i = 1:T/dt %simulating the network on CTLN model
 %    A=vrec(i,:)*(G')+theta;
 %    A(A<0)=0;
 %    vrec(i+1,:) = vrec(i,:) + dt.*(-vrec(i,:)+A);
 %end
 
          vuse = zeros(1,NS);
          noise = zeros(1,NS);
          for i = 1:T/dt %simulating the network on Wilson-cowan model
              for j = 1:NS
                  if i<=delay(j)/dt
                      vuse(j)=0;
                  else
                      vuse(j)=vrec(i-fix(delay(j)./dt),j);
                  end
                  noise(j)=variance(j)*randn(1);
              end
             input=vuse*(G')+theta + induction2(i,:)+noise;
             a1=1.5;
             a2=3;
             diff=-vrec(i,:)+1./(1+exp(-a1*(input-a2)))-1./(1+exp(a1*a2));
             vrec(i+1,:) = vrec(i,:) + dt.*diff/1;
          end


        %caculate frequency
        
        count(NS,1) = 0;
        location=zeros(2,NS);
        loc=ones(1,NS);
        fre =0;
        for j = 1:NS
            for i = T/(2*dt):T/dt
                if vrec(i-1,j)<vrec(i,j)
                    if vrec(i+1,j)<vrec(i,j)
                        count(j)=count(j)+1;
                        if loc(j)<3
                            location(loc(j),j)=i;
                            loc(j)=loc(j)+1;
                        end
                    end
                end
            end
        end
        
        if max(vrec(T/(2*dt):T/dt,2))<1.1*max(vrec(3*T/(4*dt):T/dt,2))
            if max(vrec(3*T/(4*dt):T/dt,2))-min(vrec(3*T/(4*dt):T/dt,2))>0.01
                fre=100000/(location(2,1)-location(1,1));
            else
                fre=0;
            end
        else
            fre=0;
        end

       fontSize = 18;

        figure(10)
        clf
        ylabel('Activity', 'FontSize',fontSize);
        xlabel('time (ms)', 'FontSize',fontSize);
        plot(0:dt:T,vrec(:,1:NS),'LineWidth',2), hold on  
        ylabel(fre)
        legend(DataPAR{:,1})
%    
      
        
     
        
        
     
       
        
        
%         figure(11)
%         clf
%         title('Input and output of Proto')
%         plot(0:dt:T,vrec(:,3),'LineWidth',2), hold on 
%         plot(0:dt:T,inputall1(:,3),'LineWidth',2), hold on 
%         plot(0:dt:T,inputall2(:,3),'LineWidth',2), hold on 
%         legend('output','intput','F(input)')
       
       
        

 

%save data
  










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
