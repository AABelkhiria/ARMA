function varargout = ARMA(varargin)
% ARMA MATLAB code for ARMA.fig
%      ARMA, by itself, creates a new ARMA or raises the existing
%      singleton*.
%
%      H = ARMA returns the handle to a new ARMA or the handle to
%      the existing singleton*.
%
%      ARMA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ARMA.M with the given input arguments.
%
%      ARMA('Property','Value',...) creates a new ARMA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ARMA_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ARMA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ARMA

% Last Modified by GUIDE v2.5 22-Aug-2017 11:44:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ARMA_OpeningFcn, ...
                   'gui_OutputFcn',  @ARMA_OutputFcn, ...
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


% --- Executes just before ARMA is made visible.
function ARMA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ARMA (see VARARGIN)
% clear all
try
    handles.Link([0 0 0 0]);
catch
    fprintf('Installing RVC\n')
    run('rvctools\startup_rvc.m');
end

%                       d        a        alpha
handles.L(1) = Link([0  0.24     0.01    -pi/2 ]) ;
handles.L(2) = Link([0  0        0.219    0    ]) ;
handles.L(3) = Link([0  0        0.021   -pi/2 ]) ;   
handles.L(4) = Link([0  0.14     0        pi/2 ]) ;   
handles.L(5) = Link([0  0        0       -pi/2 ]) ;
handles.L(6) = Link([0  0.11     0        0    ]) ;

handles.ARMA = SerialLink(handles.L, 'name', 'ARMA');

% handles.ARMA.tool = transl([0, 0, 0.2]);
% handles.ARMA.base = transl([0.01, 0, 0]);

global q xo yo zo  T t1 t2 t3 t4 t5 t6 Ja T01 T12;
syms t1 t2 t3 t4 t5 t6

alpha1 = -pi/2  ;
alpha2 =  0     ;
alpha3 = -pi/2  ;
alpha4 =  pi/2  ;
alpha5 = -pi/2  ;
alpha6 =  0     ;

a1 = 0.01   ; %0.01
a2 = 0.219  ; %0
a3 = 0.021  ; %0.19
a4 = 0      ; %0.021
a5 = 0      ;
a6 = 0      ;

d1 = 0.24   ;
d2 = 0      ;
d3 = 0      ;
d4 = 0.14   ;
d5 = 0      ;
d6 = 0.11   ;

% T01 = [cos(t1) -sin(t1)*cos(alpha1) sin(alpha1)*sin(t1) a1*cos(t1);sin(t1) cos(t1)*cos(alpha1) -cos(alpha1)*sin(t1) a1*sin(t1);0 sin(alpha1) cos(alpha1) d1 ;0 0 0 1];
% T12 = [cos(t2) -sin(t2)*cos(alpha2) sin(alpha2)*sin(t2) a2*cos(t2);sin(t2) cos(t2)*cos(alpha2) -cos(alpha2)*sin(t2) a2*sin(t2);0 sin(alpha2) cos(alpha2) d2 ;0 0 0 1];
% T23 = [cos(t3) -sin(t3)*cos(alpha3) sin(alpha3)*sin(t3) a3*cos(t3);sin(t3) cos(t3)*cos(alpha3) -cos(alpha3)*sin(t3) a3*sin(t3);0 sin(alpha3) cos(alpha3) d3 ;0 0 0 1];
% T34 = [cos(t4) -sin(t4)*cos(alpha4) sin(alpha4)*sin(t4) a4*cos(t4);sin(t4) cos(t4)*cos(alpha4) -cos(alpha4)*sin(t4) a4*sin(t4);0 sin(alpha4) cos(alpha4) d4 ;0 0 0 1];
% T45 = [cos(t5) -sin(t5)*cos(alpha5) sin(alpha5)*sin(t5) a5*cos(t5);sin(t5) cos(t5)*cos(alpha5) -cos(alpha5)*sin(t5) a5*sin(t5);0 sin(alpha5) cos(alpha5) d5 ;0 0 0 1];
% T56 = [cos(t6) -sin(t6)*cos(alpha6) sin(alpha6)*sin(t6) a6*cos(t6);sin(t6) cos(t6)*cos(alpha6) -cos(alpha6)*sin(t6) a6*sin(t6);0 sin(alpha6) cos(alpha6) d6 ;0 0 0 1];
    
T01 = [cos(t1) 0 -sin(t1) a1*cos(t1); sin(t1) 0 cos(t1) a1*sin(t1); 0 -1 0 d1; 0 0 0 1];
T23 = [cos(t3) 0 -sin(t3) a3*cos(t3); sin(t3) 0 cos(t3) a3*sin(t3); 0 -1 0 0; 0 0 0 1];
T12 = [cos(t2) -sin(t2) 0 a2*cos(t2); sin(t2) cos(t2) 0 a2*sin(t2); 0 0 1 0; 0 0 0 1];
T34 = [cos(t4) 0 sin(t4) 0; sin(t4) 0 -cos(t4) 0; 0 1 0 d4; 0 0 0 1];
T45 = [cos(t5) 0 -sin(t5) 0; sin(t5) 0 cos(t5) 0; 0 -1 0 0; 0 0 0 1];
T56 = [cos(t6) -sin(t6) 0 0; sin(t6) cos(t6) 0 0; 0 0 1 d6; 0 0 0 1];

T= T01*T12*T23*T34*T45*T56;

vect = [T(1,4),T(2,4),T(3,4)];
tet= [t1,t2,t3,t4,t5,t6];

T(t1,t2,t3,t4,t5,t6) = T;

K=jacobian(vect,tet);


% J(1,1) = 0.2*sin(t5)*(sin(t4)*(cos(t3)*(0*cos(t1)*cos(t2) - 0*sin(t1)*sin(t2)) - 1.0*sin(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1))) + cos(t4)*(cos(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1)) + sin(t3)*(0*cos(t1)*cos(t2) - 0*sin(t1)*sin(t2)))) - 0.01*sin(t1) - 0.021*sin(t4)*(cos(t3)*(0*cos(t1)*cos(t2) - 0*sin(t1)*sin(t2)) - 1.0*sin(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1))) + 0*cos(t1)*sin(t2) - 0.219*sin(t1)*sin(t2) - 0*sin(t5)*(0*cos(t1)*sin(t2) - 0*cos(t4)*(cos(t3)*(0*cos(t1)*cos(t2) - 0*sin(t1)*sin(t2)) - 1.0*sin(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1))) - 1.0*sin(t1)*sin(t2) + 0*sin(t4)*(cos(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1)) + sin(t3)*(0*cos(t1)*cos(t2) - 0*sin(t1)*sin(t2))) + sin(t3)*(0*cos(t1)*cos(t2) - 0*sin(t1)*sin(t2))) - 0.021*cos(t4)*(cos(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1)) + sin(t3)*(0*cos(t1)*cos(t2) - 0*sin(t1)*sin(t2))) - 0*sin(t4)*(cos(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1)) + sin(t3)*(0*cos(t1)*cos(t2) - 0*sin(t1)*sin(t2))) - 0.19*cos(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1)) + 0.029*sin(t3)*(0*cos(t1)*cos(t2) - 0*sin(t1)*sin(t2));
% J(1,2) = 0.2*sin(t5)*(sin(t4)*(cos(t3)*(0*cos(t1)*cos(t2) - 0*sin(t1)*sin(t2)) - 1.0*sin(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1))) + cos(t4)*(cos(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1)) + sin(t3)*(0*cos(t1)*cos(t2) - 0*sin(t1)*sin(t2)))) + 0.219*cos(t1)*cos(t2) - 0.021*sin(t4)*(cos(t3)*(0*cos(t1)*cos(t2) - 0*sin(t1)*sin(t2)) - 1.0*sin(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1))) + 0*cos(t2)*sin(t1) - 0.021*cos(t4)*(cos(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1)) + sin(t3)*(0*cos(t1)*cos(t2) - 0*sin(t1)*sin(t2))) - 0*sin(t4)*(cos(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1)) + sin(t3)*(0*cos(t1)*cos(t2) - 0*sin(t1)*sin(t2))) - 0*sin(t5)*(cos(t1)*cos(t2) - 0*cos(t4)*(cos(t3)*(0*cos(t1)*cos(t2) - 0*sin(t1)*sin(t2)) - 1.0*sin(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1))) + 0*cos(t2)*sin(t1) + 0*sin(t4)*(cos(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1)) + sin(t3)*(0*cos(t1)*cos(t2) - 0*sin(t1)*sin(t2))) + sin(t3)*(0*cos(t1)*cos(t2) - 0*sin(t1)*sin(t2))) - 0.19*cos(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1)) + 0.029*sin(t3)*(0*cos(t1)*cos(t2) - 0*sin(t1)*sin(t2));
% J(1,3) = 0.029*cos(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) + 0.021*sin(t4)*(sin(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) - 1.0*cos(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))) - 0*sin(t5)*(0*cos(t4)*(sin(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) - 1.0*cos(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))) + cos(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) + 0*sin(t4)*(cos(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) + sin(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2)))) - 0.021*cos(t4)*(cos(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) + sin(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))) - 0.19*sin(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2)) - 0*sin(t4)*(cos(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) + sin(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))) - 0.2*sin(t5)*(1.0*sin(t4)*(sin(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) - 1.0*cos(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))) - cos(t4)*(cos(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) + sin(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))));
% J(1,4) = 0.021*sin(t4)*(sin(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) - 1.0*cos(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))) - 0*cos(t4)*(sin(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) - 1.0*cos(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))) - 0*sin(t5)*(0*cos(t4)*(sin(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) - 1.0*cos(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))) + 0*sin(t4)*(cos(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) + sin(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2)))) - 0.021*cos(t4)*(cos(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) + sin(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))) - 0.2*sin(t5)*(1.0*sin(t4)*(sin(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) - 1.0*cos(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))) - cos(t4)*(cos(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) + sin(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))));
% J(1,5) = 0.2*cos(t5)*(cos(t4)*(sin(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) - 1.0*cos(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))) + sin(t4)*(cos(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) + sin(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2)))) - 0*cos(t5)*(cos(t1)*sin(t2) + 0*sin(t1)*sin(t2) + 0*sin(t4)*(sin(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) - 1.0*cos(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))) + sin(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) - 0*cos(t4)*(cos(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) + sin(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))));
% J(1,6) = 0;
% J(2,1) = 0.01*cos(t1) + 0.219*cos(t1)*sin(t2) - 0*sin(t5)*(cos(t1)*sin(t2) - 0*cos(t1) + 0*sin(t1)*sin(t2) + 0*sin(t4)*(sin(t3)*(cos(t1) + 0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) - 1.0*cos(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))) + sin(t3)*(cos(t1) + 0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) - 0*cos(t4)*(cos(t3)*(cos(t1) + 0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) + sin(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2)))) + 0*sin(t1)*sin(t2) - 0.021*cos(t4)*(sin(t3)*(cos(t1) + 0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) - 1.0*cos(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))) - 0*sin(t4)*(sin(t3)*(cos(t1) + 0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) - 1.0*cos(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))) + 0.029*sin(t3)*(cos(t1) + 0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) + 0.2*sin(t5)*(cos(t4)*(sin(t3)*(cos(t1) + 0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) - 1.0*cos(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))) + sin(t4)*(cos(t3)*(cos(t1) + 0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) + sin(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2)))) - 0.021*sin(t4)*(cos(t3)*(cos(t1) + 0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) + sin(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))) + 0.19*cos(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2));
% J(2,2) = 0.219*cos(t2)*sin(t1) - 0*cos(t1)*cos(t2) - 0.021*cos(t4)*(sin(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) - 1.0*cos(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))) - 0*sin(t4)*(sin(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) - 1.0*cos(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))) + 0.029*sin(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) + 0.19*cos(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2)) + 0.2*sin(t5)*(cos(t4)*(sin(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) - 1.0*cos(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))) + sin(t4)*(cos(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) + sin(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2)))) - 0.021*sin(t4)*(cos(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) + sin(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))) - 0*sin(t5)*(cos(t2)*sin(t1) - 0*cos(t1)*cos(t2) + 0*sin(t4)*(sin(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) - 1.0*cos(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2))) + sin(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) - 0*cos(t4)*(cos(t3)*(0*cos(t1)*sin(t2) + 0*cos(t2)*sin(t1)) + sin(t3)*(cos(t1)*cos(t2) - 1.0*sin(t1)*sin(t2)))); 
% J(2,3) = 0.021*sin(t4)*(sin(t3)*(sin(t1) - 0*cos(t1)*cos(t2) + 0*sin(t1)*sin(t2)) - 1.0*cos(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1))) - 0.2*sin(t5)*(1.0*sin(t4)*(sin(t3)*(sin(t1) - 0*cos(t1)*cos(t2) + 0*sin(t1)*sin(t2)) - 1.0*cos(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1))) - cos(t4)*(cos(t3)*(sin(t1) - 0*cos(t1)*cos(t2) + 0*sin(t1)*sin(t2)) + sin(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1)))) - 0*sin(t5)*(0*sin(t4)*(cos(t3)*(sin(t1) - 0*cos(t1)*cos(t2) + 0*sin(t1)*sin(t2)) + sin(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1))) + cos(t3)*(sin(t1) - 0*cos(t1)*cos(t2) + 0*sin(t1)*sin(t2)) + 0*cos(t4)*(sin(t3)*(sin(t1) - 0*cos(t1)*cos(t2) + 0*sin(t1)*sin(t2)) - 1.0*cos(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1)))) - 0.021*cos(t4)*(cos(t3)*(sin(t1) - 0*cos(t1)*cos(t2) + 0*sin(t1)*sin(t2)) + sin(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1))) - 0*sin(t4)*(cos(t3)*(sin(t1) - 0*cos(t1)*cos(t2) + 0*sin(t1)*sin(t2)) + sin(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1))) + 0.029*cos(t3)*(sin(t1) - 0*cos(t1)*cos(t2) + 0*sin(t1)*sin(t2)) - 0.19*sin(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1));
% J(2,4) = 0.021*sin(t4)*(sin(t3)*(sin(t1) - 0*cos(t1)*cos(t2) + 0*sin(t1)*sin(t2)) - 1.0*cos(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1))) - 0.2*sin(t5)*(1.0*sin(t4)*(sin(t3)*(sin(t1) - 0*cos(t1)*cos(t2) + 0*sin(t1)*sin(t2)) - 1.0*cos(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1))) - cos(t4)*(cos(t3)*(sin(t1) - 0*cos(t1)*cos(t2) + 0*sin(t1)*sin(t2)) + sin(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1)))) - 0.021*cos(t4)*(cos(t3)*(sin(t1) - 0*cos(t1)*cos(t2) + 0*sin(t1)*sin(t2)) + sin(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1))) - 0*sin(t5)*(0*sin(t4)*(cos(t3)*(sin(t1) - 0*cos(t1)*cos(t2) + 0*sin(t1)*sin(t2)) + sin(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1))) + 0*cos(t4)*(sin(t3)*(sin(t1) - 0*cos(t1)*cos(t2) + 0*sin(t1)*sin(t2)) - 1.0*cos(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1)))) - 0*cos(t4)*(sin(t3)*(sin(t1) - 0*cos(t1)*cos(t2) + 0*sin(t1)*sin(t2)) - 1.0*cos(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1)));
% J(2,5) = 0*cos(t5)*(0*sin(t1) - 0*sin(t4)*(sin(t3)*(sin(t1) - 0*cos(t1)*cos(t2) + 0*sin(t1)*sin(t2)) - 1.0*cos(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1))) + 0*cos(t1)*sin(t2) - 1.0*sin(t1)*sin(t2) + 0*cos(t4)*(cos(t3)*(sin(t1) - 0*cos(t1)*cos(t2) + 0*sin(t1)*sin(t2)) + sin(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1))) - 1.0*sin(t3)*(sin(t1) - 0*cos(t1)*cos(t2) + 0*sin(t1)*sin(t2))) + 0.2*cos(t5)*(sin(t4)*(cos(t3)*(sin(t1) - 0*cos(t1)*cos(t2) + 0*sin(t1)*sin(t2)) + sin(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1))) + cos(t4)*(sin(t3)*(sin(t1) - 0*cos(t1)*cos(t2) + 0*sin(t1)*sin(t2)) - 1.0*cos(t3)*(cos(t1)*sin(t2) + cos(t2)*sin(t1))));
% J(2,6) = 0;
% J(3,1) = 0;
% J(3,2) = 0;
% J(3,3) = 0.021*cos(t3)*cos(t4) - 0.029*cos(t3) + 0*cos(t3)*sin(t4) - 0.021*sin(t3)*sin(t4) + 0*sin(t5)*(cos(t3) + 0*cos(t3)*sin(t4) + 0*cos(t4)*sin(t3)) - 0.2*sin(t5)*(cos(t3)*cos(t4) - 1.0*sin(t3)*sin(t4));
% J(3,4) = 0.021*cos(t3)*cos(t4) + 0*cos(t4)*sin(t3) - 0.021*sin(t3)*sin(t4) + 0*sin(t5)*(0*cos(t3)*sin(t4) + 0*cos(t4)*sin(t3)) - 0.2*sin(t5)*(cos(t3)*cos(t4) - 1.0*sin(t3)*sin(t4));
% J(3,5) = 0*cos(t5)*(sin(t3) - 0*cos(t3)*cos(t4) + 0*sin(t3)*sin(t4) - 0) - 0.2*cos(t5)*(cos(t3)*sin(t4) + cos(t4)*sin(t3));
% J(3,6) = 0;

Ja = K ;
Ja(t1,t2,t3,t4,t5,t6) = K ;

xo = 0.2602 ;
yo = 0      ;
zo = 0.4798 ;

%handles.type = 'nooffset'; %'nooffset'   'puma'
% handles.ARMA.ikineType = handles.type ;

q = [0 -1.57 0 0 0 0];
handles.ARMA.plot(q);

set(handles.text46,'String','Ready to move');

% Choose default command line output for ARMA
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ARMA wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ARMA_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global q T
syms t1 t2 t3 t4 t5 t6

val = get(hObject,'Value');

q(1) = val;

val = val * 1000;
val = floor(val);
val = val / 1000;
txt = sprintf('%.3f', val);
set(handles.q1,'String',txt);

handles.ARMA.plot(q);
guidata(hObject, handles);

H = eval(T(q(1), q(2), q(3), q(4), q(5), q(6) ));

txt1 = sprintf('%.4f', H(1,4));
txt2 = sprintf('%.4f', H(2,4));
txt3 = sprintf('%.4f', H(3,4));
set(handles.text36,'String',txt1);
set(handles.text37,'String',txt2);
set(handles.text38,'String',txt3);

set(handles.slider7, 'Value', H(1,4));
set(handles.slider8, 'Value', H(2,4));
set(handles.slider9, 'Value', H(3,4));

set(handles.edit1,'String',txt1);
set(handles.edit2,'String',txt2);
set(handles.edit3,'String',txt3);

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles.q(1) = get(handles.q1,'Value');
global q T
syms t1 t2 t3 t4 t5 t6

val = get(hObject,'Value');

q(2) = val;

val = val * 1000;
val = floor(val);
val = val / 1000;
txt = sprintf('%.3f', val);
set(handles.text9,'String',txt);

handles.ARMA.plot(q);
guidata(hObject, handles);

H = eval(T(q(1), q(2), q(3), q(4), q(5), q(6) ));

txt1 = sprintf('%.4f', H(1,4));
txt2 = sprintf('%.4f', H(2,4));
txt3 = sprintf('%.4f', H(3,4));
set(handles.text36,'String',txt1);
set(handles.text37,'String',txt2);
set(handles.text38,'String',txt3);

set(handles.slider7, 'Value', H(1,4));
set(handles.slider8, 'Value', H(2,4));
set(handles.slider9, 'Value', H(3,4));

set(handles.edit1,'String',txt1);
set(handles.edit2,'String',txt2);
set(handles.edit3,'String',txt3);

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global q T
syms t1 t2 t3 t4 t5 t6

val = get(hObject,'Value');

q(3) = val;

val = val * 1000;
val = floor(val);
val = val / 1000;
txt = sprintf('%.3f', val);
set(handles.text10,'String',txt);

handles.ARMA.plot(q);
guidata(hObject, handles);

H = eval(T(q(1), q(2), q(3), q(4), q(5), q(6) ));

txt1 = sprintf('%.4f', H(1,4));
txt2 = sprintf('%.4f', H(2,4));
txt3 = sprintf('%.4f', H(3,4));
set(handles.text36,'String',txt1);
set(handles.text37,'String',txt2);
set(handles.text38,'String',txt3);

set(handles.slider7, 'Value', H(1,4));
set(handles.slider8, 'Value', H(2,4));
set(handles.slider9, 'Value', H(3,4));

set(handles.edit1,'String',txt1);
set(handles.edit2,'String',txt2);
set(handles.edit3,'String',txt3);

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global q T
syms t1 t2 t3 t4 t5 t6

val = get(hObject,'Value');

q(4) = val;

val = val * 1000;
val = floor(val);
val = val / 1000;
txt = sprintf('%.3f', val);
set(handles.text11,'String',txt);

handles.ARMA.plot(q);
guidata(hObject, handles);

H = eval(T(q(1), q(2), q(3), q(4), q(5), q(6) ));

txt1 = sprintf('%.4f', H(1,4));
txt2 = sprintf('%.4f', H(2,4));
txt3 = sprintf('%.4f', H(3,4));
set(handles.text36,'String',txt1);
set(handles.text37,'String',txt2);
set(handles.text38,'String',txt3);

set(handles.slider7, 'Value', H(1,4));
set(handles.slider8, 'Value', H(2,4));
set(handles.slider9, 'Value', H(3,4));

set(handles.edit1,'String',txt1);
set(handles.edit2,'String',txt2);
set(handles.edit3,'String',txt3);

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global q T
syms t1 t2 t3 t4 t5 t6

val = get(hObject,'Value');

q(5) = val;

val = val * 1000;
val = floor(val);
val = val / 1000;
txt = sprintf('%.3f', val);
set(handles.text12,'String',txt);

handles.ARMA.plot(q);
guidata(hObject, handles);

H = eval(T(q(1), q(2), q(3), q(4), q(5), q(6) ));

txt1 = sprintf('%.4f', H(1,4));
txt2 = sprintf('%.4f', H(2,4));
txt3 = sprintf('%.4f', H(3,4));
set(handles.text36,'String',txt1);
set(handles.text37,'String',txt2);
set(handles.text38,'String',txt3);

set(handles.slider7, 'Value', H(1,4));
set(handles.slider8, 'Value', H(2,4));
set(handles.slider9, 'Value', H(3,4));

set(handles.edit1,'String',txt1);
set(handles.edit2,'String',txt2);
set(handles.edit3,'String',txt3);

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider6_Callback(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global q T
syms t1 t2 t3 t4 t5 t6

val = get(hObject,'Value');

q(6) = val;

val = val * 1000;
val = floor(val);
val = val / 1000;
txt = sprintf('%.3f', val);
set(handles.text13,'String',txt);

handles.ARMA.plot(q);
guidata(hObject, handles);

H = eval(T(q(1), q(2), q(3), q(4), q(5), q(6) ));

txt1 = sprintf('%.4f', H(1,4));
txt2 = sprintf('%.4f', H(2,4));
txt3 = sprintf('%.4f', H(3,4));
set(handles.text36,'String',txt1);
set(handles.text37,'String',txt2);
set(handles.text38,'String',txt3);

set(handles.slider7, 'Value', H(1,4));
set(handles.slider8, 'Value', H(2,4));
set(handles.slider9, 'Value', H(3,4));

set(handles.edit1,'String',txt1);
set(handles.edit2,'String',txt2);
set(handles.edit3,'String',txt3);

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function uitable1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on slider movement.
function slider7_Callback(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.text46,'String','Moving');

xi = str2num(get(handles.text36,'String'));
yi = str2num(get(handles.text37,'String'));
zi = str2num(get(handles.text38,'String'));

xf = get(hObject,'Value')

global Ja q T
syms t1 t2 t3 t4 t5 t6

XYZinit=[xi;yi;zi];
XYZfinal=[xf;yi;zi];

xa = xi;

if xf > xa
   while xf > xa
       xabs = abs(xf-xa);

       if xabs > 0.1
           dec = 0.05;
       elseif xabs > 0.01
           dec = 0.007;
       else
           dec = 0.001;
       end
        
       xa = xa + dec;
       H=pinv(eval(Ja(q(1),q(2),q(3),q(4),q(5),q(6))));
       q = transpose( H*(XYZfinal - XYZinit) + transpose(q) );
       Ha = eval(T(q(1), q(2), q(3), q(4), q(5), q(6) ));
       XYZinit(1) = Ha(1,4);
       XYZinit(2) = Ha(2,4);
       XYZinit(3) = Ha(3,4);
       handles.ARMA.plot(q);
       
        txt1 = sprintf('%.4f', Ha(1,4));
        txt2 = sprintf('%.4f', Ha(2,4));
        txt3 = sprintf('%.4f', Ha(3,4));
        set(handles.text36,'String',txt1);
        set(handles.text37,'String',txt2);
        set(handles.text38,'String',txt3);
        
        txt = sprintf('%.3f', q(1));
        set(handles.q1,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(2));
        set(handles.text9,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(3));
        set(handles.text10,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(4));
        set(handles.text11,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(5));
        set(handles.text12,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(6));
        set(handles.text13,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        
        set(handles.edit1,'String',txt1);
        set(handles.edit2,'String',txt2);
        set(handles.edit3,'String',txt3);
              
   end
elseif xf < xa
    while xf < xa
       xabs = abs(xf-xa);

        if xabs > 0.1
            dec = 0.05;
        elseif xabs > 0.01
            dec = 0.007;
        else
            dec = 0.001;
        end
        
       xa = xa - dec ;
       H=pinv(eval(Ja(q(1),q(2),q(3),q(4),q(5),q(6))));
       q = transpose( H*(XYZfinal - XYZinit) + transpose(q) );
       Ha = eval(T(q(1), q(2), q(3), q(4), q(5), q(6) ));
       XYZinit(1) = Ha(1,4);
       XYZinit(2) = Ha(2,4);
       XYZinit(3) = Ha(3,4);
       handles.ARMA.plot(q);
       
        txt1 = sprintf('%.4f', Ha(1,4));
        txt2 = sprintf('%.4f', Ha(2,4));
        txt3 = sprintf('%.4f', Ha(3,4));
        set(handles.text36,'String',txt1);
        set(handles.text37,'String',txt2);
        set(handles.text38,'String',txt3);
        
        txt = sprintf('%.3f', q(1));
        set(handles.q1,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(2));
        set(handles.text9,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(3));
        set(handles.text10,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(4));
        set(handles.text11,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(5));
        set(handles.text12,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(6));
        set(handles.text13,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        
        set(handles.edit1,'String',txt1);
        set(handles.edit2,'String',txt2);
        set(handles.edit3,'String',txt3);
        
    end
end

set(handles.text46,'String','Done !');


% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider8_Callback(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.text46,'String','Moving');

xi = str2num(get(handles.text36,'String'));
yi = str2num(get(handles.text37,'String'));
zi = str2num(get(handles.text38,'String'));

yf = get(hObject,'Value')

global Ja q T
syms t1 t2 t3 t4 t5 t6

XYZinit=[xi;yi;zi];
XYZfinal=[xi;yf;zi];

ya = yi;

if yf > ya
   while yf > ya
       yabs = abs(yf-ya);

       if yabs > 0.1
           dec = 0.05;
       elseif yabs > 0.01
           dec = 0.007;
       else
           dec = 0.001;
       end
        
       ya = ya + dec;
       H=pinv(eval(Ja(q(1),q(2),q(3),q(4),q(5),q(6))));
       q = transpose( H*(XYZfinal - XYZinit) + transpose(q) );
       Ha = eval(T(q(1), q(2), q(3), q(4), q(5), q(6) ));
       XYZinit(1) = Ha(1,4);
       XYZinit(2) = Ha(2,4);
       XYZinit(3) = Ha(3,4);
       handles.ARMA.plot(q);
       
        txt1 = sprintf('%.4f', Ha(1,4));
        txt2 = sprintf('%.4f', Ha(2,4));
        txt3 = sprintf('%.4f', Ha(3,4));
        set(handles.text36,'String',txt1);
        set(handles.text37,'String',txt2);
        set(handles.text38,'String',txt3);
        
        txt = sprintf('%.3f', q(1));
        set(handles.q1,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(2));
        set(handles.text9,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(3));
        set(handles.text10,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(4));
        set(handles.text11,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(5));
        set(handles.text12,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(6));
        set(handles.text13,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        
        set(handles.edit1,'String',txt1);
        set(handles.edit2,'String',txt2);
        set(handles.edit3,'String',txt3);
              
   end
elseif yf < ya
    while yf < ya
       yabs = abs(yf-ya);

        if yabs > 0.1
            dec = 0.05;
        elseif yabs > 0.01
            dec = 0.007;
        else
            dec = 0.001;
        end
        
       ya = ya - dec ;
       H=pinv(eval(Ja(q(1),q(2),q(3),q(4),q(5),q(6))));
       q = transpose( H*(XYZfinal - XYZinit) + transpose(q) );
       Ha = eval(T(q(1), q(2), q(3), q(4), q(5), q(6) ));
       XYZinit(1) = Ha(1,4);
       XYZinit(2) = Ha(2,4);
       XYZinit(3) = Ha(3,4);
       handles.ARMA.plot(q);
       
        txt1 = sprintf('%.4f', Ha(1,4));
        txt2 = sprintf('%.4f', Ha(2,4));
        txt3 = sprintf('%.4f', Ha(3,4));
        set(handles.text36,'String',txt1);
        set(handles.text37,'String',txt2);
        set(handles.text38,'String',txt3);
        
        txt = sprintf('%.3f', q(1));
        set(handles.q1,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(2));
        set(handles.text9,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(3));
        set(handles.text10,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(4));
        set(handles.text11,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(5));
        set(handles.text12,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(6));
        set(handles.text13,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        
        set(handles.edit1,'String',txt1);
        set(handles.edit2,'String',txt2);
        set(handles.edit3,'String',txt3);
        
    end
end

set(handles.text46,'String','Done !');


% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider9_Callback(hObject, eventdata, handles)
% hObject    handle to slider9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.text46,'String','Moving');

xi = str2num(get(handles.text36,'String'));
yi = str2num(get(handles.text37,'String'));
zi = str2num(get(handles.text38,'String'));

zf = get(hObject,'Value')

global Ja q T
syms t1 t2 t3 t4 t5 t6

XYZinit=[xi;yi;zi];
XYZfinal=[xi;yi;zf];

za = zi;

if zf > za
   while zf > za
       zabs = abs(zf-za);

       if zabs > 0.1
           dec = 0.05;
       elseif zabs > 0.01
           dec = 0.007;
       else
           dec = 0.001;
       end
        
       za = za + dec;
       H=pinv(eval(Ja(q(1),q(2),q(3),q(4),q(5),q(6))));
       q = transpose( H*(XYZfinal - XYZinit) + transpose(q) );
       Ha = eval(T(q(1), q(2), q(3), q(4), q(5), q(6) ));
       XYZinit(1) = Ha(1,4);
       XYZinit(2) = Ha(2,4);
       XYZinit(3) = Ha(3,4);
       handles.ARMA.plot(q);
       
        txt1 = sprintf('%.4f', Ha(1,4));
        txt2 = sprintf('%.4f', Ha(2,4));
        txt3 = sprintf('%.4f', Ha(3,4));
        set(handles.text36,'String',txt1);
        set(handles.text37,'String',txt2);
        set(handles.text38,'String',txt3);
        
        txt = sprintf('%.3f', q(1));
        set(handles.q1,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(2));
        set(handles.text9,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(3));
        set(handles.text10,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(4));
        set(handles.text11,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(5));
        set(handles.text12,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(6));
        set(handles.text13,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        
        set(handles.edit1,'String',txt1);
        set(handles.edit2,'String',txt2);
        set(handles.edit3,'String',txt3);
              
   end
elseif zf < za
    while zf < za
       zabs = abs(zf-za);

        if zabs > 0.1
            dec = 0.05;
        elseif zabs > 0.01
            dec = 0.007;
        else
            dec = 0.001;
        end
        
       za = za - dec ;
       H=pinv(eval(Ja(q(1),q(2),q(3),q(4),q(5),q(6))));
       q = transpose( H*(XYZfinal - XYZinit) + transpose(q) );
       Ha = eval(T(q(1), q(2), q(3), q(4), q(5), q(6) ));
       XYZinit(1) = Ha(1,4);
       XYZinit(2) = Ha(2,4);
       XYZinit(3) = Ha(3,4);
       handles.ARMA.plot(q);
       
        txt1 = sprintf('%.4f', Ha(1,4));
        txt2 = sprintf('%.4f', Ha(2,4));
        txt3 = sprintf('%.4f', Ha(3,4));
        set(handles.text36,'String',txt1);
        set(handles.text37,'String',txt2);
        set(handles.text38,'String',txt3);
        
        txt = sprintf('%.3f', q(1));
        set(handles.q1,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(2));
        set(handles.text9,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(3));
        set(handles.text10,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(4));
        set(handles.text11,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(5));
        set(handles.text12,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(6));
        set(handles.text13,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        
        set(handles.edit1,'String',txt1);
        set(handles.edit2,'String',txt2);
        set(handles.edit3,'String',txt3);
        
    end
end

set(handles.text46,'String','Done !');

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.text46,'String','Moving');

xi = str2num(get(handles.text36,'String'));
yi = str2num(get(handles.text37,'String'));
zi = str2num(get(handles.text38,'String'));

xf = str2num(get(handles.edit1,'String'));
yf = str2num(get(handles.edit2,'String'));
zf = str2num(get(handles.edit3,'String'));

global Ja q T
syms t1 t2 t3 t4 t5 t6

XYZinit=[xi;yi;zi];
XYZfinal=[xf;yf;zf];

xabs = abs(xf-xi);
yabs = abs(yf-yi);

xa = xi;
ya = yi;
za = zi;

if zf > za
   while zf > za
       zabs = abs(zf-za);

       if zabs > 0.1
           dec = 0.05;
       elseif zabs > 0.01
           dec = 0.007;
       else
           dec = 0.001;
       end
        
       za = za + dec;
       H=pinv(eval(Ja(q(1),q(2),q(3),q(4),q(5),q(6))));
       q = transpose( H*([xi;yi;za] - XYZinit) + transpose(q) );
       Ha = eval(T(q(1), q(2), q(3), q(4), q(5), q(6) ));
       XYZinit(1) = Ha(1,4);
       XYZinit(2) = Ha(2,4);
       XYZinit(3) = Ha(3,4);
       handles.ARMA.plot(q);
       
        txt1 = sprintf('%.4f', Ha(1,4));
        txt2 = sprintf('%.4f', Ha(2,4));
        txt3 = sprintf('%.4f', Ha(3,4));
        set(handles.text36,'String',txt1);
        set(handles.text37,'String',txt2);
        set(handles.text38,'String',txt3);
        
        txt = sprintf('%.3f', q(1));
        set(handles.q1,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(2));
        set(handles.text9,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(3));
        set(handles.text10,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(4));
        set(handles.text11,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(5));
        set(handles.text12,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(6));
        set(handles.text13,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
              
   end
elseif zf < za
    while zf < za
       zabs = abs(zf-za);

        if zabs > 0.1
            dec = 0.05;
        elseif zabs > 0.01
            dec = 0.007;
        else
            dec = 0.001;
        end
        
       za = za - dec ;
       H=pinv(eval(Ja(q(1),q(2),q(3),q(4),q(5),q(6))));
       q = transpose( H*([xi;yi;za] - XYZinit) + transpose(q) );
       Ha = eval(T(q(1), q(2), q(3), q(4), q(5), q(6) ));
       XYZinit(1) = Ha(1,4);
       XYZinit(2) = Ha(2,4);
       XYZinit(3) = Ha(3,4);
       handles.ARMA.plot(q);
       
        txt1 = sprintf('%.4f', Ha(1,4));
        txt2 = sprintf('%.4f', Ha(2,4));
        txt3 = sprintf('%.4f', Ha(3,4));
        set(handles.text36,'String',txt1);
        set(handles.text37,'String',txt2);
        set(handles.text38,'String',txt3);
        
        txt = sprintf('%.3f', q(1));
        set(handles.q1,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(2));
        set(handles.text9,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(3));
        set(handles.text10,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(4));
        set(handles.text11,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(5));
        set(handles.text12,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(6));
        set(handles.text13,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        
    end
end

if xabs < yabs
    
    if xf > xa
       while xf >= xa
           xabs = abs(xf-xa);

           if xabs > 0.1
               dec = 0.03;
           elseif xabs > 0.01
               dec = 0.004;
           else
               dec = 0.001;
           end

           xa = xa + dec;
           H=pinv(eval(Ja(q(1),q(2),q(3),q(4),q(5),q(6))));
           q = transpose( H*([xa;yi;zf] - XYZinit) + transpose(q) );
           Ha = eval(T(q(1), q(2), q(3), q(4), q(5), q(6) ));
           XYZinit(1) = Ha(1,4);
           XYZinit(2) = Ha(2,4);
           XYZinit(3) = Ha(3,4);
           handles.ARMA.plot(q);
           
            txt1 = sprintf('%.4f', Ha(1,4));
            txt2 = sprintf('%.4f', Ha(2,4));
            txt3 = sprintf('%.4f', Ha(3,4));
            set(handles.text36,'String',txt1);
            set(handles.text37,'String',txt2);
            set(handles.text38,'String',txt3);

            txt = sprintf('%.3f', q(1));
            set(handles.q1,'String',txt);
            val = str2double(txt);
            set(handles.slider1,'Value',val);
            txt = sprintf('%.3f', q(2));
            set(handles.text9,'String',txt);
            val = str2double(txt);
            set(handles.slider1,'Value',val);
            txt = sprintf('%.3f', q(3));
            set(handles.text10,'String',txt);
            val = str2double(txt);
            set(handles.slider1,'Value',val);
            txt = sprintf('%.3f', q(4));
            set(handles.text11,'String',txt);
            val = str2double(txt);
            set(handles.slider1,'Value',val);
            txt = sprintf('%.3f', q(5));
            set(handles.text12,'String',txt);
            val = str2double(txt);
            set(handles.slider1,'Value',val);
            txt = sprintf('%.3f', q(6));
            set(handles.text13,'String',txt);
            val = str2double(txt);
            set(handles.slider1,'Value',val);
        
       end
    elseif xf < xa
        while xf <= xa
           xabs = abs(xf-xa);

            if xabs > 0.1
                dec = 0.03;
            elseif xabs > 0.01
                dec = 0.004;
            else
                dec = 0.001;
            end

           xa = xa - dec ;
           H=pinv(eval(Ja(q(1),q(2),q(3),q(4),q(5),q(6))));
           q = transpose( H*([xa;yi;zf] - XYZinit) + transpose(q) );
           Ha = eval(T(q(1), q(2), q(3), q(4), q(5), q(6) ));
           XYZinit(1) = Ha(1,4);
           XYZinit(2) = Ha(2,4);
           XYZinit(3) = Ha(3,4);
           handles.ARMA.plot(q);
           
            txt1 = sprintf('%.4f', Ha(1,4));
            txt2 = sprintf('%.4f', Ha(2,4));
            txt3 = sprintf('%.4f', Ha(3,4));
            set(handles.text36,'String',txt1);
            set(handles.text37,'String',txt2);
            set(handles.text38,'String',txt3);

            txt = sprintf('%.3f', q(1));
            set(handles.q1,'String',txt);
            val = str2double(txt);
            set(handles.slider1,'Value',val);
            txt = sprintf('%.3f', q(2));
            set(handles.text9,'String',txt);
            val = str2double(txt);
            set(handles.slider1,'Value',val);
            txt = sprintf('%.3f', q(3));
            set(handles.text10,'String',txt);
            val = str2double(txt);
            set(handles.slider1,'Value',val);
            txt = sprintf('%.3f', q(4));
            set(handles.text11,'String',txt);
            val = str2double(txt);
            set(handles.slider1,'Value',val);
            txt = sprintf('%.3f', q(5));
            set(handles.text12,'String',txt);
            val = str2double(txt);
            set(handles.slider1,'Value',val);
            txt = sprintf('%.3f', q(6));
            set(handles.text13,'String',txt);
            val = str2double(txt);
            set(handles.slider1,'Value',val);
           
        end
    end
    
    
    if yf > ya
       while yf >= ya
           yabs = abs(yf-ya);

           if yabs > 0.1
               dec = 0.03;
           elseif yabs > 0.01
               dec = 0.004;
           else
               dec = 0.001;
           end

           ya = ya + dec;
           H=pinv(eval(Ja(q(1),q(2),q(3),q(4),q(5),q(6))));
           q = transpose( H*([xa;ya;zf] - XYZinit) + transpose(q) );
           Ha = eval(T(q(1), q(2), q(3), q(4), q(5), q(6) ));
           XYZinit(1) = Ha(1,4);
           XYZinit(2) = Ha(2,4);
           XYZinit(3) = Ha(3,4);
           handles.ARMA.plot(q);
           
            txt1 = sprintf('%.4f', Ha(1,4));
            txt2 = sprintf('%.4f', Ha(2,4));
            txt3 = sprintf('%.4f', Ha(3,4));
            set(handles.text36,'String',txt1);
            set(handles.text37,'String',txt2);
            set(handles.text38,'String',txt3);

            txt = sprintf('%.3f', q(1));
            set(handles.q1,'String',txt);
            val = str2double(txt);
            set(handles.slider1,'Value',val);
            txt = sprintf('%.3f', q(2));
            set(handles.text9,'String',txt);
            val = str2double(txt);
            set(handles.slider1,'Value',val);
            txt = sprintf('%.3f', q(3));
            set(handles.text10,'String',txt);
            val = str2double(txt);
            set(handles.slider1,'Value',val);
            txt = sprintf('%.3f', q(4));
            set(handles.text11,'String',txt);
            val = str2double(txt);
            set(handles.slider1,'Value',val);
            txt = sprintf('%.3f', q(5));
            set(handles.text12,'String',txt);
            val = str2double(txt);
            set(handles.slider1,'Value',val);
            txt = sprintf('%.3f', q(6));
            set(handles.text13,'String',txt);
            val = str2double(txt);
            set(handles.slider1,'Value',val);
           
       end
    elseif yf < ya
        while yf <= ya
           yabs = abs(yf-ya);

            if yabs > 0.1
                dec = 0.03;
            elseif yabs > 0.01
                dec = 0.004;
            else
                dec = 0.001;
            end

           ya = ya - dec ;
           H=pinv(eval(Ja(q(1),q(2),q(3),q(4),q(5),q(6))));
           q = transpose( H*([xa;ya;zf] - XYZinit) + transpose(q) );
           Ha = eval(T(q(1), q(2), q(3), q(4), q(5), q(6) ));
           XYZinit(1) = Ha(1,4);
           XYZinit(2) = Ha(2,4);
           XYZinit(3) = Ha(3,4);
           handles.ARMA.plot(q);
           
            txt1 = sprintf('%.4f', Ha(1,4));
            txt2 = sprintf('%.4f', Ha(2,4));
            txt3 = sprintf('%.4f', Ha(3,4));
            set(handles.text36,'String',txt1);
            set(handles.text37,'String',txt2);
            set(handles.text38,'String',txt3);

            txt = sprintf('%.3f', q(1));
            set(handles.q1,'String',txt);
            val = str2double(txt);
            set(handles.slider1,'Value',val);
            txt = sprintf('%.3f', q(2));
            set(handles.text9,'String',txt);
            val = str2double(txt);
            set(handles.slider1,'Value',val);
            txt = sprintf('%.3f', q(3));
            set(handles.text10,'String',txt);
            val = str2double(txt);
            set(handles.slider1,'Value',val);
            txt = sprintf('%.3f', q(4));
            set(handles.text11,'String',txt);
            val = str2double(txt);
            set(handles.slider1,'Value',val);
            txt = sprintf('%.3f', q(5));
            set(handles.text12,'String',txt);
            val = str2double(txt);
            set(handles.slider1,'Value',val);
            txt = sprintf('%.3f', q(6));
            set(handles.text13,'String',txt);
            val = str2double(txt);
            set(handles.slider1,'Value',val);
           
        end
    end
    
    
else
    
    if yf > ya
       while yf >= ya
           yabs = abs(yf-ya);

           if yabs > 0.1
               dec = 0.03;
           elseif yabs > 0.01
               dec = 0.004;
           else
               dec = 0.001;
           end

           ya = ya + dec;
           H=pinv(eval(Ja(q(1),q(2),q(3),q(4),q(5),q(6))));
           q = transpose( H*([xi;ya;zf] - XYZinit) + transpose(q) );
           Ha = eval(T(q(1), q(2), q(3), q(4), q(5), q(6) ));
           XYZinit(1) = Ha(1,4);
           XYZinit(2) = Ha(2,4);
           XYZinit(3) = Ha(3,4);
           handles.ARMA.plot(q);
           
        txt1 = sprintf('%.4f', Ha(1,4));
        txt2 = sprintf('%.4f', Ha(2,4));
        txt3 = sprintf('%.4f', Ha(3,4));
        set(handles.text36,'String',txt1);
        set(handles.text37,'String',txt2);
        set(handles.text38,'String',txt3);
        
        txt = sprintf('%.3f', q(1));
        set(handles.q1,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(2));
        set(handles.text9,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(3));
        set(handles.text10,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(4));
        set(handles.text11,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(5));
        set(handles.text12,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(6));
        set(handles.text13,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
           
       end
    elseif yf < ya
        while yf <= ya
           yabs = abs(yf-ya);

            if yabs > 0.1
                dec = 0.03;
            elseif yabs > 0.01
                dec = 0.004;
            else
                dec = 0.001;
            end

           ya = ya - dec ;
           H=pinv(eval(Ja(q(1),q(2),q(3),q(4),q(5),q(6))));
           q = transpose( H*([xi;ya;zf] - XYZinit) + transpose(q) );
           Ha = eval(T(q(1), q(2), q(3), q(4), q(5), q(6) ));
           XYZinit(1) = Ha(1,4);
           XYZinit(2) = Ha(2,4);
           XYZinit(3) = Ha(3,4);
           handles.ARMA.plot(q);
           
        txt1 = sprintf('%.4f', Ha(1,4));
        txt2 = sprintf('%.4f', Ha(2,4));
        txt3 = sprintf('%.4f', Ha(3,4));
        set(handles.text36,'String',txt1);
        set(handles.text37,'String',txt2);
        set(handles.text38,'String',txt3);
        
        txt = sprintf('%.3f', q(1));
        set(handles.q1,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(2));
        set(handles.text9,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(3));
        set(handles.text10,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(4));
        set(handles.text11,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(5));
        set(handles.text12,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(6));
        set(handles.text13,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
           
        end
    end
    
    if xf > xa
       while xf >= xa
           xabs = abs(xf-xa);

           if xabs > 0.1
               dec = 0.03;
           elseif xabs > 0.01
               dec = 0.004;
           else
               dec = 0.001;
           end

           xa = xa + dec;
           H=pinv(eval(Ja(q(1),q(2),q(3),q(4),q(5),q(6))));
           q = transpose( H*([xa;ya;zf] - XYZinit) + transpose(q) );
           Ha = eval(T(q(1), q(2), q(3), q(4), q(5), q(6) ));
           XYZinit(1) = Ha(1,4);
           XYZinit(2) = Ha(2,4);
           XYZinit(3) = Ha(3,4);
           handles.ARMA.plot(q);
           
        txt1 = sprintf('%.4f', Ha(1,4));
        txt2 = sprintf('%.4f', Ha(2,4));
        txt3 = sprintf('%.4f', Ha(3,4));
        set(handles.text36,'String',txt1);
        set(handles.text37,'String',txt2);
        set(handles.text38,'String',txt3);
        
        txt = sprintf('%.3f', q(1));
        set(handles.q1,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(2));
        set(handles.text9,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(3));
        set(handles.text10,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(4));
        set(handles.text11,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(5));
        set(handles.text12,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(6));
        set(handles.text13,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
           
       end
    elseif xf < xa
        while xf <= xa
           xabs = abs(xf-xa);

            if xabs > 0.1
                dec = 0.03;
            elseif xabs > 0.01
                dec = 0.004;
            else
                dec = 0.001;
            end

           xa = xa - dec ;
           H=pinv(eval(Ja(q(1),q(2),q(3),q(4),q(5),q(6))));
           q = transpose( H*([xa;ya;zf] - XYZinit) + transpose(q) );
           Ha = eval(T(q(1), q(2), q(3), q(4), q(5), q(6) ));
           XYZinit(1) = Ha(1,4);
           XYZinit(2) = Ha(2,4);
           XYZinit(3) = Ha(3,4);
           handles.ARMA.plot(q);
           
        txt1 = sprintf('%.4f', Ha(1,4));
        txt2 = sprintf('%.4f', Ha(2,4));
        txt3 = sprintf('%.4f', Ha(3,4));
        set(handles.text36,'String',txt1);
        set(handles.text37,'String',txt2);
        set(handles.text38,'String',txt3);
        
        txt = sprintf('%.3f', q(1));
        set(handles.q1,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(2));
        set(handles.text9,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(3));
        set(handles.text10,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(4));
        set(handles.text11,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(5));
        set(handles.text12,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
        txt = sprintf('%.3f', q(6));
        set(handles.text13,'String',txt);
        val = str2double(txt);
        set(handles.slider1,'Value',val);
           
        end
    end
    
end

xi = str2num(get(handles.text36,'String'));
yi = str2num(get(handles.text37,'String'));
zi = str2num(get(handles.text38,'String'));

xf = str2num(get(handles.edit1,'String'));
yf = str2num(get(handles.edit2,'String'));
zf = str2num(get(handles.edit3,'String'));

XYZinit=[xi;yi;zi];
XYZfinal=[xf;yf;zf];

H=pinv(eval(Ja(q(1),q(2),q(3),q(4),q(5),q(6))));
q = transpose( H*(XYZfinal - XYZinit) + transpose(q) );
Ha = eval(T(q(1), q(2), q(3), q(4), q(5), q(6) ));

txt1 = sprintf('%.4f', Ha(1,4));
txt2 = sprintf('%.4f', Ha(2,4));
txt3 = sprintf('%.4f', Ha(3,4));
set(handles.text36,'String',txt1);
set(handles.text37,'String',txt2);
set(handles.text38,'String',txt3);

txt = sprintf('%.3f', q(1));
set(handles.q1,'String',txt);
val = str2double(txt);
set(handles.slider1,'Value',val);
txt = sprintf('%.3f', q(2));
set(handles.text9,'String',txt);
val = str2double(txt);
set(handles.slider1,'Value',val);
txt = sprintf('%.3f', q(3));
set(handles.text10,'String',txt);
val = str2double(txt);
set(handles.slider1,'Value',val);
txt = sprintf('%.3f', q(4));
set(handles.text11,'String',txt);
val = str2double(txt);
set(handles.slider1,'Value',val);
txt = sprintf('%.3f', q(5));
set(handles.text12,'String',txt);
val = str2double(txt);
set(handles.slider1,'Value',val);
txt = sprintf('%.3f', q(6));
set(handles.text13,'String',txt);
val = str2double(txt);
set(handles.slider1,'Value',val);

set(handles.text46,'String','Done !');

handles.ARMA.plot(q);


% H=pinv(eval(Ja(q(1),q(2),q(3),q(4),q(5),q(6))));
% q = transpose( H*(XYZfinal - XYZinit) + transpose(q) )

    




% pos = transl(xf*1000,-yf*1000,zf*1000) %- transl(xi*1000,-yi*1000,zi*1000);
% 
% P04K0 = [pos(1,4) - 200*pos(1,3) ; pos(2,4) - 200*pos(2,3) ; pos(3,4) - 200*pos(3,3)];
% 
% teta1_1 = (atan2(pos(2,4)- (0.2*pos(2,3)),pos(1,4)-(0.2*pos(1,3))))
% teta1_2 = (atan2(pos(2,4)- (0.2*pos(2,3)),pos(1,4)-(0.2*pos(1,3)))) + pi;
% 
% l12 = 219^2 + 21^2;
% 
% tab1 = T01
% tab2 = T12
% 
% tab1(t1) = T01
% tab2(t2) = T12
% T02 = tab1(teta1_2) * tab2(0)
% 
% P02K0 = [T02(1,4) ; T02(2,4) ; T02(3,4)];
% 
% 
% P24K0 = P04K0 - P02K0 ;
% 
% nom1 = l12 - 190^2 + (( norm(P24K0) )^2);
% den1 = 2 * norm(P24K0) * sqrt(l12);
% nom1/den1
% 
% nom2 = norm(P24K0) - ( (l12 - 190^2 + (( norm(P24K0) )^2) ) / (2 * norm(P24K0)) );
% 
% Mu = asin( nom1 / den1 ) + asin( nom2 / 190 )
% eval(Mu)
% alp = atan2(219,21);
% 
% teta3_1 = pi - Mu - alp
% teta3_2 = pi + Mu - alp
% 
% eval(teta3_1)
% eval(teta3_2)


% taylor(sin(t6))*(taylor(cos(t4))*taylor(sin(t1)) + taylor(sin(t4))*(taylor(cos(t1))*taylor(sin(t2))*taylor(sin(t3)) - taylor(cos(t1))*taylor(cos(t2))*taylor(cos(t3)))) + taylor(cos(t6))*(taylor(cos(t5))*(taylor(sin(t1))*taylor(sin(t4)) - taylor(cos(t4))*(taylor(cos(t1))*taylor(sin(t2))*taylor(sin(t3)) - taylor(cos(t1))*taylor(cos(t2))*taylor(cos(t3)))) - taylor(sin(t5))*(taylor(cos(t1))*taylor(cos(t2))*taylor(sin(t3)) + taylor(cos(t1))*taylor(cos(t3))*taylor(sin(t2)))) == 1

function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.text46,'String','Resetting ..');

global q T
syms t1 t2 t3 t4 t5 t6

q = [0 -1.57 0 0 0 0];
handles.ARMA.plot(q);

H = eval(T(q(1), q(2), q(3), q(4), q(5), q(6) ));

txt1 = sprintf('%.4f', H(1,4));
txt2 = sprintf('%.4f', H(2,4));
txt3 = sprintf('%.4f', H(3,4));
set(handles.text36,'String',txt1);
set(handles.text37,'String',txt2);
set(handles.text38,'String',txt3);

set(handles.slider7, 'Value', H(1,4));
set(handles.slider8, 'Value', H(2,4));
set(handles.slider9, 'Value', H(3,4));

set(handles.edit1,'String',txt1);
set(handles.edit2,'String',txt2);
set(handles.edit3,'String',txt3);

set(handles.slider1, 'Value', q(1));
set(handles.slider2, 'Value', q(2));
set(handles.slider3, 'Value', q(3));
set(handles.slider4, 'Value', q(4));
set(handles.slider5, 'Value', q(5));
set(handles.slider6, 'Value', q(6));

set(handles.text46,'String','Done !');
