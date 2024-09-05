% Load data
% clear all;
% clc;
% close
load('MI subj.mat');
%%%%%%%% RIGHT HAND OF SUBJECT aa %%%%%%%%

sign1=[RH_aa];
 %% test signal

x=sign1;
% x = x';                         % Ensure x is a row vector
fs = 1;                       % Sampling rate (samples/second)
[n m]=size(x);
x1 = x/max(abs(x));              % Normalize so maximum value is 1

% figure(1)
% clf
% subplot(2,1,1)
% plot((1:n)/fs,x1(:,1))
% box off
% xlim([0 n/fs])
% ylim([-1.5 1.5])
% title('TEST RIGHT HAND EEG SIGNAL OF SUBJECT aa')
% xlabel('TIME (SECONDS)')
%%
%%%%segmentation
sign11=sign1(1:20347,:);
sign12=sign1(20348:40694,:);
sign13=sign1(40695:61041,:);
sign14=sign1(61042:81389,:);


%% Set parameters

p = 5; q = 6; r = 2;                % Parameters
J = 19;                              % Number of levels
h=6;                                % p of the Autoregrssive technique
red = 1/r * 1/(1-p/q);              % Redundancy

% Add path to functions
addpath radwt_functions_v2          % use _v1 or _v2
addpath extra_functions

%% Plot frequency responses of RADWT

% subplot(2,1,2)
% Plot_FreqResps(p,q,s,J,fs);
% orient landscape
% print -dpdf figures/demo3_fig1

%% Compute wavelet representation of the first segment of signl1
%%%%%%%%%%%% 1st Segment of MI %%%%%%%%
x=sign11;
N=length(x);
x=x';
[n m]=size(x);
for i=1:n             
  w = radwt(x(i,:), J, p, q, r);
  y=iradwt(w,p,q,r);
  y=real(y(1:N));
  Y(i,:)=y;
end
Y=Y';
%%
%%%%%%%%%%%% 2nd Segment of MI %%%%%%%
x2=sign12;
N2=length(x2);
x2=x2';
[n2 m2]=size(x2);
for i=1:n2             
  w2 = radwt(x2(i,:), J, p, q, r);
  y2=iradwt(w2,p,q,r);
  y2=real(y2(1:N2));
  Y2(i,:)=y2;
end
Y2=Y2';
%%
%%%%%%%%%%%% 3rd Segment of MI %%%%%%
x3=sign13;
N3=length(x3);
x3=x3';
[n3 m3]=size(x3);
for i=1:n3             
  w3 = radwt(x3(i,:), J, p, q, r);
  y3=iradwt(w3,p,q,r);
  y3=real(y3(1:N2));
  Y3(i,:)=y3;
end
Y3=Y3';
%%
%%%%%%%%% 4th Segment of MI %%%%%%%%
x4=sign14;
N4=length(x4);
x4=x4';
[n4 m4]=size(x4);
for i=1:n4             
  w4 = radwt(x4(i,:), J, p, q, r);
  y4=iradwt(w4,p,q,r);
  y4=real(y4(1:N4));
  Y4(i,:)=y4;
end
Y4=Y4';
%% using Autoregrasive
AR1=pcov(Y,h);
AR2=pcov(Y2,h);
AR3=pcov(Y3,h);
AR4=pcov(Y4,h);
AR_RHaa=[AR1;AR2;AR3;AR4];
%%
clear all;
clc;
load('MI subj.mat');

%%%%%%%%%% RIGHT FOOT OF SUBJECT aa %%%%%%%%%

sign1=[RF_aa];
%% test signal

x=sign1;
% x = x';                         % Ensure x is a row vector
fs = 0.01;                       % Sampling rate (samples/second)
[n m]=size(x);
x1 = x/max(abs(x));              % Normalize so maximum value is 1

% figure(2)
% clf
% subplot(2,1,1)
% plot((1:n)/fs,x1(:,1))
% box off
% xlim([0 n/fs])
% ylim([-1.5 1.5])
% title('TEST RIGHT FOOT EEG SIGNAL OF SUBJECT aa')
% xlabel('TIME (SECONDS)')
%%
sign11=sign1(1:27455,:);
sign12=sign1(27456:54910,:);
sign13=sign1(54911:82365,:);
sign14=sign1(82366:109823,:);

%% Set parameters

p = 5; q = 6; r = 2;                % Parameters
J = 19;                             % Number of levels
h=6;                                % p of the Autoregrssive technique
red = 1/r * 1/(1-p/q);              % Redundancy
% Add path to functions
addpath radwt_functions_v2          % use _v1 or _v2
addpath extra_functions

%% Plot frequency responses of RADWT
% subplot(2,1,2)
% Plot_FreqResps(p,q,s,J,fs);
% orient landscape
% print -dpdf figures/demo3_fig1

%% Compute wavelet representation of the first segment of signl1

%%%%%%%%%%%% 1st Segment of MI %%%%%%%%
x=sign11;
N=length(x);
x=x';
[n m]=size(x);
for i=1:n             
  w = radwt(x(i,:), J, p, q, r);
  y=iradwt(w,p,q,r);
  y=real(y(1:N));
  Y(i,:)=y;
end
Y=Y';
%%
%%%%%%%%%%%% 2nd Segment of MI %%%%%%%
x2=sign12;
N2=length(x2);
x2=x2';
[n2 m2]=size(x2);
for i=1:n2             
  w2 = radwt(x2(i,:), J, p, q, r);
  y2=iradwt(w2,p,q,r);
  y2=real(y2(1:N2));
  Y2(i,:)=y2;
end
Y2=Y2';
%%
%%%%%%%%%%%% 3rd Segment of MI %%%%%%
x3=sign13;
N3=length(x3);
x3=x3';
[n3 m3]=size(x3);
for i=1:n3             
  w3 = radwt(x3(i,:), J, p, q, r);
  y3=iradwt(w3,p,q,r);
  y3=real(y3(1:N2));
  Y3(i,:)=y3;
end
Y3=Y3';
%%
%%%%%%%%% 4th Segment of MI %%%%%%%%
x4=sign14;
N4=length(x4);
x4=x4';
[n4 m4]=size(x4);
for i=1:n4             
  w4 = radwt(x4(i,:), J, p, q, r);
  y4=iradwt(w4,p,q,r);
  y4=real(y4(1:N4));
  Y4(i,:)=y4;
end
Y4=Y4';
%% using Autoregrasive
AR1=pcov(Y,h);
AR2=pcov(Y2,h);
AR3=pcov(Y3,h);
AR4=pcov(Y4,h);
AR_RFaa=[AR1;AR2;AR3;AR4];

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RIGHT HAND OF SUBJECT al %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;
load('MI subj.mat');

%%%%%%%% RIGHT HAND OF SUBJECT al %%%%%%%%

sign1=[RH_al];
 %% test signal

x=sign1;
% x = x';                         % Ensure x is a row vector
fs = 1;                       % Sampling rate (samples/second)
[n m]=size(x);
x1 = x/max(abs(x));              % Normalize so maximum value is 1

% figure(3)
% clf
% subplot(2,1,1)
% plot((1:n)/fs,x1(:,1))
% box off
% xlim([0 n/fs])
% ylim([-1.5 1.5])
% title('TEST RIGHT HAND EEG SIGNAL OF SUBJECT al')
% xlabel('TIME (SECONDS)')
%%
%%%% Segmentation
sign11=sign1(1:20347,:);
sign12=sign1(20348:40694,:);
sign13=sign1(40695:61041,:);
sign14=sign1(61042:81389,:);


%% Set parameters

p = 5; q = 6; r = 2;                % Parameters
J = 19;                              % Number of levels
h=6;                                % p of the Autoregrssive technique
red = 1/r * 1/(1-p/q);              % Redundancy

% Add path to functions
addpath radwt_functions_v2          % use _v1 or _v2
addpath extra_functions

%% Plot frequency responses of RADWT

% subplot(2,1,2)
% Plot_FreqResps(p,q,s,J,fs);
% orient landscape
% print -dpdf figures/demo3_fig1

%% Compute wavelet representation of the first segment of signl1
%%%%%%%%%%%% 1st Segment of MI %%%%%%%%
x=sign11;
N=length(x);
x=x';
[n m]=size(x);
for i=1:n             
  w = radwt(x(i,:), J, p, q, r);
  y=iradwt(w,p,q,r);
  y=real(y(1:N));
  Y(i,:)=y;
end
Y=Y';
%%
%%%%%%%%%%%% 2nd Segment of MI %%%%%%%
x2=sign12;
N2=length(x2);
x2=x2';
[n2 m2]=size(x2);
for i=1:n2             
  w2 = radwt(x2(i,:), J, p, q, r);
  y2=iradwt(w2,p,q,r);
  y2=real(y2(1:N2));
  Y2(i,:)=y2;
end
Y2=Y2';
%%
%%%%%%%%%%%% 3rd Segment of MI %%%%%%
x3=sign13;
N3=length(x3);
x3=x3';
[n3 m3]=size(x3);
for i=1:n3             
  w3 = radwt(x3(i,:), J, p, q, r);
  y3=iradwt(w3,p,q,r);
  y3=real(y3(1:N2));
  Y3(i,:)=y3;
end
Y3=Y3';
%%
%%%%%%%%% 4th Segment of MI %%%%%%%%
x4=sign14;
N4=length(x4);
x4=x4';
[n4 m4]=size(x4);
for i=1:n4             
  w4 = radwt(x4(i,:), J, p, q, r);
  y4=iradwt(w4,p,q,r);
  y4=real(y4(1:N4));
  Y4(i,:)=y4;
end
Y4=Y4';
%% using Autoregrasive
AR1=pcov(Y,h);
AR2=pcov(Y2,h);
AR3=pcov(Y3,h);
AR4=pcov(Y4,h);
AR_RHal=[AR1;AR2;AR3;AR4];
%%
clear all;
clc;
load('MI subj.mat');

%%%%%%%%%% RIGHT FOOT OF SUBJECT al %%%%%%%%%
sign1=[RF_al];
%% test signal

x=sign1;
% x = x';                         % Ensure x is a row vector
fs = 0.01;                       % Sampling rate (samples/second)
[n m]=size(x);
x1 = x/max(abs(x));              % Normalize so maximum value is 1

% figure(4)
% clf
% subplot(2,1,1)
% plot((1:n)/fs,x1(:,1))
% box off
% xlim([0 n/fs])
% ylim([-1.5 1.5])
% title('TEST RIGHT FOOT EEG SIGNAL OF SUBJECT al')
% xlabel('TIME (SECONDS)')
%%
sign11=sign1(1:27455,:);
sign12=sign1(27456:54910,:);
sign13=sign1(54911:82365,:);
sign14=sign1(82366:109823,:);

%% Set parameters

p = 5; q = 6; r = 2;                % Parameters
J = 19;                             % Number of levels
h=6;                                % p of the Autoregrssive technique
red = 1/r * 1/(1-p/q);              % Redundancy
% Add path to functions
addpath radwt_functions_v2          % use _v1 or _v2
addpath extra_functions

%% Plot frequency responses of RADWT
subplot(2,1,2)
% Plot_FreqResps(p,q,s,J,fs);
% orient landscape
% print -dpdf figures/demo3_fig1

%% Compute wavelet representation of the first segment of signl1

%%%%%%%%%%%% 1st Segment of MI %%%%%%%%
x=sign11;
N=length(x);
x=x';
[n m]=size(x);
for i=1:n             
  w = radwt(x(i,:), J, p, q, r);
  y=iradwt(w,p,q,r);
  y=real(y(1:N));
  Y(i,:)=y;
end
Y=Y';
%%
%%%%%%%%%%%% 2nd Segment of MI %%%%%%%
x2=sign12;
N2=length(x2);
x2=x2';
[n2 m2]=size(x2);
for i=1:n2             
  w2 = radwt(x2(i,:), J, p, q, r);
  y2=iradwt(w2,p,q,r);
  y2=real(y2(1:N2));
  Y2(i,:)=y2;
end
Y2=Y2';
%%
%%%%%%%%%%%% 3rd Segment of MI %%%%%%
x3=sign13;
N3=length(x3);
x3=x3';
[n3 m3]=size(x3);
for i=1:n3             
  w3 = radwt(x3(i,:), J, p, q, r);
  y3=iradwt(w3,p,q,r);
  y3=real(y3(1:N2));
  Y3(i,:)=y3;
end
Y3=Y3';
%%
%%%%%%%%% 4th Segment of MI %%%%%%%%
x4=sign14;
N4=length(x4);
x4=x4';
[n4 m4]=size(x4);
for i=1:n4             
  w4 = radwt(x4(i,:), J, p, q, r);
  y4=iradwt(w4,p,q,r);
  y4=real(y4(1:N4));
  Y4(i,:)=y4;
end
Y4=Y4';
%% using Autoregrasive
AR1=pcov(Y,h);
AR2=pcov(Y2,h);
AR3=pcov(Y3,h);
AR4=pcov(Y4,h);
AR_RFal=[AR1;AR2;AR3;AR4];

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% RIGHT HAND OF SUBJECT av %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
load('MI subj.mat');
%%%%%%%% RIGHT HAND OF SUBJECT av %%%%%%%%

sign1=[RH_av];
 %% test signal

x=sign1;
% x = x';                         % Ensure x is a row vector
fs = 1;                       % Sampling rate (samples/second)
[n m]=size(x);
x1 = x/max(abs(x));              % Normalize so maximum value is 1

% figure(5)
% clf
% subplot(2,1,1)
% plot((1:n)/fs,x1(:,1))
% box off
% xlim([0 n/fs])
% ylim([-1.5 1.5])
% title('TEST RIGHT HAND EEG SIGNAL OF SUBJECT av')
% xlabel('TIME (SECONDS)')
%%
%%%%segmentation
sign11=sign1(1:20347,:);
sign12=sign1(20348:40694,:);
sign13=sign1(40695:61041,:);
sign14=sign1(61042:81389,:);


%% Set parameters

p = 5; q = 6; r = 2;                % Parameters
J = 19;                              % Number of levels
h=6;                                % p of the Autoregrssive technique
red = 1/r * 1/(1-p/q);              % Redundancy

% Add path to functions
addpath radwt_functions_v2          % use _v1 or _v2
addpath extra_functions

%% Plot frequency responses of RADWT

% subplot(2,1,2)
% Plot_FreqResps(p,q,s,J,fs);
% orient landscape
% print -dpdf figures/demo3_fig1

%% Compute wavelet representation of the first segment of signl1
%%%%%%%%%%%% 1st Segment of MI %%%%%%%%
x=sign11;
N=length(x);
x=x';
[n m]=size(x);
for i=1:n             
  w = radwt(x(i,:), J, p, q, r);
  y=iradwt(w,p,q,r);
  y=real(y(1:N));
  Y(i,:)=y;
end
Y=Y';
%%
%%%%%%%%%%%% 2nd Segment of MI %%%%%%%
x2=sign12;
N2=length(x2);
x2=x2';
[n2 m2]=size(x2);
for i=1:n2             
  w2 = radwt(x2(i,:), J, p, q, r);
  y2=iradwt(w2,p,q,r);
  y2=real(y2(1:N2));
  Y2(i,:)=y2;
end
Y2=Y2';
%%
%%%%%%%%%%%% 3rd Segment of MI %%%%%%
x3=sign13;
N3=length(x3);
x3=x3';
[n3 m3]=size(x3);
for i=1:n3             
  w3 = radwt(x3(i,:), J, p, q, r);
  y3=iradwt(w3,p,q,r);
  y3=real(y3(1:N2));
  Y3(i,:)=y3;
end
Y3=Y3';
%%
%%%%%%%%% 4th Segment of MI %%%%%%%%
x4=sign14;
N4=length(x4);
x4=x4';
[n4 m4]=size(x4);
for i=1:n4             
  w4 = radwt(x4(i,:), J, p, q, r);
  y4=iradwt(w4,p,q,r);
  y4=real(y4(1:N4));
  Y4(i,:)=y4;
end
Y4=Y4';
%% using Autoregrasive
AR1=pcov(Y,h);
AR2=pcov(Y2,h);
AR3=pcov(Y3,h);
AR4=pcov(Y4,h);
AR_RHav=[AR1;AR2;AR3;AR4];
%%
clear all;
clc;
load('MI subj.mat');

%%%%%%%%%% RIGHT FOOT OF SUBJECT al %%%%%%%%%
sign1=[RF_av];
%% test signal

x=sign1;
% x = x';                         % Ensure x is a row vector
fs = 0.01;                       % Sampling rate (samples/second)
[n m]=size(x);
x1 = x/max(abs(x));              % Normalize so maximum value is 1

% figure(6)
% clf
% subplot(2,1,1)
% plot((1:n)/fs,x1(:,1))
% box off
% xlim([0 n/fs])
% ylim([-1.5 1.5])
% title('TEST RIGHT FOOT EEG SIGNAL OF SUBJECT av')
% xlabel('TIME (SECONDS)')
%%
sign11=sign1(1:27455,:);
sign12=sign1(27456:54910,:);
sign13=sign1(54911:82365,:);
sign14=sign1(82366:109823,:);

%% Set parameters

p = 5; q = 6; r = 2;                % Parameters
J = 19;                             % Number of levels
h=6;                                % p of the Autoregrssive technique
red = 1/r * 1/(1-p/q);              % Redundancy
% Add path to functions
addpath radwt_functions_v2          % use _v1 or _v2
addpath extra_functions

%% Plot frequency responses of RADWT
% subplot(2,1,2)
% Plot_FreqResps(p,q,s,J,fs);
% orient landscape
% print -dpdf figures/demo3_fig1

%% Compute wavelet representation of the first segment of signl1

%%%%%%%%%%%% 1st Segment of MI %%%%%%%%
x=sign11;
N=length(x);
x=x';
[n m]=size(x);
for i=1:n             
  w = radwt(x(i,:), J, p, q, r);
  y=iradwt(w,p,q,r);
  y=real(y(1:N));
  Y(i,:)=y;
end
Y=Y';
%%
%%%%%%%%%%%% 2nd Segment of MI %%%%%%%
x2=sign12;
N2=length(x2);
x2=x2';
[n2 m2]=size(x2);
for i=1:n2             
  w2 = radwt(x2(i,:), J, p, q, r);
  y2=iradwt(w2,p,q,r);
  y2=real(y2(1:N2));
  Y2(i,:)=y2;
end
Y2=Y2';
%%
%%%%%%%%%%%% 3rd Segment of MI %%%%%%
x3=sign13;
N3=length(x3);
x3=x3';
[n3 m3]=size(x3);
for i=1:n3             
  w3 = radwt(x3(i,:), J, p, q, r);
  y3=iradwt(w3,p,q,r);
  y3=real(y3(1:N2));
  Y3(i,:)=y3;
end
Y3=Y3';
%%
%%%%%%%%% 4th Segment of MI %%%%%%%%
x4=sign14;
N4=length(x4);
x4=x4';
[n4 m4]=size(x4);
for i=1:n4             
  w4 = radwt(x4(i,:), J, p, q, r);
  y4=iradwt(w4,p,q,r);
  y4=real(y4(1:N4));
  Y4(i,:)=y4;
end
Y4=Y4';
%% using Autoregrasive
AR1=pcov(Y,h);
AR2=pcov(Y2,h);
AR3=pcov(Y3,h);
AR4=pcov(Y4,h);
AR_RFav=[AR1;AR2;AR3;AR4];
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% RIGHT HAND OF SUBJECT aw %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
load('MI subj.mat');
%%%%%%%% RIGHT HAND OF SUBJECT aw %%%%%%%%

sign1=[RH_aw];
 %% test signal

x=sign1;
% x = x';                         % Ensure x is a row vector
fs = 1;                       % Sampling rate (samples/second)
[n m]=size(x);
x1 = x/max(abs(x));              % Normalize so maximum value is 1
% 
% figure(7)
% clf
% subplot(2,1,1)
% plot((1:n)/fs,x1(:,1))
% box off
% xlim([0 n/fs])
% ylim([-1.5 1.5])
% title('TEST RIGHT HAND EEG SIGNAL OF SUBJECT aw')
% xlabel('TIME (SECONDS)')
%%
%%%%segmentation
sign11=sign1(1:20347,:);
sign12=sign1(20348:40694,:);
sign13=sign1(40695:61041,:);
sign14=sign1(61042:81389,:);


%% Set parameters

p = 5; q = 6; r = 2;                % Parameters
J = 19;                              % Number of levels
h=6;                                % p of the Autoregrssive technique
red = 1/r * 1/(1-p/q);              % Redundancy

% Add path to functions
addpath radwt_functions_v2          % use _v1 or _v2
addpath extra_functions

%% Plot frequency responses of RADWT

% subplot(2,1,2)
% Plot_FreqResps(p,q,s,J,fs);
% orient landscape
% print -dpdf figures/demo3_fig1

%% Compute wavelet representation of the first segment of signl1
%%%%%%%%%%%% 1st Segment of MI %%%%%%%%
x=sign11;
N=length(x);
x=x';
[n m]=size(x);
for i=1:n             
  w = radwt(x(i,:), J, p, q, r);
  y=iradwt(w,p,q,r);
  y=real(y(1:N));
  Y(i,:)=y;
end
Y=Y';
%%
%%%%%%%%%%%% 2nd Segment of MI %%%%%%%
x2=sign12;
N2=length(x2);
x2=x2';
[n2 m2]=size(x2);
for i=1:n2             
  w2 = radwt(x2(i,:), J, p, q, r);
  y2=iradwt(w2,p,q,r);
  y2=real(y2(1:N2));
  Y2(i,:)=y2;
end
Y2=Y2';
%%
%%%%%%%%%%%% 3rd Segment of MI %%%%%%
x3=sign13;
N3=length(x3);
x3=x3';
[n3 m3]=size(x3);
for i=1:n3             
  w3 = radwt(x3(i,:), J, p, q, r);
  y3=iradwt(w3,p,q,r);
  y3=real(y3(1:N2));
  Y3(i,:)=y3;
end
Y3=Y3';
%%
%%%%%%%%% 4th Segment of MI %%%%%%%%
x4=sign14;
N4=length(x4);
x4=x4';
[n4 m4]=size(x4);
for i=1:n4             
  w4 = radwt(x4(i,:), J, p, q, r);
  y4=iradwt(w4,p,q,r);
  y4=real(y4(1:N4));
  Y4(i,:)=y4;
end
Y4=Y4';
%% using Autoregrasive
AR1=pcov(Y,h);
AR2=pcov(Y2,h);
AR3=pcov(Y3,h);
AR4=pcov(Y4,h);
AR_RHaw=[AR1;AR2;AR3;AR4];
%%
clear all;
clc;
load('MI subj.mat');

%%%%%%%%%% RIGHT FOOT OF SUBJECT aw %%%%%%%%%
sign1=[RF_aw];
%% test signal

x=sign1;
% x = x';                         % Ensure x is a row vector
fs = 0.01;                       % Sampling rate (samples/second)
[n m]=size(x);
x1 = x/max(abs(x));              % Normalize so maximum value is 1
% 
% figure(8)
% clf
% subplot(2,1,1)
% plot((1:n)/fs,x1(:,1))
% box off
% xlim([0 n/fs])
% ylim([-1.5 1.5])
% title('TEST RIGHT FOOT EEG SIGNAL OF SUBJECT aw')
% xlabel('TIME (SECONDS)')
%%
sign11=sign1(1:27455,:);
sign12=sign1(27456:54910,:);
sign13=sign1(54911:82365,:);
sign14=sign1(82366:109823,:);

%% Set parameters

p = 5; q = 6; r = 2;                % Parameters
J = 19;                             % Number of levels
h=6;                                % p of the Autoregrssive technique
red = 1/r * 1/(1-p/q);              % Redundancy
% Add path to functions
addpath radwt_functions_v2          % use _v1 or _v2
addpath extra_functions

%% Plot frequency responses of RADWT
% subplot(2,1,2)
% Plot_FreqResps(p,q,s,J,fs);
% orient landscape
% print -dpdf figures/demo3_fig1

%% Compute wavelet representation of the first segment of signl1

%%%%%%%%%%%% 1st Segment of MI %%%%%%%%
x=sign11;
N=length(x);
x=x';
[n m]=size(x);
for i=1:n             
  w = radwt(x(i,:), J, p, q, r);
  y=iradwt(w,p,q,r);
  y=real(y(1:N));
  Y(i,:)=y;
end
Y=Y';
%%
%%%%%%%%%%%% 2nd Segment of MI %%%%%%%
x2=sign12;
N2=length(x2);
x2=x2';
[n2 m2]=size(x2);
for i=1:n2             
  w2 = radwt(x2(i,:), J, p, q, r);
  y2=iradwt(w2,p,q,r);
  y2=real(y2(1:N2));
  Y2(i,:)=y2;
end
Y2=Y2';
%%
%%%%%%%%%%%% 3rd Segment of MI %%%%%%
x3=sign13;
N3=length(x3);
x3=x3';
[n3 m3]=size(x3);
for i=1:n3             
  w3 = radwt(x3(i,:), J, p, q, r);
  y3=iradwt(w3,p,q,r);
  y3=real(y3(1:N2));
  Y3(i,:)=y3;
end
Y3=Y3';
%%
%%%%%%%%% 4th Segment of MI %%%%%%%%
x4=sign14;
N4=length(x4);
x4=x4';
[n4 m4]=size(x4);
for i=1:n4             
  w4 = radwt(x4(i,:), J, p, q, r);
  y4=iradwt(w4,p,q,r);
  y4=real(y4(1:N4));
  Y4(i,:)=y4;
end
Y4=Y4';
%% using Autoregrasive
AR1=pcov(Y,h);
AR2=pcov(Y2,h);
AR3=pcov(Y3,h);
AR4=pcov(Y4,h);
AR_RFaw=[AR1;AR2;AR3;AR4];
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% RIGHT HAND OF SUBJECT ay %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
load('MI subj.mat');
%%%%%%%% RIGHT HAND OF SUBJECT ay %%%%%%%%

sign1=[RH_ay];
 %% test signal

x=sign1;
% x = x';                         % Ensure x is a row vector
fs = 1;                       % Sampling rate (samples/second)
[n m]=size(x);
x1 = x/max(abs(x));              % Normalize so maximum value is 1

% figure(9)
% clf
% subplot(2,1,1)
% plot((1:n)/fs,x1(:,1))
% box off
% xlim([0 n/fs])
% ylim([-1.5 1.5])
% title('TEST RIGHT HAND EEG SIGNAL OF SUBJECT ay')
% xlabel('TIME (SECONDS)')
%%
%%%%segmentation
sign11=sign1(1:20347,:);
sign12=sign1(20348:40694,:);
sign13=sign1(40695:61041,:);
sign14=sign1(61042:81389,:);


%% Set parameters

p = 5; q = 6; r = 2;                % Parameters
J = 19;                              % Number of levels
h=6;                                % p of the Autoregrssive technique
red = 1/r * 1/(1-p/q);              % Redundancy

% Add path to functions
addpath radwt_functions_v2          % use _v1 or _v2
addpath extra_functions

%% Plot frequency responses of RADWT

% subplot(2,1,2)
% Plot_FreqResps(p,q,s,J,fs);
% orient landscape
% print -dpdf figures/demo3_fig1

%% Compute wavelet representation of the first segment of signl1
%%%%%%%%%%%% 1st Segment of MI %%%%%%%%
x=sign11;
N=length(x);
x=x';
[n m]=size(x);
for i=1:n             
  w = radwt(x(i,:), J, p, q, r);
  y=iradwt(w,p,q,r);
  y=real(y(1:N));
  Y(i,:)=y;
end
Y=Y';
%%
%%%%%%%%%%%% 2nd Segment of MI %%%%%%%
x2=sign12;
N2=length(x2);
x2=x2';
[n2 m2]=size(x2);
for i=1:n2             
  w2 = radwt(x2(i,:), J, p, q, r);
  y2=iradwt(w2,p,q,r);
  y2=real(y2(1:N2));
  Y2(i,:)=y2;
end
Y2=Y2';
%%
%%%%%%%%%%%% 3rd Segment of MI %%%%%%
x3=sign13;
N3=length(x3);
x3=x3';
[n3 m3]=size(x3);
for i=1:n3             
  w3 = radwt(x3(i,:), J, p, q, r);
  y3=iradwt(w3,p,q,r);
  y3=real(y3(1:N2));
  Y3(i,:)=y3;
end
Y3=Y3';
%%
%%%%%%%%% 4th Segment of MI %%%%%%%%
x4=sign14;
N4=length(x4);
x4=x4';
[n4 m4]=size(x4);
for i=1:n4             
  w4 = radwt(x4(i,:), J, p, q, r);
  y4=iradwt(w4,p,q,r);
  y4=real(y4(1:N4));
  Y4(i,:)=y4;
end
Y4=Y4';
%% using Autoregrasive
AR1=pcov(Y,h);
AR2=pcov(Y2,h);
AR3=pcov(Y3,h);
AR4=pcov(Y4,h);
AR_RHay=[AR1;AR2;AR3;AR4];
%%
clear all;
clc;
load('MI subj.mat');

%%%%%%%%%% RIGHT FOOT OF SUBJECT ay %%%%%%%%%

sign1=[RF_ay];
%% test signal

x=sign1;
% x = x';                         % Ensure x is a row vector
fs = 0.01;                       % Sampling rate (samples/second)
[n m]=size(x);
x1 = x/max(abs(x));              % Normalize so maximum value is 1

% figure(10)
% clf
% subplot(2,1,1)
% plot((1:n)/fs,x1(:,1))
% box off
% xlim([0 n/fs])
% ylim([-1.5 1.5])
% title('TEST RIGHT FOOT EEG SIGNAL OF SUBJECT ay')
% xlabel('TIME (SECONDS)')
%%
sign11=sign1(1:27455,:);
sign12=sign1(27456:54910,:);
sign13=sign1(54911:82365,:);
sign14=sign1(82366:109823,:);

%% Set parameters

p = 5; q = 6; r = 2;                % Parameters
J = 19;                             % Number of levels
h=6;                                % p of the Autoregrssive technique
red = 1/r * 1/(1-p/q);              % Redundancy
% Add path to functions
addpath radwt_functions_v2          % use _v1 or _v2
addpath extra_functions

%% Plot frequency responses of RADWT
% subplot(2,1,2)
% Plot_FreqResps(p,q,s,J,fs);
% orient landscape
% print -dpdf figures/demo3_fig1

%% Compute wavelet representation of the first segment of signl1

%%%%%%%%%%%% 1st Segment of MI %%%%%%%%
x=sign11;
N=length(x);
x=x';
[n m]=size(x);
for i=1:n             
  w = radwt(x(i,:), J, p, q, r);
  y=iradwt(w,p,q,r);
  y=real(y(1:N));
  Y(i,:)=y;
end
Y=Y';
%%
%%%%%%%%%%%% 2nd Segment of MI %%%%%%%
x2=sign12;
N2=length(x2);
x2=x2';
[n2 m2]=size(x2);
for i=1:n2             
  w2 = radwt(x2(i,:), J, p, q, r);
  y2=iradwt(w2,p,q,r);
  y2=real(y2(1:N2));
  Y2(i,:)=y2;
end
Y2=Y2';
%%
%%%%%%%%%%%% 3rd Segment of MI %%%%%%
x3=sign13;
N3=length(x3);
x3=x3';
[n3 m3]=size(x3);
for i=1:n3             
  w3 = radwt(x3(i,:), J, p, q, r);
  y3=iradwt(w3,p,q,r);
  y3=real(y3(1:N2));
  Y3(i,:)=y3;
end
Y3=Y3';
%%
%%%%%%%%% 4th Segment of MI %%%%%%%%
x4=sign14;
N4=length(x4);
x4=x4';
[n4 m4]=size(x4);
for i=1:n4             
  w4 = radwt(x4(i,:), J, p, q, r);
  y4=iradwt(w4,p,q,r);
  y4=real(y4(1:N4));
  Y4(i,:)=y4;
end
Y4=Y4';
%% using Autoregrasive
AR1=pcov(Y,h);
AR2=pcov(Y2,h);
AR3=pcov(Y3,h);
AR4=pcov(Y4,h);
AR_RFay=[AR1;AR2;AR3;AR4];
%%
% clear all 
% clc
load('Subject AR_6.mat')
H(1:2580,1)=[1];
F(1:2580,1)=[2];
% 
% RH_aa=[AR_RHaa,H];RF_aa=[AR_RFaa,F];aa=[RH_aa;RF_aa];
% RH_al=[AR_RHal,H];RF_al=[AR_RFal,F];al=[RH_al;RF_al];
% RH_av=[AR_RHav,H];RF_av=[AR_RFav,F];av=[RH_av;RF_av];
% RH_aw=[AR_RHaw,H];RF_aw=[AR_RFaw,F];aw=[RH_aw;RF_aw];
% RH_ay=[AR_RHay,H];RF_ay=[AR_RFay,F];ay=[RH_ay;RF_ay];

RH=[AR_RHaa;AR_RHal;AR_RHav;AR_RHaw;AR_RHay];
RF=[AR_RFaa;AR_RFal;AR_RFav;AR_RFaw;AR_RFay];
aa=[RH,H;RF,F];
% aa=[AR_RHaa;AR_RFaa];
% al=[AR_RHal;AR_RFal];
% av=[AR_RHav;AR_RFav];
% aw=[AR_RHaw;AR_RFaw];
% ay=[AR_RHay;AR_RFay];

%%
X=[aa];    % 6-fold cross validation for Subject aa

a1=[X(1:430,:);X(2581:3010,:)];
a2=[X(431:860,:);X(3011:3440,:)];
a3=[X(861:1290,:);X(3441:3870,:)];
a4=[X(1291:1720,:);X(3871:4300,:)];
a5=[X(1721:2150,:);X(4301:4730,:)];
a6=[X(2151:2580,:);X(4731:5160,:)];
Xtrain=[a1;a2;a5;a6];

% y1=[1*ones(86,1);2*ones(86,1)];
% y2=[1*ones(86,1);2*ones(86,1)];
% y3=[1*ones(86,1);2*ones(86,1)];
% y4=[1*ones(86,1);2*ones(86,1)];
% y5=[1*ones(86,1);2*ones(86,1)];
% y6=[1*ones(86,1);2*ones(86,1)];
% 
% Y=[y1;y2;y3;y5;y6];

Xtest=[a3;a4];
% Xtest1=[a1;a2];%Xtest=testing set
% Xtest2=[a2;a3];
% Xtest3=[a3;a4];
% Xtest4=[a4;a5];
% Xtest5=[a5;a6];
% Xtest6=[a6;a1];
% Ytest=[y4];
%%
%%%%%%%%%%%%% WKNN classifier %%%%%%%%%

inputTable = array2table(Xtest, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89', 'column_90', 'column_91', 'column_92', 'column_93', 'column_94', 'column_95', 'column_96', 'column_97', 'column_98', 'column_99', 'column_100', 'column_101', 'column_102', 'column_103', 'column_104', 'column_105', 'column_106', 'column_107', 'column_108', 'column_109', 'column_110', 'column_111', 'column_112', 'column_113', 'column_114', 'column_115', 'column_116', 'column_117', 'column_118', 'column_119'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89', 'column_90', 'column_91', 'column_92', 'column_93', 'column_94', 'column_95', 'column_96', 'column_97', 'column_98', 'column_99', 'column_100', 'column_101', 'column_102', 'column_103', 'column_104', 'column_105', 'column_106', 'column_107', 'column_108', 'column_109', 'column_110', 'column_111', 'column_112', 'column_113', 'column_114', 'column_115', 'column_116', 'column_117', 'column_118'};
predictors = inputTable(:, predictorNames);
response = inputTable.column_119;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
classificationKNN = fitcknn(...
    predictors, ...
    response, ...
    'Distance', 'Euclidean', ...
    'Exponent', [], ...
    'NumNeighbors', 10, ...
    'DistanceWeight', 'SquaredInverse', ...
    'Standardize', true, ...
    'ClassNames', [1; 2]);

% Create the result struct with predict function
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
knnPredictFcn = @(x) predict(classificationKNN, x);
trainedClassifier.predictFcn = @(x) knnPredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.ClassificationKNN = classificationKNN;
inputTable = array2table(Xtrain, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89', 'column_90', 'column_91', 'column_92', 'column_93', 'column_94', 'column_95', 'column_96', 'column_97', 'column_98', 'column_99', 'column_100', 'column_101', 'column_102', 'column_103', 'column_104', 'column_105', 'column_106', 'column_107', 'column_108', 'column_109', 'column_110', 'column_111', 'column_112', 'column_113', 'column_114', 'column_115', 'column_116', 'column_117', 'column_118', 'column_119'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89', 'column_90', 'column_91', 'column_92', 'column_93', 'column_94', 'column_95', 'column_96', 'column_97', 'column_98', 'column_99', 'column_100', 'column_101', 'column_102', 'column_103', 'column_104', 'column_105', 'column_106', 'column_107', 'column_108', 'column_109', 'column_110', 'column_111', 'column_112', 'column_113', 'column_114', 'column_115', 'column_116', 'column_117', 'column_118'};
predictors = inputTable(:, predictorNames);
response = inputTable.column_119;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];


% Compute resubstitution accuracy
validationAccuracy1 = 1 - resubLoss(trainedClassifier.ClassificationKNN, 'LossFun', 'ClassifError');
ACC_WKNN=validationAccuracy1*100


%%%%%%%%%%% FTree classifier %%%%%%%%%

inputTable2 = array2table(Xtest, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89', 'column_90', 'column_91', 'column_92', 'column_93', 'column_94', 'column_95', 'column_96', 'column_97', 'column_98', 'column_99', 'column_100', 'column_101', 'column_102', 'column_103', 'column_104', 'column_105', 'column_106', 'column_107', 'column_108', 'column_109', 'column_110', 'column_111', 'column_112', 'column_113', 'column_114', 'column_115', 'column_116', 'column_117', 'column_118', 'column_119'});

predictorNames2 = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89', 'column_90', 'column_91', 'column_92', 'column_93', 'column_94', 'column_95', 'column_96', 'column_97', 'column_98', 'column_99', 'column_100', 'column_101', 'column_102', 'column_103', 'column_104', 'column_105', 'column_106', 'column_107', 'column_108', 'column_109', 'column_110', 'column_111', 'column_112', 'column_113', 'column_114', 'column_115', 'column_116', 'column_117', 'column_118'};
predictors2 = inputTable2(:, predictorNames2);
response2 = inputTable2.column_119;
isCategoricalPredictor2 = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
classificationTree = fitctree(...
    predictors2, ...
    response2, ...
    'SplitCriterion', 'gdi', ...
    'MaxNumSplits', 100, ...
    'Surrogate', 'off', ...
    'ClassNames', [1; 2]);

% Create the result struct with predict function
predictorExtractionFcn1 = @(x) array2table(x, 'VariableNames', predictorNames2);
treePredictFcn = @(x) predict(classificationTree, x);
trainedClassifier2.predictFcn = @(x) treePredictFcn(predictorExtractionFcn1(x));

% Add additional fields to the result struct
trainedClassifier2.ClassificationTree = classificationTree;
inputTable2 = array2table(Xtrain, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89', 'column_90', 'column_91', 'column_92', 'column_93', 'column_94', 'column_95', 'column_96', 'column_97', 'column_98', 'column_99', 'column_100', 'column_101', 'column_102', 'column_103', 'column_104', 'column_105', 'column_106', 'column_107', 'column_108', 'column_109', 'column_110', 'column_111', 'column_112', 'column_113', 'column_114', 'column_115', 'column_116', 'column_117', 'column_118', 'column_119'});

predictorNames2 = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89', 'column_90', 'column_91', 'column_92', 'column_93', 'column_94', 'column_95', 'column_96', 'column_97', 'column_98', 'column_99', 'column_100', 'column_101', 'column_102', 'column_103', 'column_104', 'column_105', 'column_106', 'column_107', 'column_108', 'column_109', 'column_110', 'column_111', 'column_112', 'column_113', 'column_114', 'column_115', 'column_116', 'column_117', 'column_118'};
predictors2 = inputTable2(:, predictorNames2);
response2 = inputTable2.column_119;
isCategoricalPredictor2 = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

% Compute resubstitution accuracy
validationAccuracy2 = 1 - resubLoss(trainedClassifier2.ClassificationTree, 'LossFun', 'ClassifError');
ACC_FTree=validationAccuracy2*100


%%%%%%%%%%%% BST classifier %%%%%%%

inputTable3 = array2table(Xtest, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89', 'column_90', 'column_91', 'column_92', 'column_93', 'column_94', 'column_95', 'column_96', 'column_97', 'column_98', 'column_99', 'column_100', 'column_101', 'column_102', 'column_103', 'column_104', 'column_105', 'column_106', 'column_107', 'column_108', 'column_109', 'column_110', 'column_111', 'column_112', 'column_113', 'column_114', 'column_115', 'column_116', 'column_117', 'column_118', 'column_119'});

predictorNames3 = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89', 'column_90', 'column_91', 'column_92', 'column_93', 'column_94', 'column_95', 'column_96', 'column_97', 'column_98', 'column_99', 'column_100', 'column_101', 'column_102', 'column_103', 'column_104', 'column_105', 'column_106', 'column_107', 'column_108', 'column_109', 'column_110', 'column_111', 'column_112', 'column_113', 'column_114', 'column_115', 'column_116', 'column_117', 'column_118'};
predictors3 = inputTable3(:, predictorNames3);
response3 = inputTable3.column_119;
isCategoricalPredictor3 = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
template = templateTree(...
    'MaxNumSplits', 20);
classificationEnsemble = fitcensemble(...
    predictors3, ...
    response3, ...
    'Method', 'AdaBoostM1', ...
    'NumLearningCycles', 30, ...
    'Learners', template, ...
    'LearnRate', 0.1, ...
    'ClassNames', [1; 2]);

% Create the result struct with predict function
predictorExtractionFcn2 = @(x) array2table(x, 'VariableNames', predictorNames);
ensemblePredictFcn = @(x) predict(classificationEnsemble, x);
trainedClassifier3.predictFcn = @(x) ensemblePredictFcn(predictorExtractionFcn2(x));

% Add additional fields to the result struct
trainedClassifier3.ClassificationEnsemble = classificationEnsemble;
inputTable3 = array2table(Xtrain, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89', 'column_90', 'column_91', 'column_92', 'column_93', 'column_94', 'column_95', 'column_96', 'column_97', 'column_98', 'column_99', 'column_100', 'column_101', 'column_102', 'column_103', 'column_104', 'column_105', 'column_106', 'column_107', 'column_108', 'column_109', 'column_110', 'column_111', 'column_112', 'column_113', 'column_114', 'column_115', 'column_116', 'column_117', 'column_118', 'column_119'});

predictorNames3 = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89', 'column_90', 'column_91', 'column_92', 'column_93', 'column_94', 'column_95', 'column_96', 'column_97', 'column_98', 'column_99', 'column_100', 'column_101', 'column_102', 'column_103', 'column_104', 'column_105', 'column_106', 'column_107', 'column_108', 'column_109', 'column_110', 'column_111', 'column_112', 'column_113', 'column_114', 'column_115', 'column_116', 'column_117', 'column_118'};
predictors3 = inputTable3(:, predictorNames3);
response3 = inputTable3.column_119;
isCategoricalPredictor3 = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

% Compute resubstitution accuracy
validationAccuracy3 = 1 - resubLoss(trainedClassifier3.ClassificationEnsemble, 'LossFun', 'ClassifError');
ACC_BST=validationAccuracy3*100

