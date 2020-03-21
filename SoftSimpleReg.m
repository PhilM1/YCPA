%2.b)This script simulates the E and H fields as it propagates through a
%medium.

%2.c.i) When commenting out the inclusion box, the wave simply propagates
%without reflecting

%2.c.ii) The bc structure is used to store the plane wave into an object
%so that it can be propagated with the PlaneWaveBC.m function.

%2.c.iii) The bc{1}.s(1) is setting up the wave parameters such as its
%initial location, type of wave and function to be used to propagate that
%wave.

%2.c.iv) The xm, xp, ym, yp are the adjacent cells which are x-minus,
%x-plus, y-minus, y-plus. By observation, changing the type from 'a' to 'e'
%it looks like the waves are scattering much more and staying around a lot
%longer

%3.b) Adding a grating made is so that the wave gets many more reflections
%from each of the gratings and creates a "trapping" effect between two
%gratings as well where the wave oscillates in between each of the
%gratings. It also takes more time for the wave to propagate to the end of
%the gratings.

%3.c) By lowering the frequency of the excitation it seems like the
%gratings have little effect on the wave propagation. When increasing the
%frequency, the waves have a very difficult time even getting through the
%first grating and almost all of the wave gets bounced off.

%4.a) changed the grating to a funnel type object
%4.b) Added another wave with slight displacement


winstyle = 'docked';
% winstyle = 'normal';

set(0,'DefaultFigureWindowStyle',winstyle)
set(0,'defaultaxesfontsize',18)
set(0,'defaultaxesfontname','Times New Roman')
% set(0,'defaultfigurecolor',[1 1 1])

% clear VARIABLES;
clear
global spatialFactor;
global c_eps_0 c_mu_0 c_c c_eta_0
global simulationStopTimes;
global AsymForcing
global dels
global SurfHxLeft SurfHyLeft SurfEzLeft SurfHxRight SurfHyRight SurfEzRight



dels = 0.75;
spatialFactor = 1;

c_c = 299792458;                  % speed of light
c_eps_0 = 8.8542149e-12;          % vacuum permittivity
c_mu_0 = 1.2566370614e-6;         % vacuum permeability
c_eta_0 = sqrt(c_mu_0/c_eps_0);


tSim = 200e-15; %simulation time
f = 230e12; %frequency
%f = 230e11; %frequency for grating testing (low)
%f = 430e12; %frequency for grating testing (high)
lambda = c_c/f; %wavelength

xMax{1} = 20e-6; %box or window size in x
nx{1} = 200; % number of divisions in x (matrix size)
ny{1} = 0.75*nx{1}; % number of divisions in y (matrix size)
 

Reg.n = 1;

mu{1} = ones(nx{1},ny{1})*c_mu_0; %permeability matrix

epi{1} = ones(nx{1},ny{1})*c_eps_0; %permittivity matrix

%epi{1}(125:150,55:95)= c_eps_0*11.3; % add box with different permitivity

%Create a grating:
% epi{1}(100:101,55:95)= c_eps_0*11.3;
% epi{1}(110:111,55:95)= c_eps_0*11.3;
% epi{1}(120:121,55:95)= c_eps_0*11.3;
% epi{1}(130:131,55:95)= c_eps_0*11.3;
% epi{1}(140:141,55:95)= c_eps_0*11.3;
% epi{1}(150:151,55:95)= c_eps_0*11.3;



%funnel structure
funnel_Length = 50;
funnel_End = 150;
funnel_ExitSize = 10;
funnel_BandThickness = 5;
funnel_Center = 75;

for i = 1:funnel_Length
    epi{1}(funnel_End+1-i, funnel_Center+(i-1)+funnel_ExitSize/2:funnel_Center+(i-1)+funnel_ExitSize/2+funnel_BandThickness)= c_eps_0*11.3; %top
    epi{1}(funnel_End+1-i, funnel_Center-(i-1)-funnel_ExitSize/2-funnel_BandThickness:funnel_Center-(i-1)-funnel_ExitSize/2)= c_eps_0*11.3;%bottom
end

sigma{1} = zeros(nx{1},ny{1});
sigmaH{1} = zeros(nx{1},ny{1});

dx = xMax{1}/nx{1}; %division size
dt = 0.25*dx/c_c; %timestep size
nSteps = round(tSim/dt*2); %number of steps for simulation
yMax = ny{1}*dx; %box or window size in y
nsteps_lamda = lambda/dx; 

movie = 1;
Plot.off = 0;
Plot.pl = 0;
Plot.ori = '13';
Plot.N = 100;
Plot.MaxEz = 1.1;
Plot.MaxH = Plot.MaxEz/c_eta_0;
Plot.pv = [0 0 90];
Plot.reglim = [0 xMax{1} 0 yMax];


bc{1}.NumS = 2;
bc{1}.s(1).xpos = nx{1}/(10) + 1;
bc{1}.s(1).type = 'ss';
bc{1}.s(1).fct = @PlaneWaveBC;

%second wave in front of first
bc{1}.s(2).xpos = nx{1}/(4) + 1;
bc{1}.s(2).type = 'e';
bc{1}.s(2).fct = @PlaneWaveBC;

% mag = -1/c_eta_0;
mag = 1;
phi = 0;
omega = f*2*pi;
betap = 0;
t0 = 30e-15;
st = 15e-15;
%st = -0.05; %for grating example
s = 0;
y0 = yMax/2;
sty = 1.5*lambda;
bc{1}.s(1).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'s'};

%second wave
bc{1}.s(2).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'s'};

Plot.y0 = round(y0/dx);

bc{1}.xm.type = 'a';
bc{1}.xp.type = 'a';
bc{1}.ym.type = 'a';
bc{1}.yp.type = 'a';

pml.width = 20 * spatialFactor;
pml.m = 3.5;

Reg.n  = 1;
Reg.xoff{1} = 0;
Reg.yoff{1} = 0;

RunYeeReg






