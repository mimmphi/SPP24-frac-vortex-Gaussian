%put hole into the phase and reconstruct

clear
close all

%% INITIALIZE
c = 299792458; %speed of light
la = 633e-9; %wavelength
k = 2*pi/la; %wavenumber
L = .008; %side length(m)%.004 dapat dx=10.8.*(1e-6)
W = L; %width
M = 2^10; %number of samples 1024, 2^10 orig
N = M;
dx = L/M; %src sample interval
dy = W/N;
x = -L/2:dx:L/2-dx; %src coords
y = -W/2:dy:W/2-dy; %src coords
[X,Y]= meshgrid(x,y);
Z=X+1i*Y;


%% MAKE LAGUERRE-GAUSS BEAM
phi = angle(Z);
wo=.0005; %wo=.002, 
l = 3;

rho = sqrt((X.^2)+(Y.^2)); %gaussian
    
circ = exp(-rho.^2./(wo.^2));     % gaussian "aperture" overlaid on CGH
phase = l.*phi+pi;
u = double(circ).*(rho.^abs(l)).*exp(1i*phase);
int_u = abs(u.*u);

%% putting holes on the phase and propagating
sz_hol = size(phase,1)*size(phase,2);
R = randperm(sz_hol);

P_list = 1:-0.05:0; %percentage of information retained
% 1 = no loss; 0 = all loss
% P1 = 1 - P_list;
Q_list = zeros(1,numel(P_list));
T_list = zeros(1,numel(P_list));
F_list = zeros(1,numel(P_list));

%original beam
z = 1.4;
[u1_prop]=propTF(u, L, la, z);
int_uo = abs(u1_prop.*u1_prop);

figure(1)
imagesc(int_uo);
colormap(hot)


for i = 1:numel(P_list)

    P = 1 - P_list(i); %percentage of information loss
    R_indx = R(1:round(P*sz_hol));
    
    for j = 1:numel(R_indx)
        phase(R_indx(j)) = 0;
    end
    
    u1 = double(circ).*(rho.^abs(l)).*exp(1i*phase);
    %u1 = double(circ).*exp(1i*phase);
    % z = L*dx/la;
    [u1_prop]=propTF(u1, L, la, z);
    
    % u1_fft = fftshift(fft2(u1));
    % figure(2)
    % imshow(phase);
    % % exportgraphics(gcf,'phase_l=3.gif','Append',true);

    int_u1 = abs(u1_prop.*u1_prop);
    figure(3)
    imagesc(int_u1);
    colormap(hot)
    % exportgraphics(gcf,'l=3.gif','Append',true);
    
    [Q,T,F] = linfoots(int_uo,int_u1);
    Q_list(i) = Q;
    T_list(i) = T;
    F_list(i) = F;

end

figure(4);
plot(P_list,Q_list,'LineWidth',2);
hold on
plot(P_list,T_list,'LineWidth',2);
plot(P_list,F_list,'LineWidth',2);
ax2 = gca; % current axes
ax2.FontSize = 20;
xlabel('Fraction of the beam')
ylabel('Linfoots Criteria');
legend('Correlation Quality','Structural Content','Fidelity');
hold off






