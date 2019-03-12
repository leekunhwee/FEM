%%%%%%%%%%%%%%%%%%%%%%
%Copyright Jianhui Li%
%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defining system parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NSpokes = 8; % number of spokes
NJoints = NSpokes + 1; % number of joints
NLinks = 2*NSpokes; % spokes _ edge links
R0 = 1; % circle radius
Dspoke = 0.01; % diameter of spoke , m
Douter = 0.01; % diameter of outer links, m
E = 210e9; %Pa, mild steel Young's modulus
sigma_Yield = 210e6; %Pa, mild steel ultimate strength
xA = 0.161 / 2; % x of the centre, last 3 digits of the student number divided by 2
yA = 0.979 / 2; % y of the centre, 6, 5, 4th from last digits of the student number divided by 2
Fapplied = 1362; % N, force applied to stress the system

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating parameters of individual links and joints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Spoke_phi = linspace(0,2*pi()*(1-1/NSpokes),NSpokes) + pi()/2; % angle of spokes' fixing locations
Joint_x0 = [R0 * cos(Spoke_phi) xA]'; % x coordinates of Joints, N outer and 1 central
Joint_y0 = [R0 * sin(Spoke_phi) yA]'; % y coordinates of Joints, N outer and 1 central

% Now, let the links 1...N are spokes, links (N+1)...2*N are the outer connections
Link_i = [1:NSpokes 1:NSpokes]'; % spokes start at outer points, connecting links start where spokes start
Link_j = [ones(1,NSpokes)*NJoints 2:NSpokes 1]';
% spokes end at the axis, outer links end at the beginning of the next spoke

Link_theta = atan2(Joint_y0(Link_j) - Joint_y0(Link_i), Joint_x0(Link_j) - Joint_x0(Link_i));
% spokes "begin" at the system edge and "end" in the centre

Link_la = cos(Link_theta); % lambda in stiffness matrix of a member
Link_mu = sin(Link_theta); % mu in stiffness matrix of a member
Link_L0 = sqrt( (Joint_x0(Link_j) - Joint_x0(Link_i)).^2 + (Joint_y0(Link_j) - Joint_y0(Link_i)).^2 );
% % length of links

Spoke_D = ones(NSpokes,1)*Dspoke; % spoke diameter, m
Outer_D = ones(NSpokes,1)*Douter; % outer link diameter, m

Link_A = pi()/4*[Spoke_D.^2 ; Outer_D.^2]; % cross-section, spokes then links
Link_Stf = Link_A*E./Link_L0;% element stifness, spokes then links

Link_I = pi()/64*[Spoke_D.^4 ; Outer_D.^4]; % second moment for buckling
Link_Fmin = - pi()^2 * E * Link_I / 4 ./ Link_L0.^2; % buckling on compression
Link_Fmax = Link_A * sigma_Yield; % breaking on tension

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modify the matrix by artificially setting its stiffness to zero (Spoke 4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Link_Stf(4)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Composing stiffness matrix of the system and force vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Joint_Force = zeros(2*NJoints,1); 
Joint_Force(10) = -Fapplied;
% force vector, negative vertical to the (NSpokes/2 + 1)st spoke
STIFFNESS = zeros(2*NJoints,2*NJoints); % stiffness matrix
for kLink = 1: NLinks % populating the global stiffness matrix
    i0 = (Link_i(kLink)-1)*2; % row after which the node i begins
    j0 = (Link_j(kLink)-1)*2; % row after which the node j begins
    la = Link_la(kLink); % cos(theta) for the current spoke
    mu = Link_mu(kLink); % sin(theta) for the current spoke
    mtrx = Link_Stf(kLink)*[la^2 la*mu -la^2 -la*mu;
                            la*mu mu^2 -la*mu -mu^2;
                            -la^2 -la*mu la^2 la*mu;
                            -la*mu -mu^2 la*mu mu^2];
    % 2 X 2 part of the stiffness matrix for the current link, see Megson (17.23) which contains an error
    
    STIFFNESS((i0+1):(i0+2),(i0+1):(i0+2)) = STIFFNESS((i0+1):(i0+2),(i0+1):(i0+2)) + mtrx(1:2,1:2);
    STIFFNESS((i0+1):(i0+2),(j0+1):(j0+2)) = STIFFNESS((i0+1):(i0+2),(j0+1):(j0+2)) + mtrx(1:2,3:4);
    STIFFNESS((j0+1):(j0+2),(i0+1):(i0+2)) = STIFFNESS((j0+1):(j0+2),(i0+1):(i0+2)) + mtrx(3:4,1:2);
    STIFFNESS((j0+1):(j0+2),(j0+1):(j0+2)) = STIFFNESS((j0+1):(j0+2),(j0+1):(j0+2)) + mtrx(3:4,3:4);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating displacements from known forces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STIFFNESS_cut = STIFFNESS(5:end,5:end); % removing rows and columns corresponding to zero displacements
F_cut = Joint_Force(5:end); % removing forces applied to joints with zero displacements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_cut = STIFFNESS_cut \ F_cut; % CORE OPERATION: SOLVING TO FIND DEFORMATIONS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = [0;0;0;0;w_cut]; % completing deformation vector with zeroes
Joint_Force = STIFFNESS*w; % calculating force including fixed joints (use w and STIFFNESS)
Joint_x = Joint_x0 + w(1:2:end); % x of joints with deformations
Joint_y = Joint_y0 + w(2:2:end); % y of joints with deformations
Link_L = sqrt( (Joint_x(Link_j) - Joint_x(Link_i)).^2 + (Joint_y(Link_j) - Joint_y(Link_i)).^2 );
% new length of links
Link_F = (Link_L-Link_L0).*Link_Stf; % force in the link

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WSTRETCH = 2000; % stretch the shown deformations to make them visible
Joint_xD = Joint_x0 + w(1:2:end)*WSTRETCH; % x with deformations stretched WSTRETCH times
Joint_yD = Joint_y0 + w(2:2:end)*WSTRETCH; % y with deformations stretched WSTRETCH times
figure(33);
clf;
plot(Joint_x0(1),Joint_y0(1),'b.','markersize',30); % fixed joint
hold on; axis equal; axis(R0*1.3*[-1 1 -1 1]);
plot(Joint_x0(2),Joint_y0(2),'b.','markersize',30); % fixed joint
plot(Joint_x0(NSpokes/2+1),Joint_y0(NSpokes/2+1),'b.','markersize',30); % loaded joint

for kLink = 1:NLinks % plotting links
    if Link_Stf(kLink)>1e-20 % only plotting existing links
    clr = [(max((Link_F))-(Link_F(kLink))) (Link_F(kLink)-min(Link_F)) 0]/(max((Link_F))-min((Link_F)));
    
    % colour coding: tension GREENer, compression REDer
    if (Link_F(kLink) > Link_Fmax(kLink)) || (Link_F(kLink) < Link_Fmin(kLink)) 
        lntype = '--'; 
    else
        lntype = '-';
    end % dashed line if the link fails
    
    % Undeformation
    plot([Joint_x0(Link_i(kLink)) Joint_x0(Link_j(kLink))], ...
    [Joint_y0(Link_i(kLink)) Joint_y0(Link_j(kLink))],'-','color',[0.7 0.7 0.7]);
    
    % Deformation
    plot([Joint_xD(Link_i(kLink)) Joint_xD(Link_j(kLink))], ...
    [Joint_yD(Link_i(kLink)) Joint_yD(Link_j(kLink))],lntype,'linewidth',2,'color',clr);
    
    text((Joint_x0(Link_i(kLink))+Joint_x0(Link_j(kLink)))/2-0.05, ...
    (Joint_y0(Link_i(kLink))+Joint_y0(Link_j(kLink)))/2,num2str(kLink),'fontsize',16);
    % captions of links
    end
end

plot(Joint_xD,Joint_yD,'k.','markersize',15);
for kJoint = 1:NJoints % signing joint numbers
    text(Joint_x0(kJoint)*1.12-0.05,Joint_y0(kJoint)*1.1,num2str(kJoint),'color','b','fontsize',20);
end
plot(Joint_xD(NSpokes/2+1),Joint_yD(NSpokes/2+1),'r.','markersize',30);
xlabel('x'); ylabel('y');
set(gca,'fontsize',16,'linewidth',1);