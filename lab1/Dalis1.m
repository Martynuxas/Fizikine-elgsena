function Dalis1
clc, close all, clear all

% ---- construction data ----
% mass constants
m1 = 1; m2 = 1.5; m3 = 2;
% stiffness of elements
k1 = 4000;        
% forces
f1=400; f2=800; 
% time moments when forces are added to the construction
tf1 = 0.05; tf2 = 0.1;

%  --- information about nodes ---
% array of node constraints
IS = [  0,   1,   0,  0,  0,   0]; 
 % displacements of node nr. 2
u2_deltaU  =  0.5;  
 % time moment, when kinematic condition starts
u2_t_start =  0.2;    
% how long the motion takes
u2_deltaT  =  0.1;    
% masses of the particles
m=   [  m1, 1, m1, m1, m2, m3]; 
% forces applied to the particles
F =  [  -f1,   0, -f2,  0,  0,  0];   
% time moments when the forces start impacting the structure
tf = [  tf1,   0, tf2,  0,  0, 0];

%  --- information about elements ---
% stiffness coefficients
k=[k1, k1, k1, k1, k1];             
 % damping coefficients
c=[10, 10, 10, 10, 10];  

% information about elements (node1, node2, visualization level in y axis)
ind=[1, 3, 3;   
     2, 4, 1;  
     3, 5, 3;
     4, 5, 1;
     5, 6, 1];
 
%  --- visualization data ---
 
% Coordinates of particles 
% x axis 2-4-5-6-3
x=[4, 1, 2, 2, 3, 4];   
% y axis
y=[3, 1, 3, 1, 2, 1];   
% Diameters of particles 
rad1=0.2; rad2=1.5;     
rads = [rad1, rad1; 
        rad1, rad1;
        rad1, rad1;
        rad1, rad1;
        rad1, rad2;
        rad1, rad1];

nmz=length(m); % total number of nodes (in this example = 6)
nel=length(k); % total number of elements (in this example = 5)
% initial displacements and velocities
U=zeros(nmz,1);  % displacements
DU=zeros(nmz,1); % velocities

% initial data visualization
visualization(x, y, ind, U,  rads, 0);


% numerical integration  - dynamic modelling
Urez=zeros(nmz,1);      % array to save the data of nodes displacements
TT = 0.5;     % numerical integration time (dynamic modelling time)
dt=0.001;     % numerical integration step
for t=0:dt:TT       % loop of numerical integration
    % updating acceleration
    DDU=acceleration(m,k,c,ind,F,tf, U,DU,t,@F_time_function); 
    % updating velocities
    DU=DU+dt*DDU;   
    % boundary condition of velocities
    DU(find(IS))=0; 
    DU(2) = du_time_function(u2_t_start, u2_deltaT, u2_deltaU, t);
    % updating displacements
    U=U+dt*DU;      
    % saving intermediate results
    jj=round(t/dt)+1; 
    Urez(:,jj)=U;  
    % animation condition of the structure
    figure(1);   visualization(x, y, ind, U,  rads,t);
    pause(0.01);
 end
 
 
% representing displacements of particles in time
figure(2);hold on; grid on;
color={'b-';'r-';'g-';'m-';'c-';'k-';};
% plotting displacements in time of each node
for i=1:nmz  
    plot([0:dt:TT],Urez(i,:),color{i});
end

plot([tf1, tf1], [-0.2, 0.2], 'g--'); 
plot([tf1, tf1], [-0.2, 0.2], 'k--');
plot([u2_t_start, u2_t_start], [-0.2, 0.2], 'r--');
legend('1 node','2 node', '3 node', '4 node', '5 node', '6 node');
xlabel('Time (s)');
ylabel('Displacements (m)');
end
% *********************************************
% function that governs application of forces
function Ft=F_time_function(F, tf, t)
    nmz=length(F);
    Ft = zeros(nmz, 1);
    for i=1:nmz,
        if(tf(i) < t) 
            Ft(i) = F(i);
        end
    end
end
% *********************************************
% time function of velocities
function du=du_time_function(t_start, deltaT, deltaU, t)
    if t >= t_start  % if addition of velocity started
        if  t <= t_start + deltaT, % if the node is moving
            % velocity is prescribed
            b = deltaT; a = deltaU;
            
            
            du = (a*pi*cos((3*pi)/2 + (pi*t)/b))/(2*b);
            
            
            if(t == t_start || t == t_start + deltaT),
               du = du/2; 
            end
        else,
            du = 0; % displacement has already been applied, velocity of particle = 0
        end
    else
        du = 0; % motion has not started, velocity of particle = 0
    end
    return
end

% *********************************************
function DDU=acceleration(m,k,c,ind,F,tf,U,DU,t,Ft)
    nel=length(k);nmz=length(m);
    DDU=zeros(nmz,1);
    DDU= Ft(F, tf, t); % adding external forces
    for iii=1:nel      % adding forces of elements
        i=ind(iii,1);j=ind(iii,2);
        T=(U(j)-U(i))*k(iii)+(DU(j)-DU(i))*c(iii);
        DDU(i)=DDU(i)+T;
        DDU(j)=DDU(j)-T;
    end
    % now variable DDU represents all forces (F) in the structure
    % solving F = m * a, and getting acceleration
    DDU=DDU./m';
return
end

% *********************************************
function visualization(x,y,ind,U,rads,t)
    % prepare figure for visualization
    clf; hold on; grid on;axis equal; 
    axis([min(x)-1 ,max(x)+1,min(y)-1 ,max(y)+1]);
    title(['t =', num2str(t), '(s)']);
    xlabel('X');
    ylabel('Y');
    % calculate number of nodes and elements
    nmz=length(x);nel=size(ind,1);     
    % visualizing elements
    for i=1:nel
        line([x(ind(i,1))+U(ind(i,1)), ...
            x(ind(i,2))+U(ind(i,2)) ], ...
            [ind(i,3), ind(i,3)],...
            'Color','blue','Linewidth',2);
    end
    % visualizing nodes
    for i=1:nmz
        rectangle('Position',[x(i)+U(i)-rads(i,1), ...
            y(i)-rads(i,2),rads(i,1)*2,rads(i,2)*2]....
            ,'Curvature',[1,1],'FaceColor',[0.4 0.6 1]);
    end
    disp('po = ');disp(U);
return
end


