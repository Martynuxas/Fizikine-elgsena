% 1 LD --------------------------------------------------------------------
function LD1
    close all; clf; clc; clear all;
  
    % ---- construction data ----
    % mass constants
    m1 = 2; m2 = 1;     
    % stiffness constrants
    k1 = 100; k2 = 2500;  
    % spring dampers
    c1 =10; c2 = 20; 
    % masses of particles   
    m = [m1;  m1;  m1; m1; m2; m1; m2]; 
    % coordinates of particles
    cords = [  1,  1;  3,  1;  5,  1;  5,  2;  3,  2; 1, 2; 3, 4 ];    
    % constraints array
    IS = logical([  0,   0,   1,   1,   0,   0,   1,   0,   0,   0, 1, 0, 1, 1]);
    % particles displacements start times
    U_t_start = [   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 0, 0,  0.1, 0.1];
    % particle displacements
    deltaU    = [   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 0, 0,    1,   1]; 
    % particles kinematics times
    U_deltaT = [    0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 0, 0,    2,   2];
    % elements (springs) array
    elm = [1, 2; 1 , 6;  2, 3; 2 , 5; 3 , 4;  4 , 5; 4 , 7;  5, 6;  5, 7; 6, 7];
    % elements (springs) stiffness
    k = [k2, k2 ,k2, k2 , k2, k1, k2, k1, k2, k2];
    % elements (springs) dampings
    c = [c2, c2 ,c2, c2 , c2, c1, c2, c1, c2, c2];    
    
    % ading m*g forces to the particles
    g = -9.8; 
    nNodes = length(m);
    F = zeros(nNodes *2,1);
    for(i = 1:nNodes)
        F(i*2) = m(i) *g;
    end
    
    % total number of particles in construction
    nNodes = length(m);
    % degree of freadom in the particles
    DOF = 2; 
    % displacements array
    U  = zeros(nNodes*DOF,1); % displacement vector
    % rendering construction data
    rendering(U,elm,cords,F,IS,1, 'Initial construction');
    pause(0.1);

    
    % numerical integration  - dynamic modelling
    TT=20; dt=0.01;   % intrgration time and step
    U  = zeros(nNodes*DOF,1); % displacement vector
    DU = zeros(nNodes*DOF,1); % velocities

    for t=0:dt:TT
          % updating acceleration
        DDU=acceleration(U,DU,t,m, F,elm,cords, k, c);
        % updating velocities
        DU=DU+dt*DDU'; 
        % boundary condition of velocities
        DU(IS)=0; 
        % updating velocities of particles of kinematic motion
        for i=1:nNodes*DOF
           if IS(i) == 1 
                [ xxx DU(i) xxx] = ...
                    du_time_function(U_t_start(i), deltaU(i), U_deltaT(i), t);
           end 
        end
        % updating displacements
        U=U+dt*DU;      
        % rendering construction
        rendering(U,elm,cords,F,IS, 3, strcat('CMS:  ', num2str(t), '(s)'));
        pause(0.01);
    end
    disp('po CMS: U = ');disp(U');
    rendering(U,elm,cords,F,IS, 3, strcat('CMS:  ', num2str(t), '(s)'));
end



% acceleration function
function DDU = acceleration(U,DU,t,mass, F,elm,cords, k, c)
    dof=2; 
    nmz=length(mass);NN=nmz*dof;
    siz=size(elm);nel=siz(1);
    T=F;  % adding external forces
    % assembling forces of elements to the global forces vectors
    for i=1:nel   % 
        r=elm(i,1);s=elm(i,2);  % nubers of nodes of spring ends
        ur=U((r-1)*dof+1);vr=U(r*dof);     % displacements of nodes of spring ends
        us=U((s-1)*dof+1);vs=U(s*dof);
        xr=cords(r,1)+ur;yr=cords(r,2)+vr;
        xs=cords(s,1)+us;ys=cords(s,2)+vs;
        dur=DU((r-1)*dof+1);dvr=DU(r*dof); % velocities of spring particles
        dus=DU((s-1)*dof+1);dvs=DU(s*dof);

        l0=sqrt((cords(s,1)-cords(r,1))^2+(cords(s,2)-cords(r,2))^2); % initial spring length
        lrs=sqrt((xs-xr)^2+(ys-yr)^2);  % current lenght of spring
        n= [xs-xr , ys-yr]/lrs;         % element forcwe normal vector

        Trs=k(i)*(lrs-l0)+c(i)*dot( n, [dus-dur,dvs-dvr ] ); % force created by element

        % adding spring forces to the global node forces array
        T((r-1)*dof+1)=T((r-1)*dof+1)+Trs*n(1); 
        T(r*dof)=T(r*dof)+Trs*n(2);
        T((s-1)*dof+1)=T((s-1)*dof+1)-Trs*n(1);
        T(s*dof)=T(s*dof)-Trs*n(2);

    end
DDU(1:dof:NN)=T(1:dof:end)./mass;   % forces dividing from masses
DDU(2:dof:NN)=T(2:dof:end)./mass;
end

% kinematic function
function [ddu du u]=du_time_function(U_t_start, deltaU,U_deltaT, t)
    % U_t_start(i), deltaU(i), U_deltaT(i)
    ddu =0; du =0;u =0;
    if t >= U_t_start  % if displacement starts
        if  t <= U_t_start + U_deltaT, % if displacement in progress
            if(U_deltaT > 0)
                 % updating acceleration, velocities and displacements
                 omega=(pi)/U_deltaT;
                 ddu = (deltaU /2) * omega^2 * sin (omega *(t - U_t_start)+(3/2)*pi);
                 du = (deltaU /2) * omega * cos(omega *(t - U_t_start)+(3/2)*pi);
                 u = (deltaU /2) * (sin(omega *(t - U_t_start)+(3/2)*pi)+1);
            else,         
                % node displacement = 0 = 0
                ddu =0;du =0;u =0;
            end
        else,
            % displacements already done
            ddu = 0; du = 0; u = deltaU;
        end
    end
    return
end

% rendering function

function rendering(U,elm,cords,F,IS, fig, strTitle)
    DOF=2;
    xx = size(cords); nNodes = xx(1);;
    xx = size(elm); nElm = xx(1);
    NN = nNodes * DOF;
    ff = figure(fig);
    clf(ff);
    axis([0 7 0 6]);
    hold on; grid on;
    title(strTitle);
    
    % geting the construction dimention, to visualize the forces as arrows
    % in the construction
    xlim=get(gca,'XLim'); ylim=get(gca,'YLim');
    xn=xlim(2)-xlim(1);yn=ylim(2)-ylim(1); % axis diapazons
    range=min(xn,yn);           % geting minimum of axis
    maxForce=max(abs(F));       % maximum foce
    % scaling coeficients
    mast= range/maxForce*0.1;   
    constrLength=range/17;      

    for i=1:nNodes
        % rendering particles
        u=U((i-1)*DOF+1);v=U(i*DOF); % displacements of i-th particle
        r=0.2;  % i-os daleles spindulys
        % plotting particle:
        rectangle('Position',[cords(i,1)+u-r,cords(i,2)+v-r,2*r,2*r],'Curvature',[1,1],'FaceColor',[0.4 0.6 1]);
        % ploting forces of exact particle
        fx=F((i-1)*DOF+1)*mast;fy=F(i*DOF)*mast; % length of ith force particle arrow 
        x1=cords(i,1)+u;x2=cords(i,1)+u+fx;y1=cords(i,2)+v;y2=cords(i,2)+v+fy;
        line([x1,x2],[y1,y2],'Color','red','LineWidth',1);
        varr=[x1-x2;y1-y2]; varr =varr/norm(varr)*range/40; % ploting arrow end
        alf=pi/6; transf = [cos(alf) sin(alf);-sin(alf) cos(alf)];
        varr1=transf*varr; line([x2, x2+varr1(1)],[y2, y2+varr1(2)],'Color','red','LineWidth',1);
        varr1=transf'*varr;line([x2, x2+varr1(1)],[y2, y2+varr1(2)],'Color','red','LineWidth',1);
        % constraints of velocities :
        ix=IS((i-1)*DOF+1)*mast;iy=IS(i*DOF)*mast; 
        if ix ~= 0, line(([cords(i,1)+u, cords(i,1)+u]),([cords(i,2)+v-constrLength/2, cords(i,2)+v+constrLength/2]),'Color',[ 0.2 0.2 0.2],'LineWidth',3);end  
        if iy ~= 0, line(([cords(i,1)+u-constrLength/2, cords(i,1)+u+constrLength/2]),([cords(i,2)+v, cords(i,2)+v]),'Color',[ 0.2 0.2 0.2],'LineWidth',3);end  

    end
    % rendering springs
    for i=1:nElm
        r=elm(i,1);s=elm(i,2); % numbers of spring nodes (particles)
        ur=U((r-1)*DOF+1);vr=U(r*DOF);  % displacements of spring particles
        us=U((s-1)*DOF+1);vs=U(s*DOF);
        xr=cords(r,1)+ur;yr=cords(r,2)+vr;
        xs=cords(s,1)+us;ys=cords(s,2)+vs;
        plot([xr,xs] , [yr,ys],'b-');     % spring rendering as segment
    end
return
end