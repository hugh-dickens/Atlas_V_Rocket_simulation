%10/11/2018
function Atlas_V_Rocket_Simulation
%Simulation assumes zero-yaw condition throughout

%Atlas booster initial conditions
A_booster_thrust_start = 3827000;   %thrust (N)
A_Isp = 311.3;                      %specific impulse/ISP (s)
A_booster_propellant = 284089;      %mass of propellant (kg)
A_booster_inert_mass = 21351;       %inert mass (kg)
    
%Solid rocket booster initial conditions
SR_booster_thrust_start = 1688400;  %thrust of each individual booster(N)
SR_booster_mass = 46697;            %total mass of each booster(kg)
SR_booster_inert_mass = 4000;       %inert mass of each booster (kg)
SR_Isp = 279.3;                     %Isp (s)
SR_boosters = 3;                    %number of solid rocket boosters
j_time = 115;                       %time of jettison/burn time of SRBs (s)

%Other mass of rocket initial conditions
payload = 0;                        %mass of payload (kg)
r_mass = 28000;                     %mass of common centaur and propellant, payload fairing and adapters (kg)

%Initial system parameters
dt = 0.01;                          %interval step (s)
t_end = 250;                        %simulation end/jettison time of Atlas booster (s)

%Initial rocket gimbal conditions
tg1 = 12;                           %start time of rocket gimbaling due to nozzle angle changing
tg2 = 14;                           %end time of rocket gimbaling
angle_gimbal = 3;                   %Nozzle cant (degrees) - angle between nozzle and rocket that causes gimbaling

%set up menu
exit_flag=0;
hit_ground = 0;

    while (exit_flag==0)
        choice = menu('Altas V rocket launching simulation', '1. Display graphs of motion', '2. Display graphs of other changing variables',  '3. Export data to Excel', '4. Change solid rocket booster parameters', '5. Change atlas booster parameters', '6. Change simulation parameters','7. Change gimbal of rocket parameters','8. Exit simulation');
        switch (choice)

            case{1} 
                %clear current graphs
                clf;
                %Plots graphs of Altitude against time, Velocity against time, Acceleration against time and Vertical displacement against horizontal displacement                
                [t_data,accel_data,vel_data,alt_data,~,~,~,~,~,~,h_horizdata,h_vertdata,~,~,hit_ground] = calcs(A_booster_thrust_start,A_Isp,A_booster_propellant,A_booster_inert_mass,SR_booster_thrust_start,SR_booster_mass,SR_booster_inert_mass,SR_Isp,SR_boosters,j_time,payload,r_mass,dt,t_end,tg1,tg2,angle_gimbal,hit_ground);
                if hit_ground == 0
                    subplot(2,2,1);
                    plot(t_data,alt_data);
                    %Title of graph
                    title('Graph of altitude against time'); 
                    %x-axis label
                    xlabel('Time/ s'); 
                    %y-axis label
                    ylabel('Altitude/ m'); 
 
                    subplot(2,2,2);
                    plot (t_data,vel_data);
                    %Title of graph
                    title ('Graph of velocity against time'); 
                    %x-axis label
                    xlabel('Time/ s');
                    %y-axis label
                    ylabel('Velocity/ m/s');

                    subplot(2,2,3);
                    plot (t_data,accel_data);
                    %Title of graph
                    title ('Graph of acceleration against time'); 
                    %x-axis label
                    xlabel('Time/ s');
                    %y-axis label
                    ylabel('Acceleration/ m/s^2');
             
                    subplot (2,2,4);
                    plot (h_horizdata, h_vertdata);
                    %Title of graph
                    title ('Graph of vertical diplacement against horizontal displacement'); 
                    %x-axis label
                    xlabel('Horizontal displacement/ m');
                    %y-axis label
                    ylabel('Vertical displacement/ m');
                end
                
            case{2} 
                %clear current graphs
                clf;
                %Plots graphs of theta (angle between rocket and horizontal) against time, Drag against time, Thrust against time and Mass against time, Temperature against altitude and Mach number agaisnt drag coefficient               
                [t_data,~,~,alt_data,~,drag_data,thrust_data, m_data,theta_data,temp_data,~,~,Mach_data,Cd_data,hit_ground] = calcs(A_booster_thrust_start,A_Isp,A_booster_propellant,A_booster_inert_mass,SR_booster_thrust_start,SR_booster_mass,SR_booster_inert_mass,SR_Isp,SR_boosters,j_time,payload,r_mass,dt,t_end,tg1,tg2,angle_gimbal,hit_ground);
                 if hit_ground == 0              
                    subplot(2,3,1);
                    plot (t_data,theta_data);
                    %Title of graph
                    title ('Graph of theta against time'); 
                    %x-axis label
                    xlabel('Time/ s');
                    %y-axis label
                    ylabel('Theta/ degrees');
   
                    subplot(2,3,2);
                    plot (t_data,drag_data);
                    %Title of graph
                    title ('Graph of drag force against time'); 
                    %x-axis label
                    xlabel('Time/ s');
                    %y-axis label
                    ylabel('Drag force/ N');
              
                    subplot(2,3,3);
                    plot (t_data,thrust_data);
                    %Title of graph
                    title ('Graph of thrust against time'); 
                    %x-axis label
                    xlabel('Time/ s');
                    %y-axis label
                    ylabel('Thrust/ N');

                    subplot(2,3,4);
                    plot (t_data,m_data);
                    %Title of graph
                    title ('Graph of mass against time'); 
                    %x-axis label
                    xlabel('Time/ s');
                    %y-axis label
                    ylabel('Mass/ kg');
                    
                    subplot(2,3,5);
                    plot (temp_data,alt_data);
                    %Title of graph
                    title ('Temperature change with altitude')
                    %x-axis label
                    xlabel('Temperature / K');
                    %y-axis label
                    ylabel('Altitude/ m');
                    
                    subplot(2,3,6);
                    plot (Mach_data, Cd_data);
                    %Title of graph
                    title ('Drag coeffiecient against mach number');
                    %x-axis label
                    xlabel('Mach number');
                    %y-axis label
                    ylabel('Drag coefficient');
                 end
                    
            case{3}
                %Create an Excel spreadsheet with all the motion data  
                [t_data,accel_data,vel_data,alt_data] = calcs(A_booster_thrust_start,A_Isp,A_booster_propellant,A_booster_inert_mass,SR_booster_thrust_start,SR_booster_mass,SR_booster_inert_mass,SR_Isp,SR_boosters,j_time,payload,r_mass,dt,t_end,tg1,tg2,angle_gimbal,hit_ground);
                filename= 'Atlas_V_Rocket_Simulation_Data.xlsx';
                col_header={'Time/s','Acceleration/ms^-2','Velocity/ms^-1','Altitude/m'};
                xlswrite(filename, col_header, 'Sheet 1');
                xlswrite(filename, [t_data',accel_data',vel_data',alt_data'], 'Sheet 1','2');
                 
            case{4}  
                %clear current graphs
                clf;
                %while loop to ensure number of SRBs is between 0 and 5
                while 1
                    %Allow user to change solid rocket booster parameters for the simulation
                    [SR_boosters, SR_booster_thrust_start, SR_booster_mass, SR_booster_inert_mass,SR_Isp] = Change_SRB_parameters(SR_boosters, SR_booster_thrust_start, SR_booster_mass, SR_booster_inert_mass,SR_Isp);
                    if (SR_boosters>=0)&&(SR_boosters<=5)
                        break
                    end
                    msgbox('Please enter a value for no. of SRBs between 0-5');
                end
                
            case{5}   
                %clear current graphs
                clf;
                %Allow user to change Atlas booster rocket parameters for the simulation
                [A_booster_thrust_start, A_Isp,A_booster_propellant,A_booster_inert_mass] = Change_AB_parameters(A_booster_thrust_start, A_Isp, A_nozzel_diameter,A_booster_propellant, A_booster_inert_mass);
                
            case{6}   
                %clear current graphs
                clf;
                %Allow user to change system parameters for the simulation 
                [dt, t_end, payload] = Change_simuparameter(dt, t_end, payload);
                
            case{7}   
                %clear current graphs
                clf;
                %Allow user to change gimbaled thrust of rocket
                [tg1, tg2, angle_gimbal] = Change_gimbal_rocket(tg1, tg2, angle_gimbal);
                
            case{8}
                %Exit simulation
                exit_flag = 1; 
                
                %close figure window
                close(gcf)
              
            otherwise
                
disp('\nExiting simulation...\n\n');
        end
    end

    
function [t_data,accel_data,vel_data,alt_data,pres_data,drag_data,thrust_data, m_data,theta_data,temp_data,h_horizdata,h_vertdata,Mach_data,Cd_data,hit_ground] = calcs(A_booster_thrust_start,A_Isp,A_booster_propellant,A_booster_inert_mass,SR_booster_thrust_start,SR_booster_mass,SR_booster_inert_mass,~,SR_boosters,j_time,payload,r_mass,dt,t_end,~,~,~,~)
hit_ground = 0;

%Initial conditions
v = 0;          %rocket velocity (m/s)
h_vert = 0;     %rocket vertical displacement (m)
h_horiz = 0;    %horizontal displacment (m)
w = 0;          %rotational velocity of the rocket (rad/s)
theta = 90;     %angle between horizontal and rocket velocity (degrees)
phi = 0;        %angle between vertical and rocket velocity (degrees)
g = 9.81;       %accelertion due to gravity (m/s^2)
Cd = 0.3;       %drag coefficient
drag = 0;       %drag force acting on rocket (N)
A = 25.47;      %frontal area (m^2)

%Initial environmental conditions
T = 281.65;     %air temperature (K)
p = 101325;     %air pressure (Pa)

%Constants
L = 0.0065;     %Lapse rate (Km^-1)
Rho = 1;        %air density (kgm^-3)
re = 6371000;   %radius of earth (m)
G=6.67e-11;     %gravitational constant (Nm^2/kg^2)
Me = 5.9722e24; %mass of earth (kg)

% set up arrays
t_data = zeros(1, (t_end/dt));
alt_data = zeros(1, (t_end/dt));
vel_data = zeros(1, (t_end/dt));
accel_data = zeros(1, (t_end/dt));
pres_data = zeros(1,(t_end/dt));
drag_data = zeros(1,(t_end/dt));
thrust_data = zeros(1,(t_end/dt));
m_data = zeros(1, (t_end/dt));
theta_data = zeros(1,(t_end/dt));
temp_data = zeros(1,(t_end/dt));
h_horizdata = zeros(1,(t_end/dt));
h_vertdata = zeros(1,(t_end/dt));
Mach_data = zeros(1,(t_end/dt));
Cd_data = zeros(1, (t_end/dt));

%remove excess zeroes from arrays (when t_end>327.7s)
t_data(t_data == 0) = nan;
h_horizdata(h_horizdata == 0) = nan;
temp_data(temp_data == 0) = nan;
Mach_data(Mach_data == 0) = nan;

%Find initial total mass
m = (SR_boosters*SR_booster_mass)+A_booster_propellant+A_booster_inert_mass+payload+r_mass;

%main
for t=0:dt:t_end
    number = int16(t/dt)+1;
    
    %%FIND CURRENT VALUES FOR velocit ,acceleration, altitude %%
   
    %Check if propellant has run out
    if m > payload + r_mass + A_booster_inert_mass
        %Find thrust and new mass from thrust function
        [thrust,m_new] = engine_force(t,g,m,dt,SR_boosters,j_time,A_booster_thrust_start,SR_booster_thrust_start,SR_booster_inert_mass,SR_booster_mass,A_Isp,t_end);
    else
        thrust = 0;
    end
    
    %Find new theta, velocity and acceleration from runge function
    [theta,v,a,phi,w] = runge(v,dt,theta,g,A,m,thrust,Rho,phi,t,w,Cd);
    
    %Calculate new velocities
    v_vert = v*sind(theta);         %calc vertical component of new velocity
    v_horiz = v*cosd(theta);        %calc horizontal component of new velocity
    
    %Calculate new heights
    h_vert = h_vert + (v_vert*dt);     %calc vert component of height
    h_horiz = h_horiz + (v_horiz*dt);  %calc horizontal component of height
    
    %Calculate new altitude using pythagoras along the rockets trajectory
    alt = sqrt((h_vert+re)^2 + h_horiz^2) - re;
    
    %Add current values to arrays
    t_data(number) = t;                             %store new time data in array 
    alt_data(number) = alt;                         %store new altitude data in array
    vel_data(number) = v;                           %store new velocity data in array
    accel_data(number) = a;                         %store new acceleration data in array
    pres_data(number) = p;                          %store new pressure data in array
    drag_data(number) = drag;                       %store new drag data in array                       
    thrust_data(number) = thrust;                   %store new thrust data in array
    m_data(number) = m;                             %store new mass data in array
    theta_data(number) = theta;                     %store new theta data in array
    temp_data(number) = T;                          %store new temperature data in array
    h_horizdata(number) = h_horiz;                  %store new horrizontal displacment data in array
    h_vertdata(number) = h_vert;                    %store new vertical displacment data in array

    
    %%FIND NEW VALUES FOR NEXT ITERATION%%
    
    %Find new rho from drag function
    [drag,p,Rho,T,Mach,Cd,A] = find_drag(alt,g,v,L,p,T,dt,t,v_vert,SR_boosters,j_time); 
    
    %Find new acceleration due to gravity
    g=(G*Me)/((re+alt)^2);
       
    %Store new mach number and drag coefficient data in arrays
    Mach_data(number) = Mach;
    Cd_data(number) = Cd;
    
    %Set new mass to current mass
    m = m_new;
    
    %Test if rocket has crashed
    if alt <= 0 
        msgbox('Rocket hit the ground - Please change variables');
        hit_ground = 1;
        break
    end
     
end

end

function [drag,p,Rho,T,Mach,Cd,A] = find_drag(alt,g,v,L,p,T,dt,t,v_vert,SR_boosters,j_time)

    %Enviromental variables
    Po = 101325;        %launch pressure in Pa
    To = 288.15;        %launch temp in K
    R1 = 8.31447;       %Universal gas constant in J/(mol*K)
    R2 = 286;           %Universal gas constant in m^2/s^2/K
    M = 0.0289644;      %Molar mass of dry air in kg/mol 
    gama = 1.4;         %specific heat ratio

    %Area variables
    rPLF = 2.7;         %radius of payload fairing (m)
    ASRB = 0.855;       %area of SRBs (m^2)
    
    %Array of lapse rates
    Lapse_rate = [0.0065 0.0 -0.001 -0.0028 0.002];
    
    %Change of lapse rate, L
    if (alt <= 11000)                           %Troposphere
        L = Lapse_rate(1);
    elseif (alt > 11000) && (alt <= 20000)      %Tropopause
        L = Lapse_rate(2);
    elseif (alt > 20000) && (alt <= 32000)      %Stratosphere 1
        L = Lapse_rate(3);
    elseif (alt > 32000) && (alt <= 47000)      %Stratosphere 2
        L = Lapse_rate(4);
    elseif (alt > 47000) && (alt <= 51000)      %Stratopause
        L = Lapse_rate(2);
    elseif (alt > 51000) && (alt <= 71000)      %Mesosphere 1
        L = -Lapse_rate(4);
    elseif (alt > 71000) && (alt <= 84852)      %Mesosphere 2
        L = Lapse_rate(5);
    elseif (alt > 84852)                        %Mesopause & beyond
        L = Lapse_rate(2);
    end
 
    %Change of frontal area, A
    if (t >= j_time)                %Area before SRB jettison
        A = pi*(rPLF.^2);
    else                            %Area after SRB jettison
        A = pi*(rPLF.^2) + ASRB*SR_boosters;
    end
    
    %Change of temperature, T
    if alt == 0                     
        T = To;                     %Ground level temp        
    else
        if (L ~= 0.0)               %During atmospheric pause layers where temperature is roughly constant
            T = T - (L*v_vert*dt);
        end
    end

    %Change of pressure, p
    %Assume pressure change with the lapse rate of Mesosphere 2 to simulate exponential decay of pressure as predicted
    if (L == 0.0)       
        if (alt <= 20000)
            p = Po*((1-0.002*alt/To).^((g*M)/(R1*0.002)));
        elseif (alt > 20000) && (alt <= 51000)
            p = Po*((1-0.002*alt/To).^((g*M)/(R1*0.002)));
        end
    else
        p = Po*((1-L*alt/To).^((g*M)/(R1*L)));
    end
    
    %Change of density, Rho
    Rho = (p*M)/(R1*T);

    %Change of speed of sound, c
    c = (gama*R2*T).^(1/2);
    
    %Change of Mach Number with velocity of rocket
    Mach = v/c;
    
    %Change in drag coefficient, Cd
    %Calculated using values from online graphs and estimated equations
    %URL to website with data- https://space.stackexchange.com/questions/12649/how-can-i-estimate-the-coefficient-of-drag-on-a-saturn-v-rocket-a-simulator-or
    if (0<=Mach) && (Mach<0.5)
        Cd = 0.3-0.08*(Mach);
    elseif (0.5<=Mach) && (Mach<1.3)
        Cd = 0.26 + 0.3625*(Mach-0.5);
    elseif (1.3<=Mach) && (Mach<4.0)
        Cd = 0.55 - 0.1259*(Mach-1.3);
    elseif (4.0<=Mach) && (Mach<5.5)
        Cd = 0.21;
    elseif (5.5<=Mach) && (Mach<8.5)
        Cd = 0.21+0.01667*(Mach-5.5);
    else
        Cd = 0.26;
    end
    
    %Calculate new drag force
    drag = Cd*(Rho*(v.^2)/2)*A;
    
  
end

function [thrust,m_new] = engine_force(t,g,m,dt,SR_boosters,j_time,A_booster_thrust_start,SR_booster_thrust_start,SR_booster_inert_mass,SR_booster_mass,A_Isp,t_end)
    
    %thrust due to throttle profile based off Atlas V 521 User Guide graph (Page 2-19)  
    if (t > 30) && (t<=70)
        SR_booster_thrust = SR_booster_thrust_start * 0.93;        
        A_booster_thrust = A_booster_thrust_start * 0.93;
    elseif (t > 130) && (t<=200)
        SR_booster_thrust = SR_booster_thrust_start * 0.95;
        A_booster_thrust = A_booster_thrust_start * 0.95;
    elseif (t > 200) && (t<=220)
        SR_booster_thrust = SR_booster_thrust_start * 0.6;
        A_booster_thrust = A_booster_thrust_start * 0.6;
    elseif (t > 220) && (t<=230)
        SR_booster_thrust = SR_booster_thrust_start * 0.8;
        A_booster_thrust = A_booster_thrust_start * 0.8;
    elseif (t > 230) && (t<=t_end)
        SR_booster_thrust = SR_booster_thrust_start*(0.8 - ((t-230)*(0.4/(t_end-230))));
        A_booster_thrust = A_booster_thrust_start*(0.8 - ((t-230)*(0.4/(t_end-230))));
    else
        SR_booster_thrust = SR_booster_thrust_start;
        A_booster_thrust = A_booster_thrust_start;
    end
    
    %%Solid Rocket Booster%%
    %check if t = jettison time for SRBs
    if t == j_time
        m = m - (SR_booster_inert_mass*SR_boosters);  %remove SRB mass
    end
 
    %test if SRBs are still burning
    if t < j_time
        %reduce the mass by amount = to propellant burnt
        SR_mdot = ((SR_booster_mass - SR_booster_inert_mass)/j_time)*SR_boosters;
        m = m-(SR_mdot*dt);
        %set value of SR_thrust
        SR_thrust = SR_booster_thrust*SR_boosters;
    else 
        SR_thrust = 0;
    end
    
    %%Atlas Booster%%
    %reduce the mass by amount equal to propellant burnt 
    A_mdot = A_booster_thrust/(A_Isp*g);
    m_new = m - (A_mdot*dt);
    
    %Calculate overall thrust
    thrust = SR_thrust + A_booster_thrust;
    
end

function [theta,v,a,phi,w] = runge(v,dt,theta,g,A,m,thrust,Rho,phi,t,w,Cd)

rocket_height = 55.4;  %in m

%convert theta to radians
theta = theta * (pi/180);

%Set functions for theta and acceleration. Modelled as circular motion with
%acceleration towards centre of that circle. This approximation is used for
%the trajectory of rocket. URL - http://mei.org.uk/files/maw/herd.pdf
dtheta = @(v,theta) (-1*g*cos(theta))/v;
dacc = @(v,theta) (thrust - (0.5*Cd*Rho*A*v.^2) - (m*g*sin(theta)))/m;

%Set case for when v = 0 
if v == 0               
    dtheta = @(v,theta) 0;
end

%Set boosters to cause the rocket to gimbal between certain times due to
%angle of thrusters.
%Alpha is angular acceleration of rocket around its centre of mass.
%Rocket modelled as a rod using first principles to calculate moment of
%inertia.
if (tg1<t)&&(tg2>t)
    dtheta = @(v,theta) 0;
    alpha = (6*thrust*sind(angle_gimbal))/(m*rocket_height);
    phi = phi + (w*(dt)) + 0.5*alpha*((dt)^2);
    w = (alpha*dt) + w;
    theta = pi/2 - phi;
end

%Apply runge kutta method to find new acceleration and theta values
kg1 = dt * dtheta(v,theta);
kf1 = dt * dacc(v,theta);

kg2 = dt * dtheta(v+(0.5*kf1),theta+(0.5*kg1));
kf2 = dt * dacc(v+(kf1*0.5),theta+(0.5*kg1));

kg3 = dt * dtheta(v+(0.5*kf2),theta+(0.5*kg2));
kf3 = dt * dacc(v+(kf2*0.5),theta+(0.5*kg2));

kg4 = dt * dtheta(v+kf3,theta+kg3);
kf4 = dt * dacc(v+kf3,theta+kg3);

theta = theta + ((1/6)*(kg1+(2*kg2)+(2*kg3)+kg4));
v = v + ((1/6)*(kf1+(2*kf2)+(2*kf3)+kf4));

%find acceleration from change in velocity
a = ((1/6)*(kf1+(2*kf2)+(2*kf3)+kf4))/dt;

%convert theta to degrees
theta = theta * (180/pi);

end

%Changes Solid Rocket Booster parameters
function [SR_boosters, SR_booster_thrust, SR_booster_mass, SR_booster_inert_mass,SR_Isp] = Change_SRB_parameters (SR_boosters, SR_booster_thrust, SR_booster_mass, SR_booster_inert_mass,SR_Isp)
    
    %creates dialogue box with current values
    prompt = {'Number of solid rocket boosters','Thrust of solid rocket boosters (N)','Total mass of solid rocket boosters (kg)','Inert mass of solid rocket boosters(kg)','Specific impulse (ISP) of solid rocket boosters (s)'};
    dlg_title = 'Current solid rocket booster parameters:';
    num_lines = 1;
    default_values = {sprintf('%.f',SR_boosters),sprintf('%.1f',SR_booster_thrust),sprintf('%.1f',SR_booster_mass),sprintf('%.1f',SR_booster_inert_mass),sprintf('%.1f',SR_Isp)};
    answer = inputdlg(prompt,dlg_title,num_lines,default_values);        
    
    %sets solid rocket booster parameters to values in dialogue box
    SR_boosters = str2double(answer(1));
    SR_booster_thrust = str2double(answer(2));
    SR_booster_mass = str2double(answer(3));
    SR_booster_inert_mass = str2double(answer(4));
    SR_Isp = str2double(answer(5));
    
end

%Changes Atlas booster rocket parameters
function [A_booster_thrust, A_Isp,A_booster_propellant,A_booster_inert_mass] = Change_AB_parameters(A_booster_thrust, A_Isp,A_booster_propellant, A_booster_inert_mass)

    %creates dialogue box with current values
    prompt = {'Atlas booster thrust at sea level (N)','Specific impulse (ISP) of Atlas booster (s)', 'Mass of the Atlas booster propellant (kg)', 'Inert mass of Atlas booster (kg)'};
    dlg_title = 'Current solid rocket booster parameters:';
    num_lines = 1;
    default_values = {sprintf('%.1f',A_booster_thrust),sprintf('%.1f',A_Isp),sprintf('%.1f',A_booster_propellant),sprintf('%.1f',A_booster_inert_mass)};
    answer = inputdlg(prompt,dlg_title,num_lines,default_values);
    
    %sets solid rocket booster parameters to values in dialogue box
    A_booster_thrust = str2double(answer(1));
    A_Isp = str2double(answer(2));
    A_booster_propellant = str2double(answer(3));
    A_booster_inert_mass = str2double(answer(4));
    
end

function [dt, t_end, payload] = Change_simuparameter(dt, t_end, payload)

    %creates dialogue box with current values
    prompt = {'Integration step (s)','Time of jettison of Atlas booster (s)','Mass of cargo (kg)'};
    dlg_title = 'Current system parameters:';
    num_lines = 1;
    default_values = {sprintf('%.2f',dt),sprintf('%.1f',t_end),sprintf('%.1f',payload)};
    answer = inputdlg(prompt,dlg_title,num_lines,default_values);
    
    %sets system parameters to values in dialogue box
    dt = str2double(answer(1));
    t_end = str2double(answer(2));
    payload = str2double(answer(3));
   
end


function [tg1, tg2, angle_gimbal] = Change_gimbal_rocket(tg1, tg2, angle_gimbal)

    %creates dialogue box with current values
    prompt = {'Time when thrusters begin to gimbal (s)','Time when thrusters stop gimbaling (s)','Angle of gimbal/Nozzle cant (degrees)'};
    dlg_title = 'Current gimbaling parameters';
    num_lines = 1;
    default_values = {sprintf('%.2f',tg1),sprintf('%.2f',tg2),sprintf('%.2f',angle_gimbal)};
    answer = inputdlg(prompt,dlg_title,num_lines,default_values);
    
    %sets system parameters to values in dialogue box
    tg1 = str2double(answer(1));
    tg2 = str2double(answer(2));
    angle_gimbal = str2double(answer(3));
end
end

