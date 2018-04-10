% Solve the system of equations to find performance metrics v. time as
%   described by Fanchiang --> Run Monte Carlo Simulation to Optimize Design
clear all
close all
%% Import the matrices defined in excel

num_of_crew = 2; % define number of crewmembers
num_performance_elements = 26; % define number of crewmember performance resource elements
num_crew_activity_type = 21; % define number of activity types for the crewmembers 
num_design_choices = 3; %define number of design spacecraft design choices to compare
num_design_parameters =43;
num_mission_days = 180; %number of mission days
design_input_length = 172;
vehicle_environment = 2; % choose vehicle environment (1 = 1 gravity, 2 = 0 gravity, 3 = hypergravity)

% Set the time units by changing these two variables (make sure they match)
time_unit =24; % Set this time unit as the scaling for the time axis (24*30=by months, 24=by days, 1 = by hours)
x_label ='t, day'; % t, hr, day, month
total_mission_time = 24*num_mission_days/time_unit; % Specify total mission length by 24hrs x # of mission days

% load the design to crew mapping matrix
Design2CrewImpact_Mapping = xlsread('Design2CrewPerformance_Mapping.xlsx',1,'A1:Z228');
% Design2CrewImpact_Mapping =10*Design2CrewImpact_Mapping;
Design2CrewImpact_Mapping(design_input_length,:) = NaN*ones(1,num_performance_elements);

% load the design to task mapping matrix
Design2Task_Mapping = xlsread('Design2Task_Mapping.xlsx');
Design2Task_Mapping(design_input_length,:) = NaN*ones(1,num_crew_activity_type);

% load the task to crew mapping matrix
Task2Crew_Mapping = xlsread('Task2CrewPerformance.xlsx');% Uncomment this!
% Task2Crew_Mapping = 10*Task2Crew_Mapping;
% Task2Crew_Mapping(228,:) = NaN*ones(1,26); % Uncomment this!
% Task2Crew_Mapping = zeros(21,26); % Specify Task Mapping % Comment this out!

% Generate task time histories
Crew_Task_Function1day = xlsread('TaskFunction.xlsx','NoTasks'); % Imports TaskFunction list for 1 CM
Crew_Task_Function = repmat(Crew_Task_Function1day,num_mission_days,1);


% Generate Crew Time Degredation due to vehicle environment

if vehicle_environment == 1;
    Crew_Time_Degredation_Array = {@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0};
    Crew_Time_Degredation = @(t)cellfun(@(f)f(t),Crew_Time_Degredation_Array);
    
elseif vehicle_environment == 2;
    Crew_Time_Degredation_Array = {@(t)-8.848303-7.824956*erf((t-3.695185)/1.473982)-0.889235*erf((t-9.548112)/2.805284),@(t)-(3.666537+5.665088*erf((t-0.164422)/0.268212)+2.469145*erf((t-0.348657)/3.991582)),@(t)-2*t,@(t)-2*t,@(t)-2*t,@(t)-2*t,@(t)-2*t,@(t)-9.323626-8.414837*erf((t-2.898860)/1.584169)-1.365004*erf((t-6.670830)/4.679852),@(t)-2*t,@(t)-2*t,@(t)-2*t,@(t)-2*t,@(t)-2*t,@(t)-4*t,@(t)-4*t,@(t)-4*t,@(t)-4*t,@(t)-4*t,@(t)-4*t,@(t)-4*t,@(t)-4*t,@(t)-4*t,@(t)-4*t,@(t)10*cos(2*pi*t*3/(4*(num_mission_days/30)))+5-4*t,@(t)10*cos(2*pi*t*3/(4*(num_mission_days/30)))+5-4*t,@(t)10*cos(2*pi*t*3/(4*(num_mission_days/30)))+5-4*t};
    Crew_Time_Degredation = @(t)cellfun(@(f)f(t),Crew_Time_Degredation_Array);
    
elseif vehicle_environemtn == 3;
    Crew_Time_Degredation_Array = {@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0,@(t)0};
    Crew_Time_Degredation = @(t)cellfun(@(f)f(t),Crew_Time_Degredation_Array);
end
    
% Crew_Time_Degradation_Array = xlsread('Time_Degradation.xlsx');% Import time degradation, assume linear:
                      % Metrics = Degredation*t
                      

% Remove Nan Values                      
a = [isnan(Design2CrewImpact_Mapping(:,1)),[1:design_input_length]'];
a(a(:,1) == 0,:) = [];
Design_Mapping_Nan_Removed = Design2CrewImpact_Mapping;
Design_Mapping_Nan_Removed(a(:,2),:) = [];

Design2Task_Mapping(a(:,2),:) = [];

%% Specify time vectors and input 

N = num_crew_activity_type; % Number of tasks
time_vec = linspace(0,24*num_mission_days,24*num_mission_days); % Specify time range (assume 1 day of 24 hrs to start)
time_step = time_vec(2) - time_vec(1); % Initialize time step



%% Random Design Generator
% Randomly pick input vectors and solve system to optimize baseline
for n = 1:num_design_choices;%1e5 % Setup Monte Carlo for the number of design choices
  
    for j=1:num_of_crew;
%     Time_Degradation = Crew_Time_Degradation(:,j);

    Task_Function = (Crew_Task_Function(:,1+num_crew_activity_type*(j-1):num_crew_activity_type*j))';

        for i = 1:num_design_parameters % Number of input design parameters
            Input = zeros(1,4); %Initialize input matrix
            %Generate two vectors
            if n == 1
                Index = 1;%randi(4);
            elseif n == 2
                Index = 2;
            else
                Index = randi(2); % Generate random
            end

                % If NaN value is present indicating no choice, re-select
                if isnan(Design2CrewImpact_Mapping(4*i-3+(Index-1),1))
                    while isnan(Design2CrewImpact_Mapping(4*i-3+(Index-1),1))
                        Index = randi(4);
                    end
                end

            Input(Index) = 1; % Set value from index to 1, leave rest zero
            Sample_Vector(4*i-3:4*i) = Input; % Note that this sample vector could be specified (this is the design)

        end

        % Remove the NaN rows %%%%
            Sample_Vector(a(:,2)) = [];
        %%%%%%%%


    %%%%%%%%%%%%%%%%%%% Run this Portion if you specify Design Vector and Comment out n for loop and random design generator (ctrl-r/ctrl-t) %%%%%%%%    
    % Solve System
    % Computes Baseline Decay due to Vehicle Design Selection
    Design_Vector(:,n) = Sample_Vector';
    Baseline_Metrics = Sample_Vector*Design_Mapping_Nan_Removed;
    Initial_Baseline_Metrics(:,n,j) = 100 + Baseline_Metrics';

    % Determine the design to task weighting coefficients
    Design_Task_Weights = (Sample_Vector*Design2Task_Mapping)'*ones(1,num_performance_elements)/num_design_parameters; % Divide by 57 to make degradation amount set to 1 if design induces no change
    Task2Crew_Mapping_Weighted = Task2Crew_Mapping.*Design_Task_Weights; % Re-Scale Task Mapping Matrix by Design to Task weighting coefficients

    Task_Decay = zeros(1,num_performance_elements); % Initialize Task Function
    Baseline_Decay = zeros(1,num_performance_elements); % Initialize Baseline Decay

    %Propogate forward in time
        for t = 1:length(time_vec);
             Baseline_Decay = Baseline_Decay + Baseline_Metrics*time_step;
             Task_Decay = Task_Decay + (Task_Function(:,t)')*Task2Crew_Mapping_Weighted*time_step; %Find Task Function effect running sum over time expired
             temp = Baseline_Decay' + Task_Decay' + Crew_Time_Degredation(time_vec(t)/(24*30))';
             %temp(temp+100>100)=0;
             
             Time_Metrics(:,t,n,j) =  temp; % Find time dependent metrics (divide time vec by hours/month)
             
        end

    Performance(:,n,j) = 100 + squeeze(Time_Metrics(:,end,n,j)); %Capture end of mission performance
    
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% Plot Results

% Examine final mission performance factor
figure('position',[50 50 1000 300])

color_vec = {'ok','xk','ob','xb','or','xr'};
for j = 1:num_of_crew
    for i = 1:n
        plot(1:num_performance_elements,squeeze(Performance(:,i,j)),char(color_vec(j*2-1))); %Plot end of mission performance
        set(gca,'XTick',[1:num_performance_elements]);
        hold on

        plot(1:num_performance_elements,squeeze(Initial_Baseline_Metrics(:,i,j)),char(color_vec(j*2)));
    end
end
set(gca,'XTickLabel',{'P_B','P_C','P_D','P_F','P_H','P_R','P_I','P_M','P_N','P_P','P_V','P_E','P_X','C_R','C_B','C_E','C_D','C_T','C_S','C_L','C_P','C_V','C_M','Y_B','Y_A','Y_M'}');
axis([0 num_performance_elements+1 -10 150])
ylabel('Baseline, %')
grid on

plot([0,num_performance_elements+1],[0,0],'r','linewidth',2)
hold off

p=legend('Crew1_F','Crew1_0', 'Crew2_F', 'Crew2_0', 'Crew3_F', 'Crew3_0');
set(p,'location','eastoutside')


%% Examine time history of single mission performance factor
% Plots physiology performance factors for all crewmembers
for i = 1:n
    figure('name','Physiolgical Performance','position',[0 0 500 1000])
    for k = 1:num_of_crew
        for j = 1:13
                subplot(3,1,k);plot(time_vec/time_unit,100 + squeeze(Time_Metrics(j,:,i,k)));
                hold on
        end

        axis([0 total_mission_time -10 150])
        xlabel(x_label)
        ylabel('Physiological Performance, %')
        grid on
        subplot(3,1,k);plot(time_vec/time_unit,100+mean(squeeze(Time_Metrics(1:13,:,i,k))),'--k','linewidth',2);
        subplot(3,1,k);plot([0,total_mission_time],[0,0],'r','linewidth',2)

        h=legend('P_B','P_C','P_D','P_F','P_H','P_R','P_I','P_M','P_N','P_P','P_V','P_E','P_X','mean');
        set(h,'location','eastoutside')
        
        Phy_design_crew_mean(k,:,i) = 100+mean(squeeze(Time_Metrics(1:13,:,i,k)));
    end
end

% Plots cognitive performance metrics for all crewmembers
for i = 1:n
    figure('name','Cognitive Performance','position',[500 0 500 1000])
    for k = 1:num_of_crew
        for j = 14:23
                subplot(3,1,k);plot(time_vec/time_unit,100 + squeeze(Time_Metrics(j,:,i,k)));
                hold on
        end
            axis([0 total_mission_time -10 150])
            xlabel(x_label)
            ylabel('Cognitive Performance, %')
            grid on
            subplot(3,1,k);plot(time_vec/time_unit,100+mean(squeeze(Time_Metrics(14:23,:,i,k))),'--k','linewidth',2);
            subplot(3,1,k);plot([0,total_mission_time],[0,0],'r','linewidth',2)
            
            h=legend('C_R','C_B','C_E','C_D','C_T','C_S','C_L','C_P','C_V','C_M','mean');
            set(h,'location','eastoutside')
            
            Cog_design_crew_mean(k,:,i) = 100+mean(squeeze(Time_Metrics(14:23,:,i,k)));
    end
end

% Plots psychological metrics for all crewmembers
for i = 1:n
    figure('name','Psychological Performance','position',[1000 0 500 1000])
    for k = 1:num_of_crew
        for j = 24:26
                subplot(3,1,k);plot(time_vec/time_unit,100 + squeeze(Time_Metrics(j,:,i,k)));
                hold on           
        end
            axis([0 total_mission_time -10 150])
            xlabel(x_label)
            ylabel('Psychological Performance, %')
            grid on
            subplot(3,1,k);plot(time_vec/time_unit,100+mean(squeeze(Time_Metrics(24:26,:,i,k))),'--k','linewidth',2);
            subplot(3,1,k);plot([0,total_mission_time],[0,0],'r','linewidth',2)
            
            h=legend('Y_B','Y_A','Y_M','mean');
            set(h,'location','eastoutside')
           
            Psy_design_crew_mean(k,:,i) = 100+mean(squeeze(Time_Metrics(24:26,:,i,k)));
    end   
end
hold on;

%% Plot Means
% if num_of_crew ==2 && num_design_choices ==3;
%         color_vec = {'sk', '-sr','ok', 'or', 'xk', 'xr'};
% 
% elseif num_of_crew ==1 && num_design_choices ==3;
%         color_vec = {'ok', 'or','ob'};
%         
% else
%         color_vec = {'-sk','-sr','-sb','-ok','-or','-ob','-k','-xr','-xb'};
% end
color_vec = {'-k','-r','-b','--k','--r','--b','-.k','-.r','-.b'};
        
figure('name','Physiological Performance Means','position',[0 0 1000 350])
for i = 1:n
    for k = 1:num_of_crew
        clean_phy_mean = smooth(Phy_design_crew_mean(k,:,i),72);
        plot(time_vec/time_unit,clean_phy_mean,char(color_vec(3*i-3+k)),'linewidth',1);
        hold on   
    end     
end
axis([0 total_mission_time -10 150]);
grid on
plot([0,total_mission_time],[0,0],'r','linewidth',2)
xlabel(x_label)
ylabel('Mean Physiological Performace')
% h = legend('Crew1_{n=1}','Crew2_{n=1}','Crew1_{n=2}','Crew2_{n=2}','Crew1_{n=3}','Crew2_{n=3}');%,'Crew1_{n=3}','Crew2_{n=3}','Crew3_{n=3}');
% h = legend('0g,optimal design, no tasks','0g,optimal design, 6hrs sleep, no tasks','0g,optimal design, high workload');
% h = legend('0g,optimal design, no tasks','0g,non-ideal design, no tasks', '0g,mixed design, no tasks');
h = legend('0g,optimal design, some tasks','0g,optimal design,  high workload', '0g,non-ideal design, some tasks','0g,non-ideal design, high workload','0g,mixed design, some tasks', '0g, mixed design, high worload');
set(h,'location','eastoutside')

figure('name','Cognitive Performance Means','position',[500 0 1000 350])
for i = 1:n
    for k = 1:num_of_crew
        clean_cog_mean = smooth(Cog_design_crew_mean(k,:,i),72);
        plot(time_vec/time_unit,clean_cog_mean,char(color_vec(3*i-3+k)),'linewidth',1);
        hold on
    end        
end
axis([0 total_mission_time -10 150]);
grid on
plot([0,total_mission_time],[0,0],'r','linewidth',2)
xlabel(x_label)
ylabel('Mean Cognitive Performance, %')
% h = legend('Crew1_{n=1}','Crew2_{n=1}','Crew1_{n=2}','Crew2_{n=2}','Crew1_{n=3}','Crew2_{n=3}');%,'Crew1_{n=3}','Crew2_{n=3}','Crew3_{n=3}');
%h = legend('0g,optimal design, no tasks','0g,optimal design, 6hrs sleep, no tasks','0g,optimal design, high workload');
% h = legend('0g,optimal design, no tasks','0g,non-ideal design, no tasks', '0g,mixed design, no tasks');
h = legend('0g,optimal design, some tasks','0g,optimal design,  high workload', '0g,non-ideal design, some tasks','0g,non-ideal design, high workload','0g,mixed design, some tasks', '0g, mixed design, high workload');
set(h,'location','eastoutside')

figure('name','Psychological Performance Means','position',[1000 0 1000 350])
for i = 1:n
    for k = 1:num_of_crew
        clean_psy_mean = smooth(Psy_design_crew_mean(k,:,i),72);
        plot(time_vec/time_unit,clean_psy_mean,char(color_vec(3*i-3+k)),'linewidth',1);
        hold on
    end
end
axis([0 total_mission_time -10 150]);
grid on
plot([0,total_mission_time],[0,0],'r','linewidth',2)
xlabel(x_label)
ylabel('Mean Psychological Performance, %')
% h = legend('Crew1_{n=1}','Crew2_{n=1}','Crew3_{n=1}','Crew1_{n=2}','Crew2_{n=2}','Crew3_{n=2}');%,'Crew1_{n=3}','Crew2_{n=3}','Crew3_{n=3}');
% h = legend('Crew1_{n=1}','Crew2_{n=1}','Crew1_{n=2}','Crew2_{n=2}','Crew1_{n=3}','Crew2_{n=3}');%,'Crew1_{n=3}','Crew2_{n=3}','Crew3_{n=3}');
% h = legend('0g,optimal design, no tasks','0g,non-ideal design,  tasks', '0g,mixed design, no tasks');
h = legend('0g,optimal design, some tasks','0g,optimal design,  high workload', '0g,non-ideal design, some tasks','0g,non-ideal design, high workload','0g,mixed design, some tasks', '0g, mixed design, high worload');
set(h,'location','eastoutside')

Task_sum=Task_Function;
Task_sum(1,:)=[];

Task_sum(Task_sum~=1)=0;

num_of_tasks =sum(Task_sum(:));

   
fid =fopen('CrewAccomodationStatus.txt','w');
fprintf(fid,'\t\t CREW ACCOMMODATION STATUS \n');
fprintf(fid,'%.1f %% Physiological Status \n', Phy_design_crew_mean(end));
fprintf(fid,'%.1f %% Cognitive Status \n', Cog_design_crew_mean(end));
fprintf(fid,'%.1f %% Psychological Status \n', Psy_design_crew_mean(end));
fprintf(fid, '\t\t CREW UTILIZATION \n');
fprintf(fid,'%.2d Physiological Utilization \n', (100-Phy_design_crew_mean(end))/num_of_tasks);
fprintf(fid,'%.2d Cognitive Utilization \n', (100-Cog_design_crew_mean(end))/num_of_tasks);
fprintf(fid,'%.2d Psychological Utilization \n', (100-Psy_design_crew_mean(end))/num_of_tasks);