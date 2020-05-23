
#Agathe Nidriche 04.20


#TO USE THE PROGRAMM


#Parameters one needs to fill : The first two parameters, which means :
#Change HP to 1,2,3,4 according to the hypothesis you want to simulate
#Change Graph to Display_Lines or Display_Homogeneous according to the way you want to display graphs
#If you want to save figures, change the current path
#Change K, nb of cells, if needed

#ADVICES :
#HP1 : k=1 otherwise some weird checkboard pattern
#HP2 : 200 cells instead of 70 to get values close to the article
#HP3 : 
#HP4 : 


#Optionnal changes :
#Physics parameters can be modified but here they correspond to the values adviced by the authors
#The number of cells is chosen such that confluence happens with cells reaching a 10 um diameter and uniform repartition, others are "washed out"
#Thus, same for the array containing cells
#Near the graph section the x axis rescaling parameters can be modified. Care about min max values


#No need to change any parameter in the additionnal functions called or in the script

#Functions to add to the folder :

#derivative_HP1, derivative_HP2, derivative_HP3, derivative_HP4
#Calc_v_s_e_HP1, Calc_v_s_e_HP2, Calc_v_s_e_HP3, Calc_v_s_e_HP4
#Display_Homogeneous, Display_Lines
#Bluewhitered (This one I downloaded from a Matlab user hub)


#-------------------------------------Global variables----------------------------------------

#Main axis
global Xi;#Positions of nodes over time
global T; #Time arrayindx_t_threshold_i=find(T(t-1)<=t_threshold_matrix(:,i) && t_threshold_matrix(:,i)<T(t),1);     

#useful matrices/arrays for V,S,E calculation
global flag_threshold_reached; # HP3,4:array of flags indicating whether epsilon_s was reached for all cells (changed over time in ODE45)
global T_eps_s;# HP3,4:array of times indicating when epsilon_s was reached for all cells (changed over time in ODE45)
global t_threshold; # HP4 : array of last times where epsilon_th was reached for each cell (changed over time in ODE45)
global K; # HP4 :array of elastic constants for each cells (changed over time in ODE45)
global t_threshold_matrix; #Binary matrix indicating couples of cell/time where epsilon_th is reached

#counts t in ODE45
global count_t; #Counts units of time


tic #Starts clock


#-------------------------------------Initialization----------------------------------------


#--------derivative fct info/graph choice : has to be filled §

HP=4; #Which hypothesis is used :1,2,3 or 4
Type='ODE45'; #solver type

#Chose one option out of the 2

#Graph='Display_Lines'; #Displays lines corresponding to nodes positions, more accurate
Graph='Display_Homogeneous'; #Colours are more homogenous for a continous 2D graph


av_nb_of_cells=300; #Number of cells, 10 um diamater per cell (quite a little size, but fits to the surface chosen
#By the authors of the paper).



#-----Physical parameters, extracted from  Supplementary doc page 12

f=4;#Self propelling force in nN
k=4;#Elastic spring, nN
n=9.2 ;#Viscosity nN.min/micrometers
epsilon_s=0.5; #for HP3 and 4 : strain threshold to reach to initiate self propelling forces
epsilon_th=1; #for HP4 : strain threshold to reach to initiate reinforcement/fluidization cycle
delta_t=250; #for HP4 : typical time of a reinforcement/fluidization cycle, min

#-----Parameters for ODE solver and plots

t_max=800 ;#Maximum time considered, min
typical_surf=100 ;#Typical surface of a cell to calculate stress, micrometers square
count_t=1;#initialisation of time counting for ODE45
T_array=0:0.5:t_max; #Time array used in ODE45

#---Parameters for cell seeding and substrate surface

length_PDMS=700 ;#Micrometers, estimated from kymographs in the paper (initial length of PDMS before removal)
av_nb_of_nodes=av_nb_of_cells+1;
Xi_array=linspace(-length_PDMS/2,length_PDMS/2,av_nb_of_nodes);  #This is the array of the positions of the cells at time t=0 min
initial_delta_x=Xi_array(2)-Xi_array(1); #Distance between cells at rest, micrometers

#----HP3
flag_threshold_reached=zeros(av_nb_of_nodes,1);# HP3,4:array of flags indicating whether epsilon_s was reached for all cells (changed over time in ODE45)
T_eps_s=zeros(av_nb_of_nodes,1);# HP3,4:array of times indicating when epsilon_s was reached for all cells (changed over time in ODE45)

#----HP4 

#----Polynomial function of reinforcemet/fluidization for HP4 (centered on delta_t/2)

Polynome = polyfit([0,delta_t/2,delta_t],[k,3*k,k],2);

#----Initialisation of matrices, for HP4

K=ones(av_nb_of_cells,1)*k;# HP4 :array of elastic constants for each cells (changed over time in ODE45)
t_threshold=ones(av_nb_of_cells,1)*t_max;# HP4 : array of last times where epsilon_th was reached for each cell (changed over time in ODE45)
t_threshold_matrix=zeros(10000,av_nb_of_cells);


#-------Folder where figures are stocked

#Change the name of the current path
cd 'C:\Users\agatn\Desktop\Version definitive codes Octave' 

Creation_date=date; 
#Change the name of folder_figs
folder_figs=strcat(Creation_date,'_HP=',num2str(HP),'_nb_of_cells=',num2str(av_nb_of_cells),'_delta_t=',num2str(delta_t),'_f=',num2str(f),'_k=',num2str(k))
mkdir(folder_figs);
cd(folder_figs);


#-----------------------------------Solving of the equations----------------------------------


#----Call to ODE45
#Hypothesis proposed by the paper

if HP==1
      #T_array=0:step_t:t_max;
      [T,Xi]=ode45(@(T,Xi)final_derivative_HP1(T,Xi,n,f,initial_delta_x,av_nb_of_nodes,k),T_array,Xi_array); #ODE solver
      [V,S,E_point]=Calc_v_s_e_HP1(T,Xi,n,f,initial_delta_x,av_nb_of_nodes,av_nb_of_cells,k,typical_surf); #Calculation V,S,E

elseif HP==2
      [T,Xi]=ode45(@(T,Xi)final_derivative_HP2(T,Xi,n,f,initial_delta_x,av_nb_of_nodes,k,typical_surf),T_array,Xi_array);#ODE solver
      [V,S,E_point]=Calc_v_s_e_HP2(T,Xi,n,f,initial_delta_x,av_nb_of_nodes,av_nb_of_cells,k,typical_surf);#Calculation V,S,E

elseif HP==3
      [T,Xi]=ode45(@(T,Xi)final_derivative_HP3(T,Xi,n,f,initial_delta_x,av_nb_of_nodes,k,typical_surf,epsilon_s,flag_threshold_reached),T_array,Xi_array,T_eps_s);#ODE solver
      [V,S,E_point]=Calc_v_s_e_HP3(T,Xi,n,f,initial_delta_x,av_nb_of_nodes,av_nb_of_cells,k,typical_surf,epsilon_s,flag_threshold_reached,T_eps_s);#Calculation V,S,E
   

elseif HP==4
      [T,Xi]=ode45(@(T,Xi)final_derivative_HP4(T,Xi,n,f,initial_delta_x,av_nb_of_nodes,av_nb_of_cells,k,epsilon_s,K,delta_t,epsilon_th,t_threshold,Polynome,t_threshold_matrix,flag_threshold_reached,T_eps_s,count_t),T_array,Xi_array);#ODE solver
      t_threshold_matrix=t_threshold_matrix(1:count_t,:);
      [V,S,E_point]=Calc_v_s_e_HP4(T,Xi,n,f,initial_delta_x,av_nb_of_nodes,av_nb_of_cells,k,typical_surf,epsilon_s,T_eps_s,t_threshold_matrix,Polynome,delta_t,epsilon_th,t_max)#Calculation V,S,E

else
       warning('Please Enter 1,2,3 or 4 for HP value')

endif


#----Keeping only the filled part of the matrices


#----------------------Plot of a few distances of cells over time------------------------#

%{
#Displays positions of cells according to the time axis
for j=1:av_nb_of_cells/10
  plot(T,Xi(:,10*j-9))
  hold on
endfor
%}

#--------------------Scaling of V,E,S to plot according to an unique X array-----------------------#


#----------Values, to change for the graph conditions


min_X_sorted_val=-900; #min x axis. Warning : If execution issue is a problem in display function , might need to be increased
max_X_sorted_val=-min_X_sorted_val; #max x axis 
delta_interval=2; #Step size


round_precision=0.01;
length_X_sorted=(abs(min_X_sorted_val)+abs(max_X_sorted_val))*(1/delta_interval); # length of X array


#----------Arrays we use to plot V,E,S

X_sorted=linspace(min_X_sorted_val,max_X_sorted_val,length_X_sorted); #X array we use to rescale
V_sorted=zeros(size(T),size(X_sorted,2));#V scaled to X_sorted
E_point_sorted=zeros(size(T),size(X_sorted,2)-1);#E_point scaled to X_sorted
S_sorted=zeros(size(T),size(X_sorted,2));#S scaled to X_sorted


#----------Rounding of axis to improve intervalls sorting

Xi(:,:)=round(Xi(:,:)*(1/round_precision))*round_precision; #All Xi values are rounded to 2nd digit
X_sorted(:,:)=round(X_sorted(:,:)*(1/round_precision))*round_precision;

#----------Linear array, size=nb of cells
idx_cells=1:av_nb_of_cells;
idx_nodes=1:av_nb_of_nodes;


#------------------------------------Graph Options-------------------------------------------

size_T=size(T,1);

  
if strcmp(Graph,'Display_Lines')
  #---------- 1) Displays lines  

  #(notwithstanding x large intervals, there are still too few cells on a line to see a uniform graph.)
  #But it is more accurate than homogeneous graph
  
      [V_sorted,S_sorted,E_point_sorted]=Display_Lines(Xi,X_sorted,idx_cells,idx_nodes,size_T,V,E_point,S,V_sorted,E_point_sorted,S_sorted);

  
elseif strcmp(Graph,'Display_Homogeneous')
#---------- OR 2) Displays homogeneous 2D graph

  #Colors the plot with the value of V,E,S of the first cell between two cells to get 
  #an homogeneous 2D graph instead of lines (less accurate, better for qualitative aspect)%{
      [V_sorted,S_sorted,E_point_sorted]=Display_Homogenous(Xi,X_sorted,idx_cells,idx_nodes,size_T,V,E_point,S,V_sorted,E_point_sorted,S_sorted);
  
else
       warning('Please Enter A Graph Mode')
       
endif



#------------------------------------2D Graphs----------------------------------

#---------------------Velocity

figure(2)
imagesc(X_sorted,T,V_sorted); #2D graph
colormap(Bluewhitered(300)); #Fonction for colormap with white centered on zero (efficient in our case)
c=colorbar;
ylabel(c, 'v, um/min','fontsize',14) #Labelling of the colorbar
set(gca,'xlim',[-600 600],'ylim',[0 800],'YDir','reverse'); #Reverses y axis
title(strcat("Velocity of the cells, um/min, HP=",num2str(HP)," , ode solver=",Type),'fontsize',18,'fontweight','b','color','b');
xlabel('Distance from the midline, um','fontweight','b','color','r'),;
ylabel('Time of observation, min','fontweight','b','color','r');
velocity=gcf;

#downloads the figure in the folder
print(velocity,strcat('velocity_nb_of_cells=',num2str(av_nb_of_cells),'_delta_t=',num2str(delta_t),".pdf"))


#--------------------Strain rate

figure(3)

imagesc(X_sorted(2:end),T,E_point_sorted); #kymograph
colormap(Bluewhitered(250));
c=colorbar;
ylabel(c,'strain rate, s^-1','fontsize',14) #Labelling of the colorbar
set(gca,'xlim',[-600 600],'ylim',[0 800],'YDir','reverse');
title(strcat("Strain rate, s^-1, HP=",num2str(HP)," , ode solver=",Type),'fontsize',18,'fontweight','b','color','b');
xlabel('Distance from the midline, um','fontweight','b','color','r');
ylabel('Time of observation, min','fontweight','b','color','r');
strain_rate=gcf; 

#downloads the figure in the folder
print(strain_rate,strcat('strain_nb_of_cells=',num2str(av_nb_of_cells),'_delta_t=',num2str(delta_t),".pdf"))

#--------------------Stress 

figure(4)

imagesc(X_sorted,T,S_sorted);
colormap(jet(250));
c=colorbar;
ylabel(c,'Elastic Stress, Pa','fontsize',14)
set(gca,'xlim',[-600 600],'ylim',[0 800],'YDir','reverse');
title(strcat("Elastic Stress, Pa, HP=",num2str(HP)," , ode solver=",Type),'fontsize',18,'color','b');
xlabel('Distance from the midline, um','color','r');
ylabel('Time of observation, min','color','r');
stress=gcf;

#downloads the figure in the folder
print(stress,strcat('stress_nb_of_cells=',num2str(av_nb_of_cells),'_delta_t=',num2str(delta_t),".pdf"))

toc