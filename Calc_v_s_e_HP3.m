
function [V,S,E_point]=Calc_v_s_e_HP3(T,Xi,n,f,initial_delta_x,av_nb_of_nodes,av_nb_of_cells,k,typical_surf,epsilon_s,flag_threshold_reached,T_eps_s)
   

global T_eps_s

#Initialisations of matrices / arrays
V=zeros(size(T),av_nb_of_nodes);
S=zeros(size(T),av_nb_of_cells); 
E_point=zeros(size(T),av_nb_of_cells);
flag_threshold_reached=zeros(av_nb_of_nodes,1);
   
for t =1:size(T) #for all times of the time array
  
     V(t,1)=((k/initial_delta_x)*(Xi(t,2)-Xi(t,1)-initial_delta_x)-f)/n;
     V(t,av_nb_of_nodes)=((-k/initial_delta_x)*(Xi(t,av_nb_of_nodes)-Xi(t,av_nb_of_nodes-1)-initial_delta_x)+f)/n;
     
     
     #---values of Xi derivative, self propelling force depending on strain of neighbouring cells
  
  if mod(av_nb_of_nodes,2)==0 #pair nb of nodes
    
     for i = 2:round(av_nb_of_nodes/2) #before the midline 
          if flag_threshold_reached(i)==1 #if strain threshold is already reached
             V(t,i)=(-f+(-k/initial_delta_x)*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+(k/initial_delta_x)*(Xi(t,i+1)-Xi(t,i)-initial_delta_x))/n;
             #self propelling force is exerced
         elseif flag_threshold_reached(i)==0  #if strain threshold is not reached yet
              if T(t)>=T_eps_s(i)  #if the current time is greater than the time when threshold is reached 
                  flag_threshold_reached(i)=1; #strain threshold gets reached
                  V(t,i)=(-f+(-k/initial_delta_x)*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+(k/initial_delta_x)*(Xi(t,i+1)-Xi(t,i)-initial_delta_x))/n;
                  #self propelling force is exerced
              else
                  V(t,i)=((-k/initial_delta_x)*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+(k/initial_delta_x)*(Xi(t,i+1)-Xi(t,i)-initial_delta_x))/n;
                  #self propelling force is not exerced
              endif
         endif
     endfor
       
       
     for i = round(av_nb_of_nodes/2)+1:av_nb_of_nodes-1 #after the midline
       
       #same behaviour as described before
       if flag_threshold_reached(i)==1
           V(t,i)=(f-(k/initial_delta_x)*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+(k/initial_delta_x)*(Xi(t,i+1)-Xi(t,i)-initial_delta_x))/n; 
       elseif flag_threshold_reached(i)==0
            if T(t)>=T_eps_s(i) 
                flag_threshold_reached(i)=1;
                V(t,i)=(f-(k/initial_delta_x)*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+(k/initial_delta_x)*(Xi(t,i+1)-Xi(t,i)-initial_delta_x))/n;
            else
                V(t,i)=(-(k/initial_delta_x)*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+(k/initial_delta_x)*(Xi(t,i+1)-Xi(t,i)-initial_delta_x))/n;
            endif
        endif
     endfor
       
       V(t,round(av_nb_of_nodes/2)+1)=(-(k/initial_delta_x)*(Xi(t,round(av_nb_of_nodes/2)+1)-Xi(t,round(av_nb_of_nodes/2))-initial_delta_x)+(k/initial_delta_x)*(Xi(t,round(av_nb_of_nodes/2)+2)-Xi(t,round(av_nb_of_nodes/2)+1)-initial_delta_x))/n;  
 endif


 #Same thing if nb of nodes is odd
 
 if mod(av_nb_of_nodes,2)==1
       for i = 2:round(av_nb_of_nodes/2)-1
          if flag_threshold_reached(i)==1
           V(t,i)=(-f +(-k/initial_delta_x)*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+(k/initial_delta_x)*(Xi(t,i+1)-Xi(t,i)-initial_delta_x))/n;
          elseif flag_threshold_reached(i)==0
              if T(t)>=T_eps_s(i) 
                  flag_threshold_reached(i)=1;
                  V(t,i)=(-f+(-k/initial_delta_x)*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+(k/initial_delta_x)*(Xi(t,i+1)-Xi(t,i)-initial_delta_x))/n;
              else
                  V(t,i)=((-k/initial_delta_x)*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+(k/initial_delta_x)*(Xi(t,i+1)-Xi(t,i)-initial_delta_x))/n;
              endif
         endif
         
     endfor
       
       
     for i = round(av_nb_of_nodes/2)+1:av_nb_of_nodes-1
       if flag_threshold_reached(i)==1
           V(t,i)=(f-(k/initial_delta_x)*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+(k/initial_delta_x)*(Xi(t,i+1)-Xi(t,i)-initial_delta_x))/n; 
        elseif flag_threshold_reached(i)==0
            if T(t)>=T_eps_s(i) 
                flag_threshold_reached(i)=1;
                V(t,i)=(f-(k/initial_delta_x)*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+(k/initial_delta_x)*(Xi(t,i+1)-Xi(t,i)-initial_delta_x))/n;
            else

                V(t,i)=(-(k/initial_delta_x)*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+(k/initial_delta_x)*(Xi(t,i+1)-Xi(t,i)-initial_delta_x))/n;
            endif
        endif
     endfor
     V(t,round(av_nb_of_nodes/2)+1)=(-k/initial_delta_x*(Xi(t,round(av_nb_of_nodes/2)+1)-Xi(t,round(av_nb_of_nodes/2))-initial_delta_x)+k/initial_delta_x*(Xi(t,round(av_nb_of_nodes/2)+2)-Xi(t,round(av_nb_of_nodes/2)+1)-initial_delta_x))/n;  
     V(t,round(av_nb_of_nodes/2))=(-k/initial_delta_x*(Xi(t,round(av_nb_of_nodes/2))-Xi(t,round(av_nb_of_nodes/2-1))-initial_delta_x)+k/initial_delta_x*(Xi(t,round(av_nb_of_nodes/2)+1)-Xi(t,round(av_nb_of_nodes/2))-initial_delta_x))/n;  

  endif
     
    
     
     #-----values of strain rate/stress for other cells
     for i = 2:av_nb_of_nodes-1
         E_point(t,i-1)=(V(t,i)-V(t,i-1))/(initial_delta_x*60);
         S(t,i)=(-k*((Xi(t,i)-Xi(t,i-1)-initial_delta_x)/initial_delta_x)+k*((Xi(t,i+1)-Xi(t,i)-initial_delta_x)/initial_delta_x))/typical_surf;
         S(t,i-1)=k*((Xi(t,i)-Xi(t,i-1)-initial_delta_x)/initial_delta_x)*10^-9/(typical_surf*10^-12);

     endfor
    
     #-----value of strain rate for the last cell
     E_point(t,av_nb_of_nodes-1)=(V(t,av_nb_of_nodes)-V(t,av_nb_of_nodes-1))/(initial_delta_x*60);
     S(t,av_nb_of_nodes-1)=k*((Xi(t,av_nb_of_nodes)-Xi(t,av_nb_of_nodes-1)-initial_delta_x)/initial_delta_x)*10^-9/(typical_surf*10^-12);

  endfor
end
