function [V,S,E_point]=Calc_v_s_e_HP2(T,Xi,n,f,initial_delta_x,av_nb_of_nodes,av_nb_of_cells,k,typical_surf)
  
V=zeros(size(T),av_nb_of_nodes);
S=zeros(size(T),av_nb_of_cells); 
E_point=zeros(size(T),av_nb_of_cells);

for t=1:size(T)

  #---values of Xi derivative/stress for first and last cells
  V(t,1)=(k*((Xi(t,2)-Xi(t,1)-initial_delta_x)/initial_delta_x)-f)/n;
  V(t,av_nb_of_nodes)=(-k*((Xi(t,av_nb_of_nodes)-Xi(t,av_nb_of_nodes-1)-initial_delta_x)/initial_delta_x)+f)/n;
  
if mod(av_nb_of_nodes,2)==0
     
     for i = 2:round(av_nb_of_nodes/2)
        V(t,i)=(-k/initial_delta_x*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+k/initial_delta_x*(Xi(t,i+1)-Xi(t,i)-initial_delta_x)-f)/n;
     endfor
     for i = round(av_nb_of_nodes/2)+2:av_nb_of_nodes-1
        V(t,i)=(-k/initial_delta_x*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+k/initial_delta_x*(Xi(t,i+1)-Xi(t,i)-initial_delta_x)+f)/n;
     endfor
  V(t,round(av_nb_of_nodes/2)+1)=(-k/initial_delta_x*(Xi(t,round(av_nb_of_nodes/2)+1)-Xi(t,round(av_nb_of_nodes/2))-initial_delta_x)+k/initial_delta_x*(Xi(t,round(av_nb_of_nodes/2)+2)-Xi(t,round(av_nb_of_nodes/2)+1)-initial_delta_x))/n;  
endif
   
   
if mod(av_nb_of_nodes,2)==1
    
     for i = 2:round(av_nb_of_nodes/2)-1
         V(t,i)=(-k/initial_delta_x*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+k/initial_delta_x*(Xi(t,i+1)-Xi(t,i)-initial_delta_x)-f)/n;
     endfor
     for i = round(av_nb_of_nodes/2)+2:av_nb_of_nodes-1
         V(t,i)=(-k/initial_delta_x*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+k/initial_delta_x*(Xi(t,i+1)-Xi(t,i)-initial_delta_x)+f)/n;
     endfor
  V(t,round(av_nb_of_nodes/2)+1)=(-k/initial_delta_x*(Xi(t,round(av_nb_of_nodes/2)+1)-Xi(t,round(av_nb_of_nodes/2))-initial_delta_x)+k/initial_delta_x*(Xi(t,round(av_nb_of_nodes/2)+2)-Xi(t,round(av_nb_of_nodes/2)+1)-initial_delta_x))/n;  
  V(t,round(av_nb_of_nodes/2))=(-k/initial_delta_x*(Xi(t,round(av_nb_of_nodes/2))-Xi(t,round(av_nb_of_nodes/2-1))-initial_delta_x)+k/initial_delta_x*(Xi(t,round(av_nb_of_nodes/2)+1)-Xi(t,round(av_nb_of_nodes/2))-initial_delta_x))/n;  
endif
   

   
    #---values of Xi derivative/strain rate/stress for other cells
for i = 2:av_nb_of_nodes-1
      E_point(t,i-1)=(V(t,i)-V(t,i-1))/(initial_delta_x*60);
      S(t,i-1)=k*((Xi(t,i)-Xi(t,i-1)-initial_delta_x)/initial_delta_x)*10^-9/(typical_surf*10^-12);
endfor 
  
  #----value of strain rate for last cell
E_point(t,av_nb_of_nodes-1)=(V(t,av_nb_of_nodes)-V(t,av_nb_of_nodes-1))/(initial_delta_x*60);
S(t,av_nb_of_nodes-1)=k*((Xi(t,av_nb_of_nodes)-Xi(t,av_nb_of_nodes-1)-initial_delta_x)/initial_delta_x)*10^-9/(typical_surf*10^-12);
   
end