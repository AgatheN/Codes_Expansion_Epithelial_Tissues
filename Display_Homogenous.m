        
function [V_sorted,S_sorted,E_point_sorted]=Display_Homogenous(Xi,X_sorted,idx_cells,idx_nodes,size_T,V,E_point,S,V_sorted,E_point_sorted,S_sorted)
     
  
   for t=1:size_T
        [n_xi, idx_xi] = histc(Xi(t,:),X_sorted); #Sorts the cells' positions in the intervals 
        #defined par X_sorted
        for i=1:size(idx_xi,2)-1 #For all cells
                   V_sorted(t,idx_xi(i):idx_xi(i+1))=V(t,idx_nodes(i));
                   E_point_sorted(t,idx_xi(i):idx_xi(i+1))=E_point(t,idx_cells(i)); 
                   S_sorted(t,idx_xi(i):idx_xi(i+1))=S(t,idx_nodes(i)); 
        endfor
        V_sorted(t,idx_xi(end):(idx_xi(end)+n_xi(end)-n_xi(end-1)))=V(t,end);
        S_sorted(t,idx_xi(end):(idx_xi(end)+n_xi(end)-n_xi(end-1)))=S(t,end);
   endfor

end