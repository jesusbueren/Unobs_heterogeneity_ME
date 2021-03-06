subroutine load_cadastral_maps()
    use cadastral_maps; use primitives
    implicit none
    character(LEN=1)::s_c1
    character(LEN=2)::s_c2
    integer::v_l,i,j,ind,a_l!,x,y
    double precision::u
    double precision,dimension(villages,plots_in_map)::areas
    integer,dimension(villages,plots_in_map)::area_type
    integer,dimension(P_max,2)::PA_stat
    
    PA_type=-9
    
    OPEN(UNIT=12, FILE=file_map//"area_type.csv")
        read(12,*),areas
    close(12)
       
    !Load map (who is connected to who)
    neighbors_map=0
    do v_l=1,villages
        if (v_l<10) then
            Write( s_c1, '(I1)' )  v_l
            OPEN(UNIT=12, FILE=file_map//"map_"//s_c1//".csv")
        else
            Write( s_c2, '(I2)' )  v_l
            OPEN(UNIT=12, FILE=file_map//"map_"//s_c2//".csv")
        end if
            read(12,*),neighbors_map(1:plots_v(v_l),1:plots_v(v_l),v_l)
        close(12)
        
        !Store which are my neighbors
        do i=1,plots_in_map;
            neighbors_map(i,i,v_l)=1
            ind=0
            do j=1,plots_in_map
                if (neighbors_map(i,j,v_l)==1 .and. ind<P_max) then !the number of neighbors cannot be greater than P_max ortherwise I select the firt P_max neighbors
                    ind=ind+1
                    neighbors(i,ind,v_l)=j
                end if
            end do
            !I need to make sure that the reference plot is a neighbor plot for plots whose
            !number of neighboring plots is larger than P_max in the data so in case it is not, I force the last neighbor to be the reference plot
            if (ind==P_max .and. ind<i) then
                neighbors(i,P_max,v_l)=i
            end if
                
        end do
    
        !number of neighbors
        PA_type(1:plots_v(v_l),1,v_l)=min(sum(neighbors_map(1:plots_v(v_l),1:plots_v(v_l),v_l),2),P_max)
                
    end do
    
    !Area Type
    do v_l=1,villages
        do i=1,plots_v(v_l)
            do a_l=1,types_a-1
                if (areas(v_l,i)<=area_lims(a_l)) then
                    area_type(v_l,i)=a_l
                    exit
                else
                    area_type(v_l,i)=a_l+1
                end if
            end do
            PA_type(i,2,v_l)=area_type(v_l,i)
        end do 
    end do 
    
    v_l=1
    PA_stat=0
    do i=1,plots_v(v_l)
        PA_stat(PA_type(i,1,v_l),1)=PA_stat(PA_type(i,1,v_l),1)+1
        PA_stat(PA_type(i,2,v_l),2)=PA_stat(PA_type(i,2,v_l),2)+1
    end do
    print*,'distribution of number of neighbors'
    print*,dble(PA_stat(1:P_max,1))/dble(plots_v(v_l))
    print*,'distribution of area types'
    print*,dble(PA_stat(1:types_a,2))/dble(plots_v(v_l))

end subroutine