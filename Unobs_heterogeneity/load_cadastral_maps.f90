subroutine load_cadastral_maps()
    use cadastral_maps
    implicit none
    character(LEN=1)::s_c1
    character(LEN=2)::s_c2
    integer::v_l,i,j,ind!,x,y
    double precision::u
    integer,dimension(villages,plots_in_map)::area_type
    
    PA_type=-9
        
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

        !Load area type
        OPEN(UNIT=12, FILE=file_map//"area_type.csv")
            read(12,*),area_type
        close(12)

        PA_type(1:plots_v(v_l),2,v_l)=area_type(v_l,1:plots_v(v_l))
    end do

end subroutine