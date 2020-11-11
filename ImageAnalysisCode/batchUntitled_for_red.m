a=zeros(4,4,4);

a(1,1,:)=colocalization_of_preclusters_cellwise_Rconstant('Coordinates_2_25_3.dattt',1,1,50,130);
a(2,1,:)=colocalization_of_preclusters_cellwise_Rconstant('Coordinates_2_25_3.dattt',2,1,50,130);
a(3,1,:)=colocalization_of_preclusters_cellwise_Rconstant('Coordinates_2_25_3.dattt',3,1,50,130);
a(4,1,:)=colocalization_of_preclusters_cellwise_Rconstant('Coordinates_2_25_3.dattt',4,1,50,130);

a(1,2,:)=colocalization_of_preclusters_cellwise_Rconstant('Coordinates_3_25_3.dattt',1,1,50,130);
a(2,2,:)=colocalization_of_preclusters_cellwise_Rconstant('Coordinates_3_25_3.dattt',2,1,50,130);
a(3,2,:)=colocalization_of_preclusters_cellwise_Rconstant('Coordinates_3_25_3.dattt',3,1,50,130);
a(4,2,:)=colocalization_of_preclusters_cellwise_Rconstant('Coordinates_3_25_3.dattt',4,1,50,130);

a(1,3,:)=colocalization_of_preclusters_cellwise_Rconstant('Coordinates_4_25_3.dattt',1,1,50,130);
a(2,3,:)=colocalization_of_preclusters_cellwise_Rconstant('Coordinates_4_25_3.dattt',2,1,50,130);
a(3,3,:)=colocalization_of_preclusters_cellwise_Rconstant('Coordinates_4_25_3.dattt',3,1,50,130);
a(4,3,:)=colocalization_of_preclusters_cellwise_Rconstant('Coordinates_4_25_3.dattt',4,1,50,130);

a(1,4,:)=colocalization_of_preclusters_cellwise_Rconstant('Coordinates_5_25_3.dattt',1,1,50,130);
a(2,4,:)=colocalization_of_preclusters_cellwise_Rconstant('Coordinates_5_25_3.dattt',2,1,50,130);
a(3,4,:)=colocalization_of_preclusters_cellwise_Rconstant('Coordinates_5_25_3.dattt',3,1,50,130);
a(4,4,:)=colocalization_of_preclusters_cellwise_Rconstant('Coordinates_5_25_3.dattt',4,1,50,130);

round(a)