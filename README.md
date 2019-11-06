# Parallel-Binary-Classification
Final Project - Parallel And Distributed Computing Course

Shalev Lazarof 

Parallel Solution for Binary Classification that combine MPI & OMP technologies


Sequential Algorithm 

1.	Set alpha = alpha0
2.	Choose initial value of all components of the vector of weights W to be equal to zero.
3.	Cycle through all given points Xi in the order as it is defined in the input file
4.	For each point Xi define a sign of discriminant function f(Xi) = WT Xi. 
    If the values of vector W is chosen properly then all points belonging to set A will have positive value of 
    f(Xi) and all points belonging to set B will have negative value of f(Xi).
    The first point P that does not satisfies this criterion will cause to stop the check and immediate redefinition of 
    the vector W: W = W + [alpha*sign(f(P))] P
5.	Loop through stages 3-4 till one of following satisfies:
a.	All given points are classified properly
b.	The number maximum iterations LIMIT is reached
6.	Find Nmis - the number of points that are wrongly classified, meaning that the value of f(Xi) is not correct 
    for those points. Calculate a Quality of Classifier q according the formula q = Nmis / N
7.	Check if the Quality of Classifier is reached (q is less than a given value QC). 
8.	Stop if q < QC.
9.	Increment the value of alpha:    alpha = alpha + alpha0.    
    Stop if alpha > alphaMAX
10.	Loop through stages 2-9
 
Parallel Algorithm

1.	Process 0 read from the data file send the data using MPI to the other slaves
2.	Each process initiates its alpha value using his rank.
3.	Each point is classified using OMP threads.
4.	Each result is checked, when finding the first wrong point classification the weights vector is calculated using OMP, 
5.	Repeat the process until the limit or all the points are correctly classified
6.	All the processes calculate wrong/n
7.	all processes Using MPI to gather all the results and check if qc was reached
8.	if one of the process reach qc/alpha max, all other processes is ending and the result send to the master
9.	the master process writes output file
 
 
why I choose OMP & MPI technologies

The OMP architecture
OpenMP is an implementation of multithreading, a method of parallelizing whereby a master thread forks a specified 
number of slave threads and the system divides a task among them. The threads then run concurrently, with the runtime 
environment allocating threads to different processors, there are no dependency between the factors and that’s why it 
is great solution.

The MPI architecture 
MPI is a communication protocol for programming parallel computers. Both point-to-point and collective communication 
are supported. MPI's goals are high performance, scalability, and portability, the transmission delay time is reduced
significantly because each process calculates different values on its own and handles the same amount of work, what 
lead to a good Load Balance.






