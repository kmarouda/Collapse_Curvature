#Constants for the initial data
a=2; #initial resolution 
s1,s2,s3=2,2,2; #gaussian initial data for the scalar field and its time derivative
R01,R02, R03=20.0,20.0,8; #position of the initial pulse
P=0.7; #0.16668   #0.52;
B=0;
P3=0;
D=0; #delta(R) wide gaussian centered away from the origin initial data
Rd=5;
sd=3;
#global epsilon=0.05; 0.05; #dissipation strength
global omicron=0; #dissipation parameter
ah=0; #adding excision set to 0 , without excision set to 2
#switch back to 0
cd()

println(pwd())

	include("./Codes/Crevolution.jl");
	include("./Codes/curvatureinv.jl");

using ProgressMeter
global ah=1; #turn off change in the RHS
global  AH=0;
global exc=0;
global tex=0;
ori=0.0; 
#Initial data , setting length of the domain and discretization
for m in 1:3
	dir="./DATA/Crconvergence/Cr_time_steps_$(m)"
	global ex=3;
	global epsilon=0.05; 
        global y0=1.0; #origin value of Cr variable 
	N=2.0^(m+a)*100.0;
	Rf=100.0;#1.0;#50.0
	global dx=Rf/N;
	global dt=dx*0.4;   
	Nt=100.0*2^(m+a);
        Tf=Nt*dt; #final time
        println("Resolution # ", m);
        println("the final time is ", Tf);
        println("the number of grid points and timesteps is ", N+1)

         #setting spatial and temporal discretization¶
	
	global R1=range(ori, stop=Rf, step=dx);
        global R=range(ori-2*dx, stop=Rf+2*dx, step=dx);
        println("step size is  ", R[2]-R[1])

	global L=length(R);
	
        global T=range(dt,stop=Tf, step=dt);
        println("the time step is ", T[2]-T[1])

	global l=length(T);  

        R1=zeros((L-4))
        for i in 1:length(R1)
                R1[i]=R[i+2]
        end

        # Solving for the constraint

        initial_data1=rungekutta4(CrconstraintRHS, y0, R1);

	initial_data=zeros(L);
        initial_data[3:L-2]=initial_data1[:] #adding the ghost positions

	#Defining the state array for the evolution
	global state_array=[delta(R) (initial_data .-1.0) scalar_timeder(R) scalar_spaceder(R) scalar_field(R)];
	state_array=ghost(state_array);
#Runge-Kutta evolution and writting the data in separate files
	AHloc=zeros(length(T),3);
	AHloc[:,1]=T;
		@showprogress "Evolution of resolution # $m " for k in 1:l
           			if length(filter(x-> x<0, 1.0 .+ state_array[:,2]))!=0 && ah!=2
            				if AH==0
         				global AH=k;
            				println(T[AH])
                			end
	    				aux=findall(1.0 .+state_array[:,2] .< 0)[length(findall(1.0 .+state_array[:,2].<0))];
	    				AHloc[k,2]=aux;
            				AHloc[k,3]=R[aux];
					if exc==0 && ah!=2 && length(filter(x-> x<0, 1.0 .+ state_array[:,2]))>0 
						global exc=(findall(1.0 .+state_array[:,2] .< 0)[length(findall(1.0 .+state_array[:,2].<0))])-2;
						global tex=k;
					println(tex);
					end
				end
					if tex!=0 && k==2^(m-1)*tex && ah!=2
                                        	if m==1
                                           		global ex=exc+2
                                           		#global epsilon=0.0;
                                        		println("excision starts now at R=$(R[ex]) and T=$(T[k])");
                                        	elseif m==2
                                                	global ex=2*exc+1
                                                	#global epsilon=0.0;
                                        		println("excision starts now at R=$(R[ex]) and T=$(T[k])");
                                        	elseif m==3
                                                	global ex=4*exc-1
                                                	#global epsilon=0.0;
                                        		println("excision starts now at R=$(R[ex]) and T=$(T[k])");
                                        	end
					end
			if ex!=3
    			state_array=rungekutta4molstep(CRRHSexcision,state_array,T,k,ex); 
    			state_array=extrapolation(state_array,ex);
    			else
    			state_array=rungekutta4molstep(CRRHS,state_array,T,k,3)
    			end    
    			state_array=ghost(state_array);
			using CSV, Tables
				if mod(k,2^(m-1))==0
					CSV.write(joinpath(dir,"time_step$k.csv"), Tables.table([state_array KrC.(R) R]), writeheader=false)
        			end
			end
	using CSV, Tables
    	CSV.write(joinpath(dir,"AHloc.csv"), Tables.table(AHloc), writeheader=false)
end


