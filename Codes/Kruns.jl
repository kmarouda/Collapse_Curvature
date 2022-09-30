#Initial data constants definition¶
s1,s2=2,2;
R01,R02,R03=20,20,10; 
global P=2.1; 
P3=0;
B=0;
D=0;
Rd=5;
sd=2;

w=0.0; #high frequency initial data, BUT keep it to zero for gaussian initial data

#Initial data , setting length of the domain and discretization¶
m=4;
y0=1.0;  #Origin value of the Cr function appearing in the constraint
N=2^m*100;
Rf=100;
dx=Rf/N;
R1=range(0, stop=Rf, step=dx);
dt=dx*0.4;
Nt=4*2^m-1;
Tf=Nt*dt; #final time
println("the final time is ", Tf)


epsilon=0.05; #dissipation strength

R =range(-2*dx, stop=Rf+2*dx, step=dx); #change in case you want to evolve a symmetric domain

L=length(R);
println("step size is  ", dx)

T=range(dt,stop=Tf,step=dt)
println("the time step is ", dt)

include("./Codes/curvatureinv.jl");
include("./Codes/Crevolution.jl");

dir = "./DATA/Kretschmann"

using ProgressMeter

for m in 1:1
global P=P+0.05;
initial_data=rungekutta4(CrconstraintRHS, y0, R1);
initial_data1=zeros(L);
initial_data1[3:L-2]=initial_data[:];

#Defining the state array for the evolution
state_array=[delta(R) (initial_data1 .-1) scalar_timeder(R) scalar_spaceder(R) scalar_field(R)];
state_array=ghost(state_array);

ah=0; #adding excision set to 0 , without excision set to 2
global AH=0;
global ex=0;
AHloc=zeros(length(T),3);
AHloc[:,1]=T;
	@showprogress for k in 1:length(T)
    		if length(filter(x-> x<0, 1 .+ state_array[:,2]))!=0 && ah!=2
            		aux=findall(1 .+state_array[:,2] .< 0)[length(findall(1 .+state_array[:,2].<0))]; #.-2;
            		AHloc[k,2]=aux;
            		AHloc[k,3]=R[aux];
			if AH==0
  		        global AH=k;
                        println(T[AH])
			end
		end
            	if ex==0 && length(filter(x-> x<0, 1 .+ state_array[:,2]))>50
			global ex=findall(1 .+state_array[:,2] .< 0)[length(findall(1 .+state_array[:,2].<0))].-4;
              	end
		if ex!=0
    		state_array=rungekutta4molstep(CRRHSexcision,state_array,T,k,ex)
    		state_array=extrapolation(state_array,ex)
    		else
    		state_array=rungekutta4molstep(CRRHS,state_array,T,k,0)
    		end
    		state_array=ghost(state_array)
		#if mod(k,5)==1  #saving less timesteps during the evolution to save space
    			using CSV, Tables
    			CSV.write(joinpath(dir, "P_$m/time_step$k.csv"), Tables.table(state_array), writeheader=false)
		#end
	end

#writing down the position of the AH in a file
	using CSV, Tables, Peaks
    		CSV.write(joinpath(dir, "P_$m/ADMmass.csv"), Tables.table([ADM(state_array)]), writeheader=false);

#writing down the position of the AH in a file
	using CSV, Tables
    		CSV.write(joinpath(dir, "P_$m/AHloc.csv"), Tables.table(AHloc), writeheader=false)


	endtime=length(T)
	v=zeros(endtime);
	t=zeros(length(v))
	loc=zeros(length(v))
	u=zeros(length(v));
	locR=zeros(length(v));
	using ProgressMeter, CSV, DataFrames, Peaks, Polynomials
	@showprogress for n in 1:endtime
    		K=zeros(length(R)-4);
		Ric=zeros(length(R)-4);
    		functions3=DataFrame(CSV.File(joinpath(dir,"P_$m/time_step$n.csv") , header=false, delim=","));
    		for i in Int(AHloc[n,2]+3):(length(R)-2) #only looking for Kmax before the horizon forms
        		K[i-2]=Kr(Matrix(functions3),i); 
    			Ric[i-2]=RicciScalar(Matrix(functions3),i);
		end
    		t[n]=T[n]
    		v[n]=findmax(K)[1];
    		loc[n]=R[findmax(K)[2]+2];
		u[n]=findmax(Ric)[1];
		locR[n]=R[findmax(Ric)[2]+2];
	end

	using CSV, Tables
	    CSV.write(joinpath(dir,"P_$m/TlocKmax.csv"), Tables.table([t loc v]), writeheader=false)
            CSV.write(joinpath(dir,"P_$m/TlocRmax.csv"), Tables.table([t locR u]), writeheader=false)

	MADM=ADM(state_array);
	MAH=AHloc[findall(AHloc[:,3] .>0)[1],3];
	MBH=round(AHloc[length(AHloc[:,1]),3]/2, digits=3);
	P1=round(P, digits=3);

            CSV.write(joinpath(dir,"P_$m/Mass.csv"), Tables.table([P1 s2 MADM MAH MBH]), writeheader=false)
end

