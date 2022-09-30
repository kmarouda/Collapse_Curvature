print("Which timestep would you like to see in the pointwise convergence? \n\n")
time= readline();
time = parse(Float64, time);
cd()
println(pwd())

#producing the plots with respect to the 1st resolution
m=1;
M=1;
a=2;
y0=1; #origin value of the Misner-Sharp Mass function appearing in the constraint
N=2^(m+a)*100;
Rf=100.0; #1.0; #50.0 
dx=Rf/N;
dt=dx*0.4;
Nt=100*2^(m+a)-1;
Tf=Nt*dt; #final time
ori=0.0; #2.0*M/(1.0+2.0*M); #2/3; #1.9;

global R1=range(ori, stop=Rf, step=dx);
R =range(ori-2*dx, stop=Rf+2*dx, step=dx);
L=length(R);

T=range(dt,stop=Tf,step=dt)

###FINDING the timestep of the requested time
if time!=0
timestep=findall( abs.(T .-time) .<dt)[1];
end

include("./Codes/curvatureinv.jl");
include("./Codes/Crevolution.jl");



dir1="./DATA/Crconvergence/Cr_time_steps_1"
dir2="./DATA/Crconvergence/Cr_time_steps_2"
dir3="./DATA/Crconvergence/Cr_time_steps_3"


using CSV, DataFrames, ProgressMeter
AHloc=Matrix(DataFrame(CSV.File(joinpath(dir1,"AHloc.csv"),header=false, delim=",")));
MBH=round(AHloc[length(T),3]/2, digits=3);
global aux=0; 
println(aux);

#pointwise convergence
using DataFrames,CSV
function c(n)
        diff=zeros(12,(length(R)-4));
            k=Int(2*n);
            l=Int(4*n);
            res1=DataFrame(CSV.File(joinpath(dir1,"time_step$n.csv"),header=false, delim=","))
            res2=DataFrame(CSV.File(joinpath(dir2,"time_step$k.csv"),header=false, delim=","))
            res3=DataFrame(CSV.File(joinpath(dir3,"time_step$l.csv"),header=false, delim=","))
            for i in 1:L-4 #findall(abs.(R .-10.0) .< dx)[1] 
		if i>Int(aux)
                #D
                diff[1,i]=R[i+2]*sqrt((res1.Column1[i+2]-res2.Column1[2*i+1])^2);
                diff[2,i]=R[i+2]*sqrt((res3.Column1[4*i-1]-res2.Column1[2*i+1])^2);
                #Cr
                diff[3,i]=R[i+2]*sqrt((res1.Column2[i+2]-res2.Column2[2*i+1])^2);
                diff[4,i]=R[i+2]*sqrt((res3.Column2[4*i-1]-res2.Column2[2*i+1])^2);
                #Pi
                diff[5,i]=R[i+2]*sqrt((res1.Column3[i+2]-res2.Column3[2*i+1])^2);
                diff[6,i]=R[i+2]*sqrt((res3.Column3[4*i-1]-res2.Column3[2*i+1])^2);
                #Phi
                diff[7,i]=R[i+2]*sqrt((res1.Column4[i+2]-res2.Column4[2*i+1])^2);
                diff[8,i]=R[i+2]*sqrt((res3.Column4[4*i-1]-res2.Column4[2*i+1])^2);
		#constraint
		diff[11,i]=R[i+2]*sqrt((constraintcheck(Matrix(res1),i+2)-constraintcheck(Matrix(res2),2*i+1))^2);
		diff[12,i]=R[i+2]*sqrt((constraintcheck(Matrix(res3),4*i-1)-constraintcheck(Matrix(res2),2*i+1))^2);
            	#R*phi
                diff[9,i]=R[i+2]*sqrt((res1.Column5[i+2]-res2.Column5[2*i+1])^2);
                diff[10,i]=R[i+2]*sqrt((res3.Column5[4*i-1]-res2.Column5[2*i+1])^2);
		end
	     end
            return diff
end

#L2 Norm convergence
using DataFrames,CSV
function normc(n)
        norm1=0;
        norm2=0;
        norm3=0;
        norm4=0;
        norm5=0;
        norm6=0;
        norm7=0;
        norm8=0;
	norm9=0;
	norm10=0;
	norm11=0;
	norm12=0;
	norm13=0;
	norm14=0;
            k=Int(2*n);
            l=Int(4*n);
            res1=DataFrame(CSV.File(joinpath(dir1,"time_step$n.csv"),header=false, delim=","))
            res2=DataFrame(CSV.File(joinpath(dir2,"time_step$k.csv"),header=false, delim=","))
            res3=DataFrame(CSV.File(joinpath(dir3,"time_step$l.csv"),header=false, delim=","))
	    for i in 1:L-4    #(Int(floor(length(R)/2)))
	        #if i>Int(aux)
                # delta
                norm1=norm1+R[i+2]^2*(res1.Column1[i+2]-res2.Column1[2*i+1])^2;
                norm2=norm2+R[i+2]^2*(res3.Column1[4*i-1]-res2.Column1[2*i+1])^2;
                # Cr
                norm3=norm3+R[i+2]^2*(res1.Column2[i+2]-res2.Column2[2*i+1])^2;
                norm4=norm4+R[i+2]^2*(res3.Column2[4*i-1]-res2.Column2[2*i+1])^2;
                # Pi
                norm5=norm5+R[i+2]^2*(res1.Column3[i+2]-res2.Column3[2*i+1])^2;
                norm6=norm6+R[i+2]^2*(res3.Column3[4*i-1]-res2.Column3[2*i+1])^2;
                # Phi
                norm7=norm7+R[i+2]^2*(res1.Column4[i+2]-res2.Column4[2*i+1])^2;
                norm8=norm8+R[i+2]^2*(res3.Column4[4*i-1]-res2.Column4[2*i+1])^2;
		# phi
                norm11=norm11+R[i+2]^2*(res1.Column5[i+2]-res2.Column5[2*i+1])^2;
                norm12=norm12+R[i+2]^2*(res3.Column5[4*i-1]-res2.Column5[2*i+1])^2;
		# Kretsch
		norm13=norm13+R[i+2]^2*(res3.Column6[4*i-1]-res2.Column6[2*i+1])^2;
		norm14=norm14+R[i+2]^2*(res3.Column6[4*i-1]-res2.Column6[2*i+1])^2;
		if i>Int(aux+1)
		#Constraint check
			norm9=norm9+R[i+2]^2*(constraintcheck(Matrix(res1),i+2)-constraintcheck(Matrix(res2),2*i+1))^2;
			norm10=norm10+R[i+2]^2*(constraintcheck(Matrix(res3),4*i-1)-constraintcheck(Matrix(res2),2*i+1))^2;
            	end
	     end
	sum1=norm1+norm3+norm5+norm7+norm11;
	sum2=norm2+norm4+norm6+norm8+norm12;
        quotient=[sqrt(norm1)/sqrt(norm2) sqrt(norm3)/sqrt(norm4) sqrt(norm5)/sqrt(norm6) sqrt(norm7)/sqrt(norm8) sqrt(norm9)/sqrt(norm10) sqrt(norm11)/sqrt(norm12) sqrt(norm13)/sqrt(norm14) sqrt(sum1/sum2)]
    return quotient
end

print("\n \n Plots are being processed")
if time!=0

using Plots, ProgressMeter
#for n in 1:20:length(T)
plot_array = Any[]
legend=["delta" "Cr" "Pi" "Phi" "Constraint"];
#legend=["delta" "Cr" "Psi" "Psibar" "R*phi"];
n=timestep;
plot()
x=R[3:L-2];
y=zeros((L-4));
z=zeros((L-4));
 @showprogress "Producing pointwise convergence plot..." for k in 1:5
    for i in 1:(length(R)-4);
        y[i]=c(n)[1+2*(k-1),i];
        z[i]=4*c(n)[2+2*(k-1),i];
    end
    push!(plot_array,plot(x, [y z], marker=:cross,xlabel="R",ylabel=legend[k]))
end
d=round(T[n],digits=2)
fig=plot!(plot_array[1],plot_array[2], plot_array[3], plot_array[4], plot_array[5], title="T=$d")
#fig=plot!(plot_array[1],plot_array[2], plot_array[3], title="T=$d")
savefig(fig,"./Plots/pointwise$n.png");
end

print("\n the pointwise convergence plot has been produced properly for all functions\n")


if true
using Plots, ProgressMeter
D=zeros((length(T)));
M=zeros((length(T)));
Phi=zeros((length(T)));
Pi=zeros((length(T)));
phi=zeros((length(T)));
con=zeros((length(T)));
full=zeros(length(T));
KK=zeros(length(T));
plot()
@showprogress "Procuding Norm convergence plot..." for n in 1:length(T)
        #D[n]=log(2,normc(n)[1]); #log with base two since r=2 in our case. we are doubling the precision
        #M[n]=log(2,normc(n)[2]);
        #Pi[n]=log(2,normc(n)[3]);
        #Phi[n]=log(2,normc(n)[4]);
	#con[n]=log(2,normc(n)[5]);
	#phi[n]=log(2,normc(n)[6]);
	#full[n]=log(2,normc(n)[7]);
	KK[n]=log(2,normc(n)[8]); #Kretschmann convergence up to collapse
end
#plot!(T, full, labels="full",ylims=(0,4))
#plot!(T, [Pi Phi], labels=["Pi" "Phi"],ylims=(0,5))
#plot!(T, [full D M Pi Phi], labels=["full" "delta" "Cr" "Pi" "Phi"],ylims=(0,5)) #comment out for norm convergence of different variables or of the full solution
plot!(T, KK, labels="Kretschmann", ylims=(0,4))

#plot!(T, [Pi Phi phi full], labels=["Psi" "Psibar" "R*phi" "full"], ylims=(0,4))
fig1=plot!(title="L2 Norm convergence",xlabel="T",ylabel="log_2_(c(T))")
savefig(fig1,"./Plots/Krnormconvergence$a.png")
print("\n the norm convergence plot has been produced properly\n")
end

