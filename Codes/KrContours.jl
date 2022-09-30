m=5;
N=2^m*100;
Rf=100;
dx=Rf/N;
R1=range(0, stop=Rf, step=dx);
dt=dx*0.4;
Nt=100*2^m;
Tf=Nt*dt; #final time

R =range(-2*dx, stop=Rf+2*dx, step=dx);
L=length(R);
println("step size is  ", dx)

T=range(dt,stop=Tf,step=dt)
println("the time step is ", dt)

include("./Codes/Crevolution.jl");
include("/Codes/curvatureinv.jl");

dir="./DATA/Kretschmann/P_1"

using CSV, DataFrames, LaTeXStrings
AHloc=Matrix(DataFrame(CSV.File(joinpath(dir,"AHloc.csv"),header=false, delim=",")));
MADM=ADM(Matrix(DataFrame(CSV.File(joinpath(dir,"time_step1.csv"),header=false, delim=","))));

inds=1:Int(floor(length(T)/2.0));

MBH=AHloc[length(AHloc[inds,1]),3]/2;

MAH=AHloc[findall(AHloc[inds,3] .>0)[1],3]/2.0;
println(MAH);

r=R1[1:800]./MAH;
t=AHloc[inds,1]./MAH;

#dir="/home/kri/Documents/DATA/KretschmannO4/"
using LaTeXStrings

using CSV, DataFrames, ProgressMeter
Z=-1.5 .*ones(length(r),length(t));#length(T)) #neglecting the ghost points and having the same length for T as R

global n=1;
@showprogress for k in inds
    functions3=DataFrame(CSV.File(joinpath(dir,"time_step$k.csv"),header=false, delim=","));
    for i in 1:length(r) #length(R1)
        if AHloc[n,3]!=0 && i>Int(AHloc[n,2]+2)
            val1=4.0/3.0*MAH^(4.0)*Kr(Matrix(functions3),i+2);
            #val2=MAH^(2.0)*RicciScalar(Matrix(functions3),i+2);
	    #val=R[i+2]*(functions3.Column3[i+2]+(1.0 +functions3.Column2[i+2])/(1.0-functions3.Column2[i+2])*functions3.Column4[i+2]);     #MADM*functions3.Column5[i+2];
            Z[i,n]=val1; #R[i+2]*functions3.Column5[i+2];#MADM^4*Kr(Matrix(functions3),i+2);
        elseif AHloc[n,3]==0
            val1=4.0/3.0*MAH^(4.0)*Kr(Matrix(functions3),i+2);
            #val2=MAH^(2.0)*RicciScalar(Matrix(functions3),i+2);
	    #val=R[i+2]*(functions3.Column3[i+2]+(1.0 +functions3.Column2[i+2])/(1.0-functions3.Column2[i+2])*functions3.Column4[i+2]); #MADM*functions3.Column5[i+2];
            Z[i,n]=val1; #log(abs(val)); #R[i+2]*functions3.Column5[i+2];#MADM^4*Kr(Matrix(functions3),i+2);
        end
    end
	global  n=n+1;
end



#need to activate venv for this plot
using PlotlyJS, WebIO, LaTeXStrings
colorscale = [[0, "lightsalmon"], [0.5, "mediumturqoise"], [1, "gold"]];
#IJulia.clear_output(true)
#plot(AHloc[:,3], AHloc[:,1])
fig=PlotlyJS.plot(contour(
        
    x=r, #r horizontal axis

    y=t, # t , or -log.(t .-Tauo), #T vertical axis

    z=Z' ,colorscale="Hot",

    contours=attr(

        coloring ="Heatmap",

        showlabels = true, # show labels on contours

        labelfont = attr( # label font properties

            size = 12,

            color = "white", 
        )),
    
    colorbar=attr(

        title="K1/KS", # title here

        titleside="right",

        titlefont=attr(

            size=14,

            family="computer modern"

        )

    ),

), Layout(

    #title="Kretschmann invariant, excised BH evolution",

    xaxis_title="r/M",

    yaxis_title="t/M",

    legend_title="K1/KS",

    font=attr(

        family="computer modern",

        size=18,

        color="Black"

    )))


#using SyncPlot
savefig(fig ::Union{PlotlyJS.SyncPlot, Plot},"./Plots/ContourK.png");
