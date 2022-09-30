#This file makes plots of max Kretschmann for all runs, its location in the domain and the mass of the horizon with respect to the ADM mass

include("./curvatureinv.jl");
include("./Crevolution.jl");
#using Plots
#plot()
using PyCall, PyPlot
const plt=PyPlot
using CSV, DataFrames, ProgressMeter

TAH1=0;

dir="./DATA"

using LaTeXStrings

for m in 1:1 #30:34
	AHloc=Matrix(DataFrame(CSV.File(joinpath(dir,"P_$m/AHlocP_$m.csv"),header=false, delim=",")));   #P_$m/AHloc.csv before
	TlocKmax=Matrix(DataFrame(CSV.File(joinpath(dir,"P_$m/TlocKmaxP_$m.csv"),header=false, delim=","))); #P_$m/TlocKmax.csv before
	len=length(TlocKmax[:,1]);
	t=AHloc[1:len,1];
	loc=TlocKmax[1:len,2];
	v=TlocKmax[1:len,3];
        #MBH=round(AHloc[length(AHloc[:,1]),3]/2, digits=4); #final BH mass
       	MBH=round(AHloc[findall(AHloc[:,3].>0)[1],3]/2, digits=2); #newly formed BH mass
	
	AH=findall(AHloc[:,3] .>0)[1];
	if TAH1==0
		global TAH1=t[AH]/MBH; #keeping the first one
	end
	dTAH=t[AH]/MBH-TAH1;
	
	tah=round(TAH1, digits=2); 
        if abs(dTAH)<10^(-7)
                plt.axvline(x= TAH1, label=L"t_{AH}/M_{AH}"*"=$tah", ls="--")
        end
        #plot!(AHloc[:,1]/MBH, AHloc[:,3]/MBH, labels="R_AH/M",xlabel="T/M", ylabel="R/M")
        
	### LOCATION OF MAXIMUM
	#plt.plot((t ./MBH .-dTAH), loc/MBH, ls="--", label=L"M_{AH}"*"=$MBH")
		plt.xlabel(L"t/M_{AH}", fontsize=16)
	#plt.ylabel(L"r/M_{AH}", fontsize=16)
	plt.xlim(0, TAH1/MBH+1)

        ####MAXIMUM KRETCHMANN
        
        #plt.plot((t ./MBH .-dTAH), 4 ./3 .*v.*MBH^4, xlims=(0,TAH1/MBH+7), labels=L"M_{AH}"*"= $MBH");
	
	#plt.xlabel(L"T/M_{AH}");
	plt.ylabel(L"K_1/K_{S}", fontsize=16);
	#plt.plot((t ./MBH .-dTAH)[1:2:length(v)], (4 ./3 .*v.*MBH^4)[1:2:length(v)], ls="--", label=L"M_{AH}"*"= $MBH")
	#plt.legend(loc=1)
	#plt.plot((t ./MBH)[1:2:length(t)], (4 ./3 .*v.*MBH^4)[1:2:length(t)], ls="--", label=L"M_{AH}"*"= $MBH")
	plt.plot((t ./MBH .-dTAH), 4 ./3 .*v.*MBH^4, ls="--", label=L"M_{AH}"*"= $MBH")
        plt.legend(loc=1)
	u=(4 ./3 .*v.*MBH^4)[1];
	println(u);
end


###ATRIBUTES
ax = plt.gca()
        for axis in ["top","bottom","left","right"]
        ax.spines[axis].set_linewidth(1.5)
        end
        ax.xaxis.set_tick_params(width=1, length=7, labelsize=12, direction="in")
        ax.yaxis.set_tick_params(width=1, length=7, labelsize=12, direction="in")
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(2))
        ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(2))
        ax.yaxis.set_ticks_position("both")
        ax.xaxis.set_ticks_position("both")
        ax.xaxis.set_tick_params(which="minor",width=1, length=4,direction="in")
        ax.yaxis.set_tick_params(which="minor",width=1, length=4,direction="in")

plt.savefig("./Plots/Kmaxout.png")
plt.show()

