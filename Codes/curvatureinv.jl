#Curvature invariants

Kretschmann(y,i)=5/3/exp(y[i,1])^2*(y[i,3]-y[i,4])^2*((-1 .+y[i,2])*y[i,3]-(1 .+y[i,2])*y[i,4])^2+1/(3*R[i]^4)/exp(y[i,1])^2*(6*exp(y[i,1])-2*R[i]^2*y[i,3]*(y[i,3]-y[i,4])+(1 .+y[i,2])*(-6 .+R[i]^2*(y[i,3]-y[i,4])^2));


Weyl(y,i)=1/(3*R[i]^4)/exp(y[i,1])^2*(6*exp(y[i,1])-2*R[i]^2*y[i,3]*(y[i,3]-y[i,4])+(1 .+y[i,2])*(-6 .+R[i]^2*(y[i,3]-y[i,4])^2));

#Der(y,i,k)=(y[i+1,k]-y[i-1,k])/(r[i+1]-r[i-1]);
DDer(y,i,k)=(y[i+1,k]-2*y[i,k]+y[i-1,k])/(R[i+1]-R[i-1])^2;

limit(y,i)=-4*y[i, 3]^2*y[i, 4]^2 + (y[i, 3]^2 + y[i, 4]^2)^2 + (-2*DDer(y,i,2)*(-1 + y[i, 2])*y[i, 2] - 4*Der(y,i,4)*R[i]*y[i, 2]*y[i, 3] - 2*y[i, 3]^2 + 5*y[i, 2]*y[i, 3]^2 + y[i, 2]^2*y[i, 3]^2 - 2*Der(y,i,4)*R[i]*y[i, 4] + 4*Der(y,i,4)*R[i]*y[i, 2]*y[i, 4] + 2*Der(y,i,4)*R[i]*y[i, 2]^2*y[i, 4] - 10*y[i, 2]*y[i, 3]*y[i, 4] - 2*y[i, 2]^2*y[i, 3]*y[i, 4] + 2*y[i, 4]^2 + 5*y[i, 2]*y[i, 4]^2 + y[i, 2]^2*y[i, 4]^2 - 2*Der(y,i,3)*R[i]*((1 - 4*y[i, 2] + y[i, 2]^2)*y[i, 3] + 2*y[i, 2]*y[i, 4]))^2/(4*(-1 + y[i, 2])^2);


Kretschv2(y,i)=((16*(-1 + (1 + y[i, 2])/exp(y[i,1]))^2)/R[i]^4 - (16*y[i, 3]*y[i, 4]*((-1 + y[i, 2])*y[i, 3] - y[i, 2]*y[i, 4])* (y[i, 2]*y[i, 3] - (1 + y[i, 2])*y[i, 4]))/exp(2*y[i,1]) +  (-2*DDer(y,i,2)*y[i, 2] + (1 + exp(y[i,1]) - 3*y[i, 2])*y[i, 3]^2 + 2*Der(y,i,4)*R[i]*(1 + y[i, 2])*y[i, 4] + (-3 + exp(y[i,1]) - 3*y[i, 2])*y[i, 4]^2 - 2*y[i, 3]*(Der(y,i,3)*R[i]*(-1 + y[i, 2]) + (-1 + exp(y[i,1]) - 3*y[i, 2])*y[i, 4]))^2/exp(2*y[i,1]) + (2*(2*(exp(y[i,1]) - R[i]^2*y[i, 3]^2) + (1 + y[i, 2])*(-2 + R[i]^2*(y[i, 3]^2 - y[i, 4]^2)))^2)/(exp(2*y[i,1])*R[i]^4) + (2*(2*(exp(y[i,1]) + R[i]^2*y[i, 3]^2) + (1 + y[i, 2])*(-2 + R[i]^2*(-y[i, 3]^2 + y[i, 4]^2)))^2)/(exp(2*y[i,1])*R[i]^4))/4;

function Kr(y,i)
        if i<4
                return limit(y,i)
        else
                return Kretschv2(y,i)
        end
end

function KrC(r)
	ind=Int(findall(abs.(R .-r) .<dx/2.0)[1])
	if ind<3 || ind>L-2
		return 0
	else 
		return	Kr(state_array,ind)
	end
end


RicciScalar(y,i)=1/exp(y[i,1])*(y[i,3]-y[i,4])*((-1 .+y[i,2])*y[i,3]-(1 .+y[i,2])*y[i,4]);


ADM(y)=(R1 .*(1 .-exp.(-y[3:L-2,1]./2)./sqrt.(1 .-y[3:L-2,2])))[length(R1)];
HypADM(y)=(r1 .*(1 .-exp.(-y[3:L-2,1]./2)./sqrt.(2 .-y[3:L-2,2])))[length(r1)];

function constraintcheck(g,i)
if i==3
	return 0
else
	return (g[i+1,2]-g[i-1,2])/(g[i+1,6]-g[i-1,6])-(exp(g[i,1])-(1+g[i,2]))/g[i,6]-4*pi*g[i,6]*(g[i,3]^2+g[i,4]^2)+(g[i+1,1]-g[i-1,1])/(g[i+1,6]-g[i-1,6])+(g[i,2])*((g[i+1,1]-g[i-1,1])/(g[i+1,6]-g[i-1,6])+4*pi*g[i,6]*(g[i,3]^2-g[i,4]^2));
end
end
