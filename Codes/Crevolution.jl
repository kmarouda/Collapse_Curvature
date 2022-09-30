# Definition of gaussian initial data functions

function scalar_timeder(R)
    n=length(R);
    if n==1
        z=(P* (-((2 *exp(-((R - R02)^2/s2^2)) *(R - R02))/s2^2) + (2*exp(-((R + R02)^2/s2^2))* (R + R02))/s2^2))*cos(w*R)/(sqrt(2*pi)* s2)/R-(P* (exp(-((R - R02)/s2)^2) + exp(-((R + R02)/s2)^2))/(sqrt(2*pi)*s2))*w*sin(w*R)/R;
    else
    z=zeros(n);
    for i in 1:n
        z[i]=(P* (-((2 *exp(-((R[i] - R02)^2/s2^2)) *(R[i] - R02))/s2^2) + (2*exp(-((R[i] + R02)^2/s2^2))* (R[i] + R02))/s2^2))*cos(w*R[i])/(sqrt(2*pi)* s2)/R[i]-(P* (exp(-((R[i] - R02)/s2)^2) + exp(-((R[i] + R02)/s2)^2))/(sqrt(2*pi)*s2))*w*sin(w*R[i])/R[i];
    end
    end
    return z

end

function scalar_spaceder(R)
    n=length(R);
    if n==1
        z=(P* (-((2 *exp(-((R - R02)^2/s2^2)) *(R - R02))/s2^2) - (2*exp(-((R + R02)^2/s2^2))* (R + R02))/s2^2))/(sqrt(2*pi)* s2)*cos(w*R)/R-(P* (exp(-((R - R02)/s2)^2) + exp(-((R + R02)/s2)^2))/(sqrt(2*pi)*s2))*w*sin(w*R)/R-(P* (exp(-((R - R02)/s2)^2) + exp(-((R + R02)/s2)^2))/(sqrt(2*pi)*s2))*cos(w*R)/R^2.0; 
    else
    z=zeros(n);
    for i in 1:n
        z[i]=(P* (-((2 *exp(-((R[i] - R02)^2/s2^2)) *(R[i] - R02))/s2^2) - (2*exp(-((R[i] + R02)^2/s2^2))* (R[i] + R02))/s2^2))*cos(w*R[i])/(sqrt(2*pi)* s2)/R[i]-(P* (exp(-((R[i] - R02)/s2)^2) + exp(-((R[i] + R02)/s2)^2))/(sqrt(2*pi)*s2))*w*sin(w*R[i])/R[i]-(P* (exp(-((R[i] - R02)/s2)^2) + exp(-((R[i] + R02)/s2)^2))/(sqrt(2*pi)*s2))*cos(w*R[i])/R[i]^2.0; 
    end
    end
    return z
end

function scalar_field(R)
    n=length(R);
    if n==1
        z= (P* (exp(-((R - R02)/s2)^2) + exp(-((R + R02)/s2)^2))*cos(w*R)/(sqrt(2*pi)*s2))/R
    else
        z=zeros(n);
        for i in 1:n
            z[i]= (P* (exp(-((R[i] - R02)/s2)^2) + exp(-((R[i] + R02)/s2)^2))*cos(w*R[i])/(sqrt(2*pi)*s2))/R[i]
        end
    end
    return z
end

#ghosts

function ghost(y)
    L=length(y[:,1])
    y[1,1]=y[5,1]; #delta is even
    y[2,1]=y[4,1];
    y[1,2]=y[5,2]; #Cr is even
    y[2,2]=y[4,2];
    y[1,3]=y[5,3]; #Pi is even
    y[2,3]=y[4,3];
    y[1,4]=-y[5,4]; # Phi is odd
    y[2,4]=-y[4,4];
    y[1,5]=y[5,5]; #phi is even
    y[2,5]=y[4,5];    
    y[L-1,:]=2*y[L-2,:]-y[L-3,:];
    y[L,:]=2*y[L-1,:]-y[L-2,:];
    return y
end

#fourth order  dissipation

function dissipation4(y,i)
	delta4=(y[i+2,:]-4*y[i+1,:]+6*y[i,:]-4*y[i-1,:]+y[i-2,:]);
    return (-1)^4*epsilon*1/(dx)*delta4
end

#upwinded dissipation operator at the excision surface, Commented Out 
function updissipationN(y,i) #applied to i=ex
	delta4=(y[i+4,:]-4*y[i+3,:]+6*y[i+2,:]-4*y[i+1,:]+y[i,:]);
    return epsilon*1/(dx)*delta4
end

#upwinded dissipation operator 1 gridpoint away from the excision surface (in the computational domain), Commented Out
function updissipationNN(y,i) #i=ex+1
	delta4=(y[i+3,:]-4*y[i+2,:]+6*y[i+1,:]-4*y[i,:]+y[i-1,:]);
    return epsilon*1/(dx)*delta4
end

# Definition of the RHS of the equation to be integrated for the initial data

Crconstraint(y,R1)=1/R1-y/R1-8*pi*R1*scalar_timeder(R1)^2+y*(4*pi*R1*(scalar_timeder(R1)^2-scalar_spaceder(R1)^2));

function CrconstraintRHS(y,R1)
        if R1<10^(-7)
                return 0
        else
                return Crconstraint(y,R1)
        end
end

#Building initial data with a  Runge-Kutta integrator for the constraint

function rungekutta4(f,y0::Float64,T)
    n = length(T)
    y = zeros((n, length(y0)))
    y[1] = y0;
    for i in 1:n-1
        h = T[2] - T[1]
        k1 = f(y[i], T[i])
        k2 = f(y[i] + k1 * h/2, T[i] + h/2)
        k3 = f(y[i] + k2 * h/2, T[i] + h/2)
        k4 = f(y[i] + k3 * h, T[i] + h)
        y[i+1] = y[i] + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
    end
    return y
end


# Runge Kutta integrator used for the method of lines

function rungekutta4molstep(f,y0,T,w::Int64,ex)
    y = zeros((length(R),length(y0[1,:])));
    y[:,:] = y0[:,:]
        h = dt; # only for equally spaced grid in time and space, otherwise (T[w+1] - T[w])
	k1 = f(y[:,:], T[w])
        k1=ghost(k1)
        k2 = f(y[:,:] + k1 * h/2, T[w] + h/2)
        k2=ghost(k2)
        k3 = f(y[:,:] + k2 * h/2, T[w] + h/2)
        k3=ghost(k3)
        k4 = f(y[:,:] + k3 * h, T[w] + h)
        k4=ghost(k4)
        y[:,:] = y[:,:] + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
    return ghost(y[:,:])
end

function boundary(y,i)
    dk=zeros(length(y[1,:]));
    dk[1]=DeltaRHS(y,i);
    dk[2]=CrRHS(y,i);
    dk[3]=y[L-2,2]/(-1 .+y[i,2])*(Der(y,i,3)+y[i,3]/R[i])-1/(-1 .+y[i,2])^2*dk[2]*(Der(y,i,5)+y[i,5]/R[i]);
    dk[4]=Der(y,i,3);
    dk[5]=y[i,3];
    return dk-dissipation4(y,i)
end

###Setting up the functions for the bulk of different rescaling trials

Divergent_terms(y,i)=3*((y[i+1,3]*(y[i+1,2])-y[i+1,4]*(1.0 .+y[i+1,2]))*R[i+1]^2-(y[i-1,3]*(y[i-1,2])-y[i-1,4]*(1.0 .+y[i-1,2]))*R[i-1]^2)/(R[i+1]^3-R[i-1]^3);

Good_terms(y,i)=(y[i,2])*(Der(y,i,3))-(y[i,3]-y[i,4])*CrRHS(y,i);


ScalarRHSCr(y,i)=(Divergent_terms(y,i)+Good_terms(y,i))/(-1 .+y[i,2]);
Scalarorigin(y,i)=(Divergent_terms(y,i)+(y[i,2])*(Der(y,i,3)))/(-1 .+y[i,2]);

Evans1(y,i,k)=2*((y[i+1,k])*R[i+1]-(y[i-1,k])*R[i-1])/(R[i+1]^2-R[i-1]^2);
Evans2(y,i,k)=3*((y[i+1,k])*R[i+1]^2-(y[i-1,k])*R[i-1]^2)/(R[i+1]^3-R[i-1]^3);

Der(y,i,k)= (y[i+1,k]-y[i-1,k])/(R[i+1]-R[i-1]);


#DD(y,i)=(1.0 .+y[i,2]-exp(y[i,1]))/R[i];
#CrRHS(y,i)=DD(y,i)+Der(y,i,2);

CrRHS(y,i)=4*pi*(-1.0 .+y[i,2])*R[i]*y[i,3]^2-(1.0 .+y[i,2])*(4.0*pi*R[i]*y[i,4]^2-Der(y,i,1));

DeltaRHS(y,i)=Der(y,i,1)-4*pi*R[i]*(y[i,3]-y[i,4])^2;


#Defining the bulk of the evolution

function bulkcr(y,i)
    #y=ghost(y);
    dy=zeros(length(y[1,:]));
    dy[1]=DeltaRHS(y,i)-dissipation4(y,i)[1];
    dy[2]=CrRHS(y,i)-dissipation4(y,i)[2];
    dy[3]=ScalarRHSCr(y,i)-dissipation4(y,i)[3];
    dy[4]=Der(y,i,3)-dissipation4(y,i)[4];
    dy[5]=y[i,3]-dissipation4(y,i)[5];
    return dy
end

# Defining the function in the RHS of the evolution equation system

function CRRHS(y,T)
    #dx=R[i+1]-R[i]
    L=length(R)
    dy=zeros((L,length(y[1,:])));
	for i in 3:(L-3)
        	dy[i,:]=bulkcr(y,i);
        end
    dy[L-2,:]=boundary(y,L-2);
    #dy[:,1:2]=zeros(L,2); #just scalar field in curved background, Cowling approximation
    return dy
end

upDer(y,i,k)=-(5*y[i,k]-11*y[i+1,k]+10*y[i+2,k]-5*y[i+3,k]+y[i+4,k])/(2*dx);

#switch=[1.0 0.0 1.0 1.0 1.0]';

function CRRHSexcision(y,T)
    dx=R[2]-R[1]
    L=length(R)
    dy=zeros((L,length(y[1,:])));
    for i in ex:L-2
        if i==L-2
        dy[i,:]=boundary(y,i)-dissipation4(y,i);
        else
        dy[i,:]=bulkcr(y,i)-dissipation4(y,i);
	if false #upwinding instead of centered finite differences at the excision boundary, Keep this switch off
	if i==ex 
			dy[i,1]=upDer(y,i,1)-4.0*pi*R[i]*(y[i,3]-y[i,4])^2-dissipation4(y,i)[1];
			dy[i,2]=4.0*pi*(-1.0 .+y[i,2])*R[i]*y[i,3]^2-(1.0 .+y[i,2])*(4.0*pi*R[i]*y[i,4]^2-upDer(y,i,1))-dissipation4(y,i)[2];
				#dy[i,3]=upDer(y,i,3)*y[i,2]+y[i,3]*upDer(y,i,2)-upDer(y,i,4)*(1.0 .+y[i,2])-y[i,4]*upDer(y,i,2)+2*(y[i,3]*(y[i,2])-y[i,4]*(1.0 .+y[i,2]))/R[i]-updissipationN(y,i)[3];
				#dy[i,4]=upDer(y,i,3)-updissipationN(y,i)[4]
		elseif i==ex+1
			#dy[i,1]=upDer(y,i,1)-4.0*pi*R[i]*(y[i,3]-y[i,4])^2-updissipationNN(y,i)[1];
                	#dy[i,2]=4.0*pi*(-1.0 .+y[i,2])*R[i]*y[i,3]^2-(1.0 .+y[i,2])*(4.0*pi*R[i]*y[i,4]^2-upDer(y,i,1))-updissipationNN(y,i)[2];
                        	#dy[i,3]=upDer(y,i,3)*y[i,2]+y[i,3]*upDer(y,i,2)-upDer(y,i,4)*(1.0 .+y[i,2])-y[i,4]*upDer(y,i,2)+2*(y[i,3]*(y[i,2])-y[i,4]*(1.0 .+y[i,2]))/R[i]-updissipationNN(y,i)[3];
                        	#dy[i,4]=upDer(y,i,3)-updissipationNN(y,i)[4]
			dy[i,1:2]=bulkcr(y,i)[1:2]-dissipation4(y,i)[1:2];
		end
    	end
	end
    end
    #dy[:,1:2]=zeros(L,2); #just scalar field in curved background, Cowling approximation
   	return dy
end

using Interpolations

function extrapolation(y,i)
	indi=i:(i+3);
	indi2=(i-1):(i+2)
	for j in 1:5
		y[i-1,j]=CubicSplineInterpolation(indi, y[indi,j], extrapolation_bc = Line())(i-1);
	end
	for j in 1:5
		y[i-2,j]=CubicSplineInterpolation(indi2, y[indi2,j], extrapolation_bc = Line())(i-2);
	end
return y
end

function outextrapolation(y,i)
	 indi=(i-3):i;
         indi2=(i-2):(i+1);
        for j in 1:5
                y[i+1,j]=CubicSplineInterpolation(indi, y[indi,j], extrapolation_bc = Line())(i+1);
        end
	for j in 1:5
                y[i+2,j]=CubicSplineInterpolation(indi, y[indi,j], extrapolation_bc = Line())(i+2);
        end
	return y
end




