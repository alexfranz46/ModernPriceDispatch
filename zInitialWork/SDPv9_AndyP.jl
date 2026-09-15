# using CSV
# using DataFrames

#=
This file is a copy of what Andy used to write section 5 of the ADR paper.
Only modification are stylistic and output (print instead of saving to files.)
=#


# Stochastic battery operation to meet duck curve using dynmaic programming
numperiods = 24;
numthermal = 1;
numbatteries = 1;
incomingbattery = 0;
incominggeneration = 35;
batterycapacity = 30;
capacity = 70;
rampup = 10;
voll = 1000;
eta = zeros(2);
eta[1] = 1.0;
r = 15;
s = 15;

# Demand scenarios
# DemandDF = CSV.read("demand.csv", DataFrame)
# demand = zeros(numperiods)
demand = [40, 41, 42, 43, 35, 40, 40, 25, 10, 8, 6, 5, 5, 6, 8, 10, 20, 30, 55, 72, 75, 70, 64, 60]
O = Vector{Vector{Float64}}(undef, numperiods)
P = Vector{Vector{Float64}}(undef, numperiods)
for t in 1:numperiods
    # demand[t] = DemandDF[t,2]
    O[t] = [-4.0, -2.0, 0, 2.0, 4.0]
    P[t] = [0.2, 0.2, 0.2, 0.2, 0.2]
end

# Thermal generation 

numtranches = 10
tranchecap = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 10.0, 10.0, 10.0, 10.0]
tranchecost = [10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 90.0, 110.0, 150.0, 200.0]

function gencost(g::Float64)
    mycost = 0.0
    i = 0
    notdone = true
    cumulative = 0.0
    while ((cumulative < g) && notdone)
        i = i + 1
        if g < cumulative + tranchecap[i]
            mycost = mycost + (g - cumulative) * tranchecost[i]
            notdone = false
        else
            mycost = mycost + tranchecap[i] * tranchecost[i]
            cumulative = cumulative + tranchecap[i]
        end
    end
    return mycost
end

numscenarios = 5;

CostToGo = zeros(25, 71, 31);            # CTG at start of period t
mystagecost = zeros(24, 71, 31, 5);
mygenerationcost = zeros(24, 71, 31, 5);
myoptgen = zeros(24, 71, 31, 5);
myoptcharge = zeros(24, 71, 31, 5);
myslack = zeros(24, 71, 31, 5);
mystorage = zeros(24, 71, 31, 5);        # storage at end of period t

for time = 1:numperiods

    global t = numperiods + 1 - time
    @info("t = ", t)

    for i = 1:71    # generation level+1
        for j = 1:31    # battery level +1
            AveCost = 0

            for omega = 1:numscenarios
                # Find best decision for outcome omega	
                MinCost = 1000000
                stagecost = 0.0
                generationcost = 0.0
                optgen = 0
                optcharge = 0.0
                storage = 0
                theslack = 0.0
                mydemand = demand[t] + O[t][omega]


                for newgeneration = 0:(i-1+rampup)
                    if (newgeneration <= capacity)
                        for charge = -s:r
                            #		net = Int(eta[1]*charge+discharge)
                            net = Int(charge)
                            # net is charge going in to battery
                            slack = 0
                            if (j - 1 + net >= 0) && (j - 1 + net <= batterycapacity) # charge is feasible
                                if (newgeneration - net < mydemand)
                                    slack = mydemand - (newgeneration - net)
                                end
                                if (newgeneration - net <= mydemand)  # TODO: THIS IS WHERE THE MODELS DIFFER! ANDY DISALLOWS SPILL/CURTAILMENT
                                    CTG = CostToGo[t+1, newgeneration+1, j+net]
                                    if slack * voll + gencost(Float64(newgeneration)) + CTG < MinCost
                                        MinCost = slack * voll + gencost(Float64(newgeneration)) + CTG
                                        stagecost = slack * voll + gencost(Float64(newgeneration))
                                        generationcost = gencost(Float64(newgeneration))
                                        optgen = Float64(newgeneration)
                                        optcharge = Float64(net)
                                        storage = Float64(j - 1 + net)
                                        theslack = slack
                                    end
                                end
                            end
                        end   # charge
                    end  # if (newgeneration <= capacity)
                end
                mystagecost[t, i, j, omega] = stagecost
                mygenerationcost[t, i, j, omega] = generationcost
                myoptgen[t, i, j, omega] = optgen
                myoptcharge[t, i, j, omega] = optcharge
                myslack[t, i, j, omega] = theslack
                mystorage[t, i, j, omega] = storage  # storage at end of period t
                AveCost = AveCost + MinCost * P[t][omega]
            end     # of scenarios
            CostToGo[t, i, j] = AveCost  # Expected CTG at start of period t
        end
    end

end

println(CostToGo[1,35,0])  # objetive.

if false

	for t = 1:numperiods
		open("CTG" * string(t) * ".csv", "w") do f
			for i in 1:71
				for j in 1:31
					print(f, CostToGo[t, i, j])
					print(f, ", ")
				end
				println(f, "")
			end
		end
	end

	for t = 1:numperiods
		open("optgen" * string(t) * ".csv", "w") do f
			for omega = 1:1
				println(f, "Omega = ", omega)
				for i in 1:71
					for j in 1:31
						print(f, myoptgen[t, i, j, omega])
						print(f, ", ")
					end
					println(f, "")
				end
			end
		end
	end

	for t = 1:numperiods
		open("loadshed" * string(t) * ".csv", "w") do f
			for omega = 1:1
				println(f, "Omega = ", omega)
				for i in 1:71
					for j in 1:31
						print(f, myslack[t, i, j, omega])
						print(f, ", ")
					end
					println(f, "")
				end
			end
		end
	end

	for t = 1:numperiods
		open("mystorage" * string(t) * ".csv", "w") do f
			for omega = 1:1
				println(f, "Omega = ", omega)
				for i in 1:71
					for j in 1:31
						print(f, mystorage[t, i, j, omega])
						print(f, ", ")
					end
					println(f, "")
				end
			end
		end
	end
end

if false  # skip this stuff
    sequence = zeros(24)
    optchargein = zeros(24)
    shed = zeros(24)
    state = zeros(24)  # storage at end of period t
    state[1] = mystorage[1, Int(incominggeneration)+1, Int(incomingbattery)+1, 1]
    #incomingbattery;
    sequence[1] = myoptgen[1, Int(incominggeneration)+1, Int(incomingbattery)+1, 1]
    shed[1] = myslack[1, Int(incominggeneration)+1, Int(incomingbattery)+1, 1]
    for t = 2:24
        state[t] = mystorage[t, Int(sequence[t-1])+1, Int(state[t-1])+1, 1]
        sequence[t] = myoptgen[t, Int(sequence[t-1])+1, Int(state[t-1])+1, 1]
        optchargein[t] = myoptcharge[t, Int(sequence[t-1])+1, Int(state[t-1])+1, 1]
        shed[t] = myslack[t, Int(sequence[t-1])+1, Int(state[t-1])+1, 1]
    end
    #state[25] = mystorage[25,Int(sequence[24]+1),Int(state[24]+1),1]

    open("prod.csv", "w") do f
        for t in 1:24
            println(f, sequence[t])
        end
    end

    open("chargein.csv", "w") do f
        for t in 1:24
            println(f, optchargein[t])
        end
    end

    open("shed.csv", "w") do f
        for t in 1:24
            println(f, shed[t])
        end
    end

    open("storeAtend.csv", "w") do f
        for t in 1:24
            println(f, state[t])
        end
    end

    global testcost = 0.0
    for t = 1:24
        global testcost = testcost + gencost(Float64(sequence[t]))
    end

end
