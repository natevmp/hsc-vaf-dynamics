using JLD2

@load "data/SinglePatientPipeline/Nf10000/singlePatientFullSim_Nf10000_sim01.jld2"


##
using Plots
fscale = 0.8
pyplot(legendfontsize=10, guidefontsize=12, tickfontsize=10,  size=(fscale*600,fscale*400))

##
freqsS_f = (0:paramsTrue["sample size"]) / paramsTrue["sample size"]

fig1 = scatter(freqsS_f[2:end], (nVarSimS_f/nVarSimS_f[2])[2:end],
    yscale=:log10,
    ylims=(1.2*10^-5,1.3*10^0),
    markerstrokealpha=0,
    markersize=5,
    label="occupied states"
)

freqs0_f = freqsS_f[nVarSimS_f.==0]
scatter!(freqs0_f, 1.5E-5*ones(length(freqs0_f)),
    markersize=5,
    markerstrokealpha=0,
    label="unoccupied states"
)

xlabel!("variant frequency")
ylabel!("number of variants")

display(fig1)

savefig(fig1, "singleSimVafSpectrumPlot.pdf")