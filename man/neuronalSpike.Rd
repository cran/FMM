\docType{data}
\name{neuronalSpike}
\alias{neuronalSpike}
\title{Neuronal Spike Data simulated with Hodgkin-Huxley model}
\format{
	A numeric vector.
}
\source{
	'NeuroDynex' Documentation, <https://lcn-neurodynex-exercises.readthedocs.io/en/latest/#>
}
\usage{
data(neuronalSpike)
}
\description{
	Voltage data in mV simulated with Hodgkin Huxley model (parameters: C=1,
	gNa=260, gK=30, gL=0.31, vK=-12, vNa=115, vL=10.6, bar(alphaN)=1.15,
	bar(betaN)=0.85, bar(alphaM)=0.9, bar(betaM)=1.3, bar(alphaH)=1,
	bar(betaH)=1) and applied current of 12 microA 1 millisecond. The
	simulation has been done with a modified 'NeuroDynex Python' module.
}
\examples{
data(neuronalSpike)
str(neuronalSpike)
}
\references{
	Wulfram Gerstner, Werner M. Kistler, Richard Naud,
	and Liam Paninski (2014). Neuronal Dynamics: From Single Neurons to
	Networks and Models of Cognition.
	([Online Book](https://neuronaldynamics.epfl.ch/online/))
}
\keyword{datasets}
