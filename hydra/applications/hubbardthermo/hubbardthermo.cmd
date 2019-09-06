Opt std::string outfile "" o outfile outfile "name of outfile"
Opt std::string latticefile "" l latticefile latticefile "name of latticefile"
Opt std::string couplingfile "" c couplingfile couplingfile "name of couplingfile (default: latticefile)"
Opt std::string ensemble "canonical" e ensemble ensemble "ensemble to use, either canonical or grandcanonical (default: canonical)" 
Opt std::string temperaturefile "" t temperaturefile temperaturefile "name of temperaturefile"
Opt int np -1 p np np "number of particles for canonical ensemble (optional, default half filling)"
Opt double mu 0. m mu mu "chemical potential for grandcanonical ensemble (optional, default 0.)"
