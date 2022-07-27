function simulation = generateSimulation(Nmolecules,Ntraps)
    for i = 1:Nmolecules; molecules{i} = moleculeClass(); end
    for i = 1:Ntraps; traps{i} = trapClass(); end
    simulation = simulationClass();
    simulation.addObjects( molecules );
    simulation.addObjects( traps );
end