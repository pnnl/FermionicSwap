namespace FermionicSwap
{
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Chemistry.JordanWigner;
    open Microsoft.Quantum.Simulation;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Arrays;

    operation FermionicSwap(a : Qubit, b : Qubit) : Unit is Adj + Ctl {
        SWAP(a,b);
        CZ(a,b);
    }

    operation FermionicSwapLayer( swaps : (Int,Int)[], qubits : Qubit[]) : Unit is Adj + Ctl {
        for (a,b) in swaps {
            FermionicSwap(qubits[a],qubits[b]);
        }
    }

    operation FermionicSwapEvolveUnderGenerator(
        generator : EvolutionGenerator,
        trotterStepSize : Double,
        time : Double,
        register : Qubit[]
    ) : Unit is Adj + Ctl {
        let evolveFor = (FermionicSwapSimulationAlgorithm(trotterStepSize))!;
        evolveFor(time, generator, register);
    }

    operation FermionicSwapTrotterStep(
        swap_network : (Int,Int)[][],
        local_evolutions : JWOptimizedHTerms[][],
        time : Double,
        qubits : Qubit[]) : Unit
    {
        let nTerms = Length(qubits);
        for i in 0 .. Length(swap_network) {
            for ops in local_evolutions[i] {
                mutable empty = true;
                let (opa,opb,opc,opd) = ops!;
                if Length(opa) > 0 or Length(opb) > 0 or Length(opc) > 0 or Length(opd) > 0 {
                    set empty = false;
                }
                if (not empty) {
                    let generatorSystem = JordanWignerGeneratorSystem(ops);
                    let evolutionGenerator = EvolutionGenerator(JordanWignerFermionEvolutionSet(), generatorSystem);
                    TrotterStep(evolutionGenerator, 1, time)(qubits);
                }
            }
            if i < Length(swap_network) {
                FermionicSwapLayer(swap_network[i], qubits);
            }

            // not sure why this is here. leaving for the moment
            // (index:int, stepsize: Double, Qubits[]) -> Unit
            //let op_fun = (norb, time, op) -> 
            //let op = (nTerms, _ApplyHubbardTerm(nSites, tCoefficient, uCoefficient, _, _, _));
            // return DecomposedIntoTimeStepsCA(op, trotterOrder)(trotterStepSize, _);
        }
    }

    operation FermionicSwapEvolutionImpl(
        swap_network : (Int,Int)[][],
        local_evolutions : JWOptimizedHTerms[][],
        generatorIndex : GeneratorIndex,
        time : Double,
        qubits : Qubit[]
    ) : Unit is Adj + Ctl{
        body (...) {
            let ((indices, _), _) = generatorIndex!;
            let index = indices[0];
            let gi = (index-1) / 2;
            if index % 2 != 0 {
                for ops in local_evolutions[gi] {
                    mutable empty = true;
                    let (opa,opb,opc,opd) = ops!;
                    if (Length(opa) > 0 or Length(opb) > 0 or Length(opc) > 0 or Length(opd) > 0) {
                        let generatorSystem = JordanWignerGeneratorSystem(ops);
                        let evolutionGenerator = EvolutionGenerator(JordanWignerFermionEvolutionSet(), generatorSystem);
                        TrotterStep(evolutionGenerator, 1, time)(qubits);
                    }
                }
            } else {
                FermionicSwapLayer(swap_network[gi], qubits);
            }
        }
        // adjoint invert;
        // controlled distribute;
        // controlled adjoint distribute;
    }

    function FermionicSwapEvolutionFunction(
        swap_network : (Int,Int)[][],
        local_evolutions : JWOptimizedHTerms[][],
        generatorIndex : GeneratorIndex
    ) : EvolutionUnitary {
        return EvolutionUnitary(FermionicSwapEvolutionImpl(swap_network, local_evolutions, generatorIndex, _, _));
    }

    function FermionicSwapEvolutionSet(
        swap_network : (Int,Int)[][],
        local_evolutions : JWOptimizedHTerms[][]
    ) : EvolutionSet {
        return EvolutionSet(FermionicSwapEvolutionFunction(swap_network, local_evolutions, _));
    }

    function FermionicSwapGeneratorSystem(
        size : Int
    ) : GeneratorSystem {
        return GeneratorSystem(size, s -> GeneratorIndex(([s], []), []));
    }

    operation FermionicSwapSimulationAlgorithmImpl(
        step_size : Double,
        time : Double,
        evo_gen : EvolutionGenerator,
        qubits : Qubit[]
    ) : Unit is Adj + Ctl {
        let (evo_set, gen_sys) = evo_gen!;
        let (num_terms, term_dict) = gen_sys!;
        let time_steps = Ceiling(time / step_size);
        for i in 1..time_steps {
            let this_time = (i<time_steps ? step_size | time - (step_size * IntAsDouble(time_steps-1))); 
            for s in i%2==1 ? (1..num_terms) | (num_terms..-1..1) {
                evo_set!(term_dict(s))!(this_time,qubits);
            }
        }
    }
    function FermionicSwapSimulationAlgorithm(
        step_size : Double
    ) : SimulationAlgorithm {
        return SimulationAlgorithm(FermionicSwapSimulationAlgorithmImpl(
            step_size, _,_,_));
    }

    function FermionicSwapEvolutionGenerator(
        swap_network : (Int,Int)[][],
        local_evolutions : JWOptimizedHTerms[][]
    ) : EvolutionGenerator {
        return EvolutionGenerator(
            FermionicSwapEvolutionSet(swap_network, local_evolutions),
            FermionicSwapGeneratorSystem(Length(swap_network) + Length(local_evolutions))
        );
    }
}