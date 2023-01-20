namespace FermionicSwap
{
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Chemistry.JordanWigner;
    open Microsoft.Quantum.Simulation;

    operation FermionicSwap(a : Qubit, b : Qubit) : Unit is Adj + Ctl {
        SWAP(a,b);
        CZ(a,b);
    }

    operation FermionicSwapLayer( swaps : (Int,Int)[], qubits : Qubit[]) : Unit is Adj + Ctl {
        for (a,b) in swaps {
            FermionicSwap(qubits[a],qubits[b]);
        }
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
}