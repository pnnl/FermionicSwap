namespace FermionicSwap.Tests {
    open Microsoft.Quantum.Simulation;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Chemistry;
    open Microsoft.Quantum.Chemistry.JordanWigner;
    open Microsoft.Quantum.Arrays;
    open FermionicSwap;

    // We open the diagnostics namespace under an alias to help avoid
    // conflicting with deprecation stubs in Microsoft.Quantum.Canon.
    open Microsoft.Quantum.Diagnostics as Diag;

    // This test as currently written will only work for Hamiltonians with a
    // single summand due to Trotter summand reordering issues.
    operation SwapNetworkTestOp(swap_network : (Int,Int)[][],
                                qsharp_network_data : JWOptimizedHTerms[][],
                                qsharp_hamiltonian : JWOptimizedHTerms,
                                num_qubits : Int
                                ) : Unit {
        let time = 1.0;
        Diag.AssertOperationsEqualInPlace(num_qubits,
            FixedOrderFermionicSwapTrotterStep(swap_network, qsharp_network_data, time, _),
            JordanWignerApplyTrotterStep(qsharp_hamiltonian, time, _ )
        );
    }

    operation FixedOrderFermionicSwapTrotterStep(swap_network : (Int,Int)[][],
                                qsharp_network_data : JWOptimizedHTerms[][],
                                time : Double, qubits : Qubit[]) : Unit {
        FermionicSwapTrotterStep(swap_network, qsharp_network_data, time, qubits);
        let empty = new JWOptimizedHTerms[][Length(swap_network)+1];
        FermionicSwapTrotterStep(Reversed(swap_network), empty, 0.0, qubits);
    }

    operation JordanWignerApplyTrotterStep (data : JWOptimizedHTerms, trotterStepSize : Double, qubits :
Qubit[])
    : Unit is Adj + Ctl {
        let generatorSystem = JordanWignerGeneratorSystem(data);
        let evolutionGenerator = EvolutionGenerator(JordanWignerFermionEvolutionSet(), generatorSystem);
        let trotterOrder = 1;
        TrotterStep(evolutionGenerator, trotterOrder, trotterStepSize)(qubits);
    }
}