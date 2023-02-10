// Copyright Battelle Memorial Institute 2022. All rights reserved.

using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.LadderOperators;
using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;
using Microsoft.Quantum.Chemistry.JordanWigner;
using static FermionicSwap.FSTools;
using Microsoft.Quantum.Chemistry.QSharpFormat;
using static FermionicSwap.SwapNetwork;
namespace FermionicSwap.Tests
{

    using System.Linq;
    using System.Collections.Immutable;
    using SwapLayer = List<(int,int)>;
    using OperatorLayer = List<(HermitianFermionTerm, DoubleCoeff)>;
    using OperatorNetwork = List<List<(HermitianFermionTerm, DoubleCoeff)>>;

    public class TestFermionicSwap
{

    [Theory]
    [MemberData(nameof(Data))]
    public void TestEvenOddSwap(int[] start_order, int[] end_order, SwapNetwork swap_network)
    {
        var result = FSTools.evenOddSwap(start_order, end_order);
        Assert.True(result.Count == swap_network.Count, $"Need swap layers {layers_string(swap_network)}, but got {layers_string(result)}.");
        // for each swap layer, check that the list of swaps matches the test list
        foreach (var (first,second) in result.Zip(swap_network, (f,s)=> (f,s))) {
            Assert.True(first.SequenceEqual(second), $"Swap network {layers_string(result)} should be {layers_string(swap_network)}.");
        }
    }
 
    public static IEnumerable<object[]> Data => new List<object[]>
    {
        // empty layer
        new object[] { new int[] {}, new int[] {}, new SwapNetwork {} },
        // trivial swapping
        new object[] { new int[] {0,1}, new int[] {0,1}, new SwapNetwork {} },
        // nontrivial
        new object[] { new int[] {1,0}, new int[] {0,1}, new SwapNetwork {new SwapLayer {(0,1)}}},
        // site numbering does not start from zero
        new object[] { new int[] {1,2}, new int[] {2,1}, new SwapNetwork {new SwapLayer {(0,1)}}},
        // three items, trivial swapping
        new object[] { new int[] {0,1,2}, new int[] {0,1,2}, new SwapNetwork {}},
        // four items, trivial swapping
        new object[] { new int[] {0,1,2,3}, new int[] {0,1,2,3}, new SwapNetwork {}},
        // five items, trivial swapping
        new object[] { new int[] {0,1,2,3,4}, new int[] {0,1,2,3,4}, new SwapNetwork {}},
        // three items, nontrivial swapping
        new object[] { new int[] {0,1,2}, new int[] {2,1,0}, new SwapNetwork {
            new SwapLayer {(0,1)},
            new SwapLayer {(1,2)},
            new SwapLayer {(0,1)}
            }
        },
        // three items, no (even) swaps in initial swap layer
        new object[] { new int[] {0,1,2}, new int[] {0,2,1}, new SwapNetwork {new SwapLayer {(1,2)}}},
        // 7 items, a single item moves in each layer
        new object[] { new int[] {0,1,2,3,4,5,6}, new int[] {6,0,1,2,3,4,5}, new SwapNetwork {
            new SwapLayer {(5,6)},
            new SwapLayer {(4,5)},
            new SwapLayer {(3,4)}, new SwapLayer {(2,3)},
            new SwapLayer {(1,2)}, new SwapLayer {(0,1)},
            }},
        // odd larger number of items, full reverse
        new object[] { new int[] {0,1,2,3,4,5,6}, new int[] {6,5,4,3,2,1,0}, new SwapNetwork {
            new SwapLayer {(0,1),(2,3),(4,5)},
            new SwapLayer {(1,2),(3,4),(5,6)},
            new SwapLayer {(0,1),(2,3),(4,5)},
            new SwapLayer {(1,2),(3,4),(5,6)},
            new SwapLayer {(0,1),(2,3),(4,5)},
            new SwapLayer {(1,2),(3,4),(5,6)},
            new SwapLayer {(0,1),(2,3),(4,5)},
        }},
        // even larger number of items, full reverse
        new object[] { new int[] {0,1,2,3,4,5}, new int[] {5,4,3,2,1,0}, new SwapNetwork {
            new SwapLayer {(0,1),(2,3),(4,5)},
            new SwapLayer {(1,2),(3,4)},
            new SwapLayer {(0,1),(2,3),(4,5)},
            new SwapLayer {(1,2),(3,4)},
            new SwapLayer {(0,1),(2,3),(4,5)},
            new SwapLayer {(1,2),(3,4)},
        }}
    };

    [Theory]
    [MemberData(nameof(OneBodyDenseNetworkData))]
    public void TestOneBodyDenseNetwork(int num_sites, SwapNetwork swap_network)
    {
        var result = FSTools.OneBodyDenseNetwork(num_sites);
        Assert.True(result.Count == swap_network.Count, $"Need swap layers {layers_string(swap_network)}, but got {layers_string(result)}.");
        // for each swap layer, check that the list of swaps matches the test list
        foreach (var (first,second) in result.Zip(swap_network, (f,s)=> (f,s))) {
            Assert.True(first.SequenceEqual(second), $"Swap network {layers_string(result)} should be {layers_string(swap_network)}.");
        }
    }
    public static IEnumerable<object[]> OneBodyDenseNetworkData => new List<object[]>
    {
        // Small one body swap networks of various sizes. Small, larger odd and larger even.
        new object[] {0, new SwapNetwork {}},
        new object[] {1, new SwapNetwork {}},
        new object[] {2, new SwapNetwork {new SwapLayer {(0,1)}}},
        new object[] {3, new SwapNetwork {new SwapLayer {(0,1)}, new SwapLayer {(1,2)}, new SwapLayer {(0,1)}}},
        new object[] {6, new SwapNetwork {
            new SwapLayer {(0,1),(2,3),(4,5)},
            new SwapLayer {(1,2),(3,4)},
            new SwapLayer {(0,1),(2,3),(4,5)},
            new SwapLayer {(1,2),(3,4)},
            new SwapLayer {(0,1),(2,3),(4,5)},
            new SwapLayer {(1,2),(3,4)},
        }},
        new object[] {7, new SwapNetwork {
            new SwapLayer {(0,1),(2,3),(4,5)},
            new SwapLayer {(1,2),(3,4),(5,6)},
            new SwapLayer {(0,1),(2,3),(4,5)},
            new SwapLayer {(1,2),(3,4),(5,6)},
            new SwapLayer {(0,1),(2,3),(4,5)},
            new SwapLayer {(1,2),(3,4),(5,6)},
            new SwapLayer {(0,1),(2,3),(4,5)},
        }},
    };

    [Theory]
    [MemberData(nameof(ReorderedFermionTermData))]
    public void TestReorderedFermionTerm(HermitianFermionTerm term,
                                        Dictionary<int,int> desired_order,
                                        List<int> correct_term,
                                        int coefficient) {
        var new_term = ReorderedFermionTerm(term, desired_order);
        var result = new_term.Sequence.Select(o => o.Index);
        var result_string = String.Join(", ", result);
        Assert.True(result.SequenceEqual(correct_term), $"Incorrect order {result_string}.");
        Assert.Equal(coefficient, new_term.Coefficient);
    }

    // Note: the reordered terms are returned in QDK's canonical ladder operator order.
    public static IEnumerable<object[]> ReorderedFermionTermData => new List<object[]>
    {
        // Leave a correctly ordered object alone.
        new object[] {
            new HermitianFermionTerm(new int[] {0,1}),
            PositionDictionary(new int[] {0,1}),
            new List<int> {0,1},
            1
            },
        new object[] {
            new HermitianFermionTerm(new int[] {0,1,3,2}),
            PositionDictionary(new int[] {0,1,2,3}),
            new List<int> {0,1,3,2},
            1
            },
        new object[] {
            new HermitianFermionTerm(new int[] {0,5,6,4}),
            PositionDictionary(new int[] {0,1,2,3,4,5,6}),
            new List<int> {0,5,6,4},
            1
            },

        // Permute to new positions correctly in the absence of canonical reordering
        new object[] {
            new HermitianFermionTerm(new int[] {0,1}),
            PositionDictionary(new int[] {1,0}),
            // Hermitian reordering occurs here
            new List<int> {0,1},
            1
            },
        new object[] {
            new HermitianFermionTerm(new int[] {1,3,2,0}),
            PositionDictionary(new int[] {0,1,3,2}),
            new List<int> {0,3,2,1},
            1
            },
        new object[] {
            new HermitianFermionTerm(new int[] {2,4,3,1}),
            PositionDictionary(new int[] {0,1,2,4,3}),
            new List<int> {1,4,3,2},
            1
            },
        // Permute to new positions with pre/post canonical reordering
        // differing by even/even permutations from given ordering
        new object[] {
            new HermitianFermionTerm(new int[] {0,1,2,3,4,5,6,7}),
            PositionDictionary(new int[] {1,0,3,2,7,6,5,4}),
            new List<int> {0,1,2,3,7,6,5,4},
            1
            },
        // Same, but even/odd        
        new object[] {
            new HermitianFermionTerm(new int[] {0,1,2,3,4,5,6,7}),
            PositionDictionary(new int[] {0,1,3,2,7,6,5,4}),
            new List<int> {0,1,2,3,7,6,5,4},
            -1
            },       
        // Same, but odd/even        
        new object[] {
            new HermitianFermionTerm(new int[] {0,1,2,3,4,5,7,6}),
            PositionDictionary(new int[] {1,0,3,2,7,6,5,4}),
            new List<int> {0,1,2,3,7,6,5,4},
            -1
            },
        // Same, but odd/odd        
        new object[] {
            new HermitianFermionTerm(new int[] {0,1,2,3,4,5,7,6}),
            PositionDictionary(new int[] {0,1,3,2,7,6,5,4}),
            new List<int> {0,1,2,3,7,6,5,4},
            1
            },
    };

    [Theory]
    [MemberData(nameof(ProcessNetworkLayerData))]
    public void TestProcessNetworkLayer(
            TermsDictionary term_dict,
            int[] order,
            OperatorLayer correct)
    {
        var result = ProcessNetworkLayer(term_dict, order);
        Assert.True(result.Count() == correct.Count(), $"Result has {result.Count()} elements but correct result has {correct.Count()}");
        foreach (var (r,c) in result.Zip(correct)) {
            Assert.True(r == c, $"{r} does not equal {c}.");
        }
    }

    // Note: the reordered terms are returned in QDK's canonical ladder operator order.
    public static IEnumerable<object[]> ProcessNetworkLayerData => new List<object[]>
    {
        // Produce no operators if there are no terms.
        new object[] {
            new TermsDictionary(),
            new int[] {0,1,2,3},
            new OperatorLayer {}
        },
        // Produce an operator from a term
        new object[] {
            new TermsDictionary() {{ImmutableArray.Create(new int[] {0,1}), new OperatorLayer {
                (new HermitianFermionTerm(new int[] {0,1}),3.0)
            }}},
            new int[] {0,1},
            new OperatorLayer{(new HermitianFermionTerm(new int[] {0,1}),3.0)}
        },
        // Produce an operator from a misordered term
        new object[] {
            new TermsDictionary() {{ImmutableArray.Create(new int[] {0,1}), new OperatorLayer {
                (new HermitianFermionTerm(new int[] {0,1}),3.0)
            }}},
            new int[] {0,1},
            new OperatorLayer{(new HermitianFermionTerm(new int[] {0,1}),3.0)}
        },
        // Apply the greedy algorithm to produce multiple operators
        new object[] {
            new TermsDictionary() {
                {ImmutableArray.Create(new int[] {0,1}), new OperatorLayer {
                    (new HermitianFermionTerm(new int[] {0,1}),1.0)
                }},
                {ImmutableArray.Create(new int[] {3,4}), new OperatorLayer {
                    (new HermitianFermionTerm(new int[] {3,4}),1.0)
                }},

            },
            new int[] {0,1,2,3,4},
            new OperatorLayer{
                (new HermitianFermionTerm(new int[] {0,1}),1.0),
                (new HermitianFermionTerm(new int[] {3,4}),1.0)
            }
        },
        // Process multiple operators with the same indices
        // double check time application
        new object[] {
            new TermsDictionary() {
                {ImmutableArray.Create(new int[] {0,1}), new OperatorLayer {
                    (new HermitianFermionTerm(new int[] {0,1}),1.0),
                    (new HermitianFermionTerm(new int[] {0,1,1,0}),2.0)
                }},
                {ImmutableArray.Create(new int[] {3,4}), new OperatorLayer {
                    (new HermitianFermionTerm(new int[] {3,4}),3.0)
                }},

            },
            new int[] {0,1,2,3,4},
            new OperatorLayer{
                (new HermitianFermionTerm(new int[] {0,1}),1.0),
                (new HermitianFermionTerm(new int[] {3,4}),3.0),
                (new HermitianFermionTerm(new int[] {0,1,1,0}),2.0)
            }
        },
        // Handle operator overlap correctly
        new object[] {
            new TermsDictionary() {
                {ImmutableArray.Create(new int[] {0,1}), new OperatorLayer {
                    (new HermitianFermionTerm(new int[] {0,1}),1.0),
                    (new HermitianFermionTerm(new int[] {0,1,1,0}),1.0)
                }},
                {ImmutableArray.Create(new int[] {1,2,3,4}), new OperatorLayer {
                    (new HermitianFermionTerm(new int[] {2,4,3,1}),1.0)
                }},
                {ImmutableArray.Create(new int[] {3,4}), new OperatorLayer {
                    (new HermitianFermionTerm(new int[] {3,4}),1.0)
                }},

            },
            new int[] {0,1,2,3,4},
            new OperatorLayer{
                (new HermitianFermionTerm(new int[] {0,1}),1.0),
                (new HermitianFermionTerm(new int[] {3,4}),1.0),
                (new HermitianFermionTerm(new int[] {0,1,1,0}),1.0),
                (new HermitianFermionTerm(new int[] {2,4,3,1}),1.0)
            }
        },
        // Correctly reorder the terms
        new object[] {
            new TermsDictionary() {
                {ImmutableArray.Create(new int[] {0,1}), new OperatorLayer {
                    (new HermitianFermionTerm(new int[] {0,1}),1.0),
                    (new HermitianFermionTerm(new int[] {0,1,1,0}),1.0)
                }},
                {ImmutableArray.Create(new int[] {1,2,3,4}), new OperatorLayer {
                    (new HermitianFermionTerm(new int[] {2,4,3,1}),1.0)
                }},
                {ImmutableArray.Create(new int[] {3,4}), new OperatorLayer {
                    (new HermitianFermionTerm(new int[] {3,4}),1.0)
                }},

            },
            new int[] {4,3,2,1,0},
            new OperatorLayer{
                (new HermitianFermionTerm(new int[] {1,0}),1.0), // 3,4
                (new HermitianFermionTerm(new int[] {4,3}),1.0), // 0,1
                (new HermitianFermionTerm(new int[] {0,2,3,1}),1.0), // 2,4,3,1
                (new HermitianFermionTerm(new int[] {4,3,3,4}),1.0), // 0,1,1,0
            }
        }

    };

    [Theory]
    [MemberData(nameof(TrotterStepDataData))]
    public void TestTrotterStepData(
            FermionHamiltonian H,
            SwapNetwork swap_network,
            int[] start_order,
            OperatorNetwork correct_network,
            int[] correct_order
            )
    {
        var (operator_network, end_order) = TrotterStepData(H, swap_network, start_order);
        Assert.True(end_order.SequenceEqual(correct_order),
            $"Resulting order {String.Join(", ", end_order)} differs from correct order {String.Join(", ", correct_order.Select(o=>o.ToString()))}.");
        Assert.True(operator_network.Count() == swap_network.Count() + 1,
            $"Resulting operator network has {operator_network.Count()} layers instead of {swap_network.Count()}.");
        foreach (var (r,c) in operator_network.Zip(correct_network)) {
            Assert.True(r.SequenceEqual(c), $"Resulting layer {String.Join(", ", r)} differs from correct layer {String.Join(", ", c)}.");
        }
    }

    // Note: the reordered terms are returned in QDK's canonical ladder operator order.
    public static IEnumerable<object[]> TrotterStepDataData() {
        var result = new List<object[]> {};

        // An empty Hamiltonian and swap network produce an empty operator network.
        var H = new FermionHamiltonian {};
        var swap_network = new SwapNetwork {};
        var start_order = new int[]{};
        var correct_network = new OperatorNetwork {};
        int[] correct_order = start_order.ToArray();
        result.Add(new object[] {H, swap_network, start_order, correct_network, correct_order});

        // An empty Hamiltonian and any swap network produce an empty operator network.
        H = new FermionHamiltonian {};
        swap_network = OneBodyDenseNetwork(3);
        start_order = new int[] {0,1,2};
        correct_network = new OperatorNetwork {};
        correct_order = start_order.Reverse().ToArray();
        result.Add(new object[] {H, swap_network, start_order, correct_network, correct_order});

        // correct networks for some dense hopping term hamiltonians
        var correct_networks = new List<OperatorNetwork> {
            // 3 sites
            new OperatorNetwork {
                new OperatorLayer {
                    (new HermitianFermionTerm(new int[] {0,1}), 1.0),
                    (new HermitianFermionTerm(new int[] {1,2}), 1.0)
                },
                new OperatorLayer {
                    (new HermitianFermionTerm(new int[] {1,2}),1.0)
                },
                new OperatorLayer {}
            },
            // 4 sites
            new OperatorNetwork {
                new OperatorLayer {
                    (new HermitianFermionTerm(new int[] {0,1}), 1.0),
                    (new HermitianFermionTerm(new int[] {2,3}), 1.0),
                    (new HermitianFermionTerm(new int[] {1,2}), 1.0)
                },
                new OperatorLayer {
                    (new HermitianFermionTerm(new int[] {1,2}), 1.0),
                },
                new OperatorLayer {
                    (new HermitianFermionTerm(new int[] {0,1}), 1.0),
                    (new HermitianFermionTerm(new int[] {2,3}), 1.0),
                },
                new OperatorLayer {}
            },
            // 5 sites
            new OperatorNetwork {
                new OperatorLayer {
                    (new HermitianFermionTerm(new int[] {0,1}), 1.0),
                    (new HermitianFermionTerm(new int[] {2,3}), 1.0),
                    (new HermitianFermionTerm(new int[] {1,2}), 1.0),
                    (new HermitianFermionTerm(new int[] {3,4}), 1.0)
                },
                new OperatorLayer {
                    (new HermitianFermionTerm(new int[] {1,2}), 1.0),
                    (new HermitianFermionTerm(new int[] {3,4}), 1.0),
                },
                new OperatorLayer {
                    (new HermitianFermionTerm(new int[] {0,1}), 1.0),
                    (new HermitianFermionTerm(new int[] {2,3}), 1.0),
                },
                new OperatorLayer {
                    (new HermitianFermionTerm(new int[] {1,2}), 1.0),
                    (new HermitianFermionTerm(new int[] {3,4}), 1.0),
                },
                new OperatorLayer {}
            }
        };

        // Dense Hopping terms on 3,4,5 site orbitals
        int num_sites = 3;
        for (; num_sites < 6; num_sites++) {
            H = new FermionHamiltonian {};
            for (int i = 0; i < num_sites; i++) {
                for (int j = i+1; j < num_sites; j++) {
                    H.Add(new HermitianFermionTerm(new int[] {i, j}), 1.0);
                }
            }
            swap_network = OneBodyDenseNetwork(num_sites);
            start_order = Enumerable.Range(0,num_sites).ToArray();
            correct_order = start_order.Reverse().ToArray();
            result.Add(new object[] {H, swap_network, start_order, correct_networks[num_sites-3], correct_order});
        }

        // verify that weights transfer correctly
        H = new FermionHamiltonian {};
        num_sites = 5;
        for (int i = 0; i < num_sites; i++) {
            for (int j = i+1; j < num_sites; j++) {
                H.Add(new HermitianFermionTerm(new int[] {i, j}), (double)(10*i+j));
            }
        }
        swap_network = OneBodyDenseNetwork(num_sites);
        start_order = Enumerable.Range(0,num_sites).ToArray();
        correct_order = start_order.Reverse().ToArray();
        correct_network = new OperatorNetwork {
            new OperatorLayer {
                (new HermitianFermionTerm(new int[] {0,1}), 1.0),
                (new HermitianFermionTerm(new int[] {2,3}), 23.0),
                (new HermitianFermionTerm(new int[] {1,2}), 12.0),
                (new HermitianFermionTerm(new int[] {3,4}), 34.0)
            },
            // order 10324
            new OperatorLayer {
                (new HermitianFermionTerm(new int[] {1,2}), 3.0),
                (new HermitianFermionTerm(new int[] {3,4}), 24.0),
            },
            // order 13042
            new OperatorLayer {
                (new HermitianFermionTerm(new int[] {0,1}), 13.0),
                (new HermitianFermionTerm(new int[] {2,3}), 4.0),
            },
            //order 31402
            new OperatorLayer {
                (new HermitianFermionTerm(new int[] {1,2}), 14.0),
                (new HermitianFermionTerm(new int[] {3,4}), 2.0),
            },
            //order 34120
            new OperatorLayer {}
            //order 43210
        };

        result.Add(new object[] {H, swap_network, start_order, correct_network, correct_order});

        return result;
    }

    // Check the operation of the following functions:
    //     ToQSharpFormat,
    //     FermionicSwapTrotterStep (qsharp),
    //     FixedOrderFermionicSwapTrotterStep (qsharp),
    // by constructing a Hamiltonian for a Trotter step and checking for
    // equality with the corresponding JordanWigner trotter step. Since
    // the two methods do not agree on the order in which terms are
    // evaluated, which results in unequal Trotter steps in general, this
    // test uses Hamiltonians consisting of single PQ terms.
    [Theory]
    [MemberData(nameof(OneTermHamiltonianData))]
    public void TestOneTermHamiltonian(
            FermionHamiltonian H,
            SwapNetwork swap_network,
            int[] start_order
            )
    {
        var (op_network, end_order) = TrotterStepData(H, swap_network, start_order);
        //we use 32 bit ints until the point of injection into q#, which requires 64 bit ints.
        var qsharp_swap_network = swap_network.ToQSharpFormat();
        var qsharp_data = ToQSharpFormat(op_network, false);
        var (_,_,qsharp_hamiltonian) = H.ToPauliHamiltonian().ToQSharpFormat();

        Assert.Equal(qsharp_swap_network.Length+1, qsharp_data.Length);
        using (var qsim = new QuantumSimulator())
        {
            SwapNetworkOneSummandTestOp.Run(qsim, qsharp_swap_network, qsharp_data, qsharp_hamiltonian,
                                    (long)start_order.Length)
                             .Wait();
        }
    }

    public static IEnumerable<object[]> OneTermHamiltonianData() {
        var result = new List<object[]> {};
        var num_sites = 5;
        for (int i = 0; i < num_sites; i++) {
            for (int j = i+1; j < num_sites; j++) {
                var H = new FermionHamiltonian {};
                H.Add(new HermitianFermionTerm(new int[] {i, j}), (double)(10*i+j));
                var swap_network = OneBodyDenseNetwork(num_sites);
                result.Add(new object[] {H, swap_network, Enumerable.Range(0,num_sites).ToArray()});
        }
        }
    
        return result;
    }

    [Theory]
    [MemberData(nameof(HamiltonianData))]
    public void TestHamiltonian(
        FermionHamiltonian H,
        int num_sites,
        SwapNetwork swap_network,
        double step_size,
        double time
    ) {
        var start_order = Enumerable.Range(0,num_sites).ToArray();
        var (op_network, end_order) = TrotterStepData(H, swap_network, start_order);
        //we use 32 bit ints until the point of injection into q#, which requires 64 bit ints.
        var qsharp_swap_network = swap_network.ToQSharpFormat();
        var qsharp_data = ToQSharpFormat(op_network, false);
        var (_,_,qsharp_hamiltonian) = H.ToPauliHamiltonian().ToQSharpFormat();

        Assert.Equal(qsharp_swap_network.Length+1, qsharp_data.Length);
        using (var qsim = new QuantumSimulator())
        {
            SwapNetworkEvolutionTestOp.Run(qsim, qsharp_swap_network, qsharp_data, qsharp_hamiltonian,
                                    (long)num_sites, step_size, time)
                             .Wait();
        }        
    }

    // One term hamiltonians, similar to previous test
    public static IEnumerable<object[]> HamiltonianData() {
        var result = new List<object[]> {};
        var num_sites = 5;
        for (int i = 0; i < num_sites; i++) {
            for (int j = i+1; j < num_sites; j++) {
                var H = new FermionHamiltonian {};
                H.Add(new HermitianFermionTerm(new int[] {i, j}), (double)(10*i+j));
                var swap_network = OneBodyDenseNetwork(num_sites);
                var step_size = 1;
                var time = 2;
                result.Add(new object[] {H, num_sites, swap_network, step_size, time});
            }
        }
        // A dense hamiltonian
        var H2 = new FermionHamiltonian {};
        var swap_network2 = OneBodyDenseNetwork(num_sites);
        var step_size2 = .0005;
        var time2 = .1;
        for (int i = 0; i < num_sites; i++) {
            for (int j = i+1; j < num_sites; j++) {
                H2.Add(new HermitianFermionTerm(new int[] {i, j}), (double)(10*i+j));
            }
        }
        result.Add(new object[] {H2, num_sites, swap_network2, step_size2, time2});

    
        return result;
    }

    public string layers_string(SwapNetwork swaps) {
            var result = "{";
            var swaps_occupied = false;
            foreach (var layer in swaps) {
                if (swaps_occupied) {
                    result += ", ";
                } else {
                    result += "{";
                    swaps_occupied = true;
                }
                var layer_occupied = false;
                foreach (var (a,b) in layer) {
                    if (layer_occupied) {
                        result += ", ";
                    } else {
                        result += "{";
                        layer_occupied = true;
                    }
                    result += $"({a}, {b})";
                }
                result += "}";
            }
            result += "}";
            return result;
        }

    }
}