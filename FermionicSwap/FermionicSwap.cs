// Copyright Battelle Memorial Institute 2022. All rights reserved.

using System.Linq;
using System.Collections.Immutable;

using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.LadderOperators;
using Microsoft.Quantum.Chemistry;

namespace FermionicSwap
{


    using System;

    public static class FSTools {
        // fixme: What is the naming convention for a library that includes a
        // single static class?

        // fixme: Is this the correct way to specify type aliases that are
        // globally accessible? I originally was using `using` directives but
        // they seem to have some quirky limitations. This way, am I stuck
        // providing boilerplate because constructors aren't inherited?

        public class FermionOperator : LadderOperator<int> {
            public FermionOperator(RaisingLowering t, int index):base(t,index) {}
        }
        public class SwapLayer : List<(int,int)> {}
        public class SwapNetwork : List<List<(int,int)>> {}
        public class OperatorLayer : List<(HermitianFermionTerm,DoubleCoeff)> {}
        public class OperatorNetwork : List<OperatorLayer> {}
        public class TermsDictionaryComparer : IEqualityComparer<ImmutableArray<int>>
        {
            public bool Equals(ImmutableArray<int> x, ImmutableArray<int> y)
            {
                return x!.SequenceEqual(y!);
            }

            public int GetHashCode(ImmutableArray<int> obj)
            {
                return obj.Aggregate(0, (ob1,ob2) => HashCode.Combine(ob1,ob2));
            }
        }
        public class TermsDictionary : Dictionary<ImmutableArray<int>, List<(HermitianFermionTerm,DoubleCoeff)>> {
            public TermsDictionary() : base(new TermsDictionaryComparer {}) {}
        }

        /// # Summary
        /// Returns the first swap layer of the swap network that takes
        /// elements of start_order to desired_positions
        ///
        /// # Input
        /// ## start_order
        /// A position-indexed array indicating the index of the site orbital corresponding to each position.
        /// ## desired_positions
        /// A map from orbital sites to desired final positions
        /// ## even_parity
        /// True if this layer should consist of even-odd swaps, false if odd-even swaps.
        ///
        /// # Output
        /// A tuple (next_order, layer) consisting of
        /// ## next_order
        /// A position indexed array indicating site orbital positioning after the swap layer is applied 
        /// ## layer
        /// A list of mutually disjoint (n,n+1) transpositions.
        private static (int[],SwapLayer) evenOddSwapLayer(int[] start_order, Dictionary<int,int> desired_positions, bool even_parity) {

            int start = even_parity?0:1;
            var new_order = start_order.ToArray();
            var swaps = new SwapLayer();

            for(int i=start;i<start_order.Length-1; i+=2) {
                if (desired_positions[start_order[i]] < desired_positions[start_order[i+1]]) {
                    new_order[i] = start_order[i];
                    new_order[i+1] = start_order[i+1];
                } else {
                    new_order[i+1] = start_order[i];
                    new_order[i] = start_order[i+1];
                    swaps.Add((i,i+1));            
                }
            }
            Console.Write("new order: " + string.Join(", ", new_order) + "\n");
            return (new_order, swaps);
        }

        /// # Summary
        /// Return a map from site orbital indices to positions, given a position-indexed array of site orbital indices.
        ///
        /// # Input
        /// ## order
        /// A position-indexed array indicating the index of the site orbital corresponding to each position.
        ///
        /// # Output
        /// A map from site orbital indices to position indices.
        public static Dictionary<int,int> PositionDictionary(int[] order) {
            var desired_positions = new Dictionary<int,int>();
                for (int i = 0; i < order.Length; i++)
                {
                    desired_positions[order[i]] = i;
                }
            return desired_positions;
        
        }

        /// # Summary
        /// Return a network of swaps that converts an initial ordering of site
        /// orbitals to the desired final ordering, by implementing an even-odd
        /// sort algorithm.
        ///
        /// Even-odd sort produces a network with minimal swaps and at most one
        /// more than minimal circuit depth. A greedy circuit-packing algorithm
        /// will eliminate the one-extra circuit depth.
        ///
        /// # Input
        /// ## start_order
        /// A position-indexed array of site orbital indices, indicating their
        /// starting order.
        /// ## end_order
        /// A position-indexed array of site orbital indices, indicating their
        /// desired order
        ///
        /// # Output
        /// A SwapNetwork describing layers of disjoint (n,n+1) transpositions
        /// which convert start_order to end_order.
        public static SwapNetwork evenOddSwap(int[] start_order, int[] end_order) {
            var result = new SwapNetwork();
            var this_order = start_order;

            bool at_least_once = false;
            bool done = false;
            bool even_parity = true;
            var desired_positions = PositionDictionary(end_order);
            while (!done) {
                var (next_order,swaps) = evenOddSwapLayer(this_order, desired_positions, even_parity);
                this_order = next_order;
                if (swaps.Count > 0) {
                    result.Add(swaps);
                    even_parity = !even_parity;
                    at_least_once = true;
                } else {
                    if (at_least_once) {
                        done = true;
                    } else {
                        at_least_once = true;
                        even_parity = !even_parity;
                    }

                }
            }
            return result;
        }

        /// # Summary
        /// Return a swap network suitable for evaluating a Trotter step for a
        /// one-body dense Hamiltonian. The resulting network fully reverses
        /// the order of the site-orbitals. The Trotter step can be evolved
        /// without using the last two of the swap layers, but we do not
        /// assume that optimization here.
        ///
        /// # Input
        /// ## num_sites
        /// The number of site orbitals in the Hamiltonian.
        /// # Output
        /// A swap network that reverses the order of the site orbitals. 
        public static SwapNetwork OneBodyDenseNetwork(int num_sites) {
            var start_order = Enumerable.Range(0,num_sites).ToArray();
            var final_order = start_order.Select(x => num_sites - x-1).ToArray();
            return evenOddSwap(start_order, final_order);
        }

        /// # Summary
        /// Return an n-body fermionic Hamiltonian term, re-indexed to be
        /// evaluated in the specified Jordan-Wigner ordering.
        ///
        /// # Input
        /// ## term
        /// A Hamiltonian term (and by implication, its Hermitian conjugate).
        /// ## actual positions
        /// A map from site orbital indices to positions, specifying a
        /// Jordan-Wigner ordering
        ///
        /// # Output
        /// A new Hamiltonian term, with reordered indices and possibly
        /// opposite sign.
        /// 
        /// ## Note
        /// Because we return the reordered sequence as a FermionTerm, the
        /// listed order of the operators will be shuffled (and the sign
        /// adjusted) to match QDK's canonical order. Because it is a
        /// HermitionFermionTerm, the reordered sequence gets replaced with its
        /// adjoint when that results in lower the canonical ordering.
        public static HermitianFermionTerm ReorderedFermionTerm(HermitianFermionTerm term, Dictionary<int,int> actual_positions) {
            return new HermitianFermionTerm(term.Sequence.Select(o=>new FermionOperator(o.Type, actual_positions[o.Index])),
                            term.Coefficient);
        }

        /// # Summary
        /// Return a plan which Q# can use to evaluate a Trotter step for a
        /// given fermionic swap network, evolving Jordan-Wigner-reordered
        /// terms between swap layers, as they become local.
        ///
        /// # Input
        /// ## H
        /// A Hamiltonian.
        /// ## swap_network
        /// The network of swaps to be applied.
        /// ## start_order
        /// A position-indexed array of site orbital positions, indicating
        /// the Jordan-Wigner ordering prior to any swaps being performed.
        /// ## time
        /// The amount of evolution time for each Hamiltonian term.
        ///
        /// # Output
        /// A tuple (network, end_order), consisting of the following:
        /// ## network
        /// A list, one layer longer than swap_network, containing the local
        /// operators to evaluate between each swap layer. Operators are
        /// ordered so that a greedy circuit-packing algorithm will produce
        /// a reasonably low-depth circuit.
        /// ## end_order
        /// A position-indexed list of site orbital indices, indicating the
        /// Jordan-Wigner ordering after all swaps are performed. 
        public static (OperatorNetwork, int[] ) TrotterStepData(
            FermionHamiltonian H,
            SwapNetwork swap_network,
            int[] start_order,
            DoubleCoeff time
            )
        {
            var op_network = new OperatorNetwork {};
            var end_order = start_order.ToArray();
            var terms = new TermsDictionary();
            foreach (var (term_type, term_list) in H.Terms) {
                foreach (var (term,term_value) in term_list) {
                    var term_sites = ImmutableArray.Create(
                        term.Sequence.OrderBy(o => o.Index).Select(o => o.Index).Distinct().ToArray()
                    );
                    if (!terms.ContainsKey(term_sites)) {
                        terms[term_sites] = new OperatorLayer {};
                    }
                    terms[term_sites].Add((term,term_value));
                }
            }

            op_network.Add(ProcessNetworkLayer(terms, end_order, time));
            // Apply each swap layer to the ordering and add layer interactions
            foreach (var layer in swap_network) {
                foreach (var (oldpos,newpos) in layer) {
                    (end_order[oldpos], end_order[newpos]) = (end_order[newpos], end_order[oldpos]);
                }
                op_network.Add(ProcessNetworkLayer(terms, end_order, time));
            }
            return (op_network, end_order);
        }

        /// # Summary
        /// Return a layer of operators for a network, updating the dictionary
        /// of already-applied terms.
        ///
        /// # Input
        /// ## term_dict
        /// A map from ordered indexed lists to lists of unapplied terms having
        /// those indices.
        /// ## order
        /// A position-indexed array of site orbital indices, indicating the
        /// current Jordan-Wigner ordering.
        /// ## time
        /// The amount of time to evolve each term.
        ///
        /// # Output
        /// A list of local HermetianFermionTerms to evaluate, expressed in
        /// the local Jordan-Wigner ordering.
        ///
        /// # Side effects
        /// Terms are removed from term_dict as they are applied.
        public static OperatorLayer ProcessNetworkLayer(
            TermsDictionary term_dict,
            int[] order,
            DoubleCoeff time)
        {
            var result = new OperatorLayer {};
            bool productive = true;
            while (productive) {
                productive = false;
                for (var start = 0; start < order.Count(); start++) {
                    for (var end = start+1; end <= order.Count(); end++) {
                        var key = ImmutableArray.Create(order[start..end].OrderBy(o=>o).Distinct().ToArray());
                        if (term_dict.ContainsKey(key)) {
                            var (term,coeff) = term_dict[key][0];
                            result.Add((ReorderedFermionTerm(term, PositionDictionary(order[0..end])), time*coeff));
                            term_dict[key].RemoveAt(0);
                            if (term_dict[key].Count() == 0) {
                                term_dict.Remove(key);
                            }
                            productive = true;
                            (start,end) = (end,end+1); // try to greedily add parallel evolutions
                        }
                    }
                }
            }
            return result;
        }
    }
}