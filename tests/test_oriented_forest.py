#
# MIT License
#
# Copyright (c) 2019 Tskit Developers
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
import discsim
import ercs
import tskit

import tsconvert


class OrientedForestBaseTest:
    def verify_tree(self, tree, n, pi, tau):
        """
        Verifies the specified oriented tree representation is equivalent to the
        specified tskit tree.
        """
        assert sorted(tree.samples()) == list(range(n))
        for j in range(n):
            p1 = []
            u = j
            while u != tskit.NULL:
                p1.append(tree.time(u))
                u = tree.parent(u)
            p2 = []
            # Nodes are numbered from 1 in oriented forests.
            u = j + 1
            while u != 0:
                p2.append(tau[u])
                u = pi[u]
            assert p1 == p2


class TestSingleTree(OrientedForestBaseTest):
    """
    Tests that we can extract an oriented tree topology from a single locus simulation.
    """

    def verify(self, n, pi, tau):
        ts = tsconvert.from_oriented_forest(n, pi, tau)
        assert ts.num_trees == 1
        tree = ts.first()
        assert ts.num_samples == n
        self.verify_tree(tree, n, pi[0], tau[0])

    def test_discsim_n3(self):
        sim = discsim.Simulator(10)
        sim.sample = [None, (3, 2), (6, 4), (7, 0)]
        sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
        sim.run()
        pi, tau = sim.get_history()
        self.verify(3, pi, tau)

    def test_ercs_n3(self):
        sim = ercs.Simulator(10)
        sim.sample = [None, (3, 2), (6, 4), (7, 0)]
        sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
        pi, tau = sim.run(1)
        self.verify(3, pi, tau)

    def test_discsim_n10(self):
        sim = discsim.Simulator(10)
        sim.sample = [None] + [(j, j) for j in range(10)]
        sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
        sim.run()
        pi, tau = sim.get_history()
        self.verify(10, pi, tau)

    def test_ercs_n10(self):
        sim = ercs.Simulator(10)
        sim.sample = [None] + [(j, j) for j in range(10)]
        sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
        pi, tau = sim.run(1)
        self.verify(10, pi, tau)

    def test_discsim_n10_nonbinary(self):
        sim = discsim.Simulator(10)
        sim.sample = [None] + [(j, j) for j in range(10)]
        sim.event_classes = [ercs.DiscEventClass(u=0.99, r=2.5)]
        sim.run()
        pi, tau = sim.get_history()
        self.verify(10, pi, tau)

    def test_ercs_n10_nonbinary(self):
        sim = ercs.Simulator(10)
        sim.sample = [None] + [(j, j) for j in range(10)]
        sim.event_classes = [ercs.DiscEventClass(u=0.99, r=5)]
        pi, tau = sim.run(1)
        self.verify(10, pi, tau)


class TestManyTrees(OrientedForestBaseTest):
    """
    Tests that we can extract an oriented tree topology from a multi locus simulation.
    """

    def verify(self, n, pi, tau):
        ts = tsconvert.from_oriented_forest(n, pi, tau)
        assert len(tau) == len(pi)
        assert ts.sequence_length == len(pi)
        assert ts.num_samples == n
        for tree in ts.trees():
            for i in range(*map(int, tree.interval)):
                self.verify_tree(tree, n, pi[i], tau[i])

    def test_discsim_n3(self):
        sim = discsim.Simulator(10)
        sim.sample = [None, (3, 2), (6, 4), (7, 0)]
        sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
        sim.run()
        pi, tau = sim.get_history()
        self.verify(3, pi, tau)

    def test_ercs_n3(self):
        sim = ercs.Simulator(10)
        sim.sample = [None, (3, 2), (6, 4), (7, 0)]
        sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
        num_loci = 10
        sim.recombination_probabilities = [0.1] * (num_loci - 1)
        pi, tau = sim.run(1)
        assert len(pi) == num_loci
        self.verify(3, pi, tau)

    def test_discsim_n10(self):
        sim = discsim.Simulator(10)
        sim.sample = [None] + [(j, j) for j in range(10)]
        sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
        sim.run()
        pi, tau = sim.get_history()
        self.verify(10, pi, tau)

    def test_ercs_n10(self):
        sim = ercs.Simulator(10)
        sim.sample = [None] + [(j, j) for j in range(10)]
        sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
        num_loci = 3
        sim.recombination_probabilities = [0.1] * (num_loci - 1)
        pi, tau = sim.run(1)
        assert len(pi) == num_loci
        self.verify(10, pi, tau)

    def test_discsim_n10_nonbinary(self):
        sim = discsim.Simulator(10)
        sim.sample = [None] + [(j, j) for j in range(10)]
        sim.event_classes = [ercs.DiscEventClass(u=0.99, r=2.5)]
        sim.run()
        pi, tau = sim.get_history()
        self.verify(10, pi, tau)

    def test_ercs_n10_nonbinary(self):
        sim = ercs.Simulator(10)
        sim.sample = [None] + [(j, j) for j in range(10)]
        sim.event_classes = [ercs.DiscEventClass(u=1, r=5)]
        num_loci = 3
        sim.recombination_probabilities = [0.1] * (num_loci - 1)
        pi, tau = sim.run(1)
        assert len(pi) == num_loci
        self.verify(10, pi, tau)
