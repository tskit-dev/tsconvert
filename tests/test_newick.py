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

import unittest

import msprime
import numpy as np

import tsconvert


def get_nonbinary_example(sample_size=20, recombination_rate=0, random_seed=42):
    ts = msprime.simulate(
        sample_size=sample_size, recombination_rate=recombination_rate,
        random_seed=random_seed,
        demographic_events=[
            msprime.SimpleBottleneck(time=0.5, population=0, proportion=1)])
    # Make sure this really has some non-binary nodes
    found = False
    for e in ts.edgesets():
        if len(e.children) > 2:
            found = True
            break
    assert found
    return ts


class TestSingleTreeRoundTrip(unittest.TestCase):
    """
    Tests that we can successfully roundtrip trees with various topologies
    from a tree sequence containing a single tree.
    """
    def verify(self, ts):
        self.assertEqual(ts.num_trees, 1)
        source_tree = ts.first()
        newick = tsconvert.to_newick(source_tree)
        conv_ts = tsconvert.from_newick(newick)
        self.assertEqual(conv_ts.num_trees, 1)
        conv_tree = conv_ts.first()
        source_str = source_tree.draw(format="unicode", node_labels={})
        conv_str = conv_tree.draw(format="unicode", node_labels={})
        # The tree sequences are canonical, so the nodes are allocated in
        # time order. We should be identical other than the leaf labels.
        # There are surely more efficient ways of comparing topologies though.
        self.assertEqual(source_str, conv_str)
        self.assertTrue(
            np.allclose(conv_ts.tables.nodes.time, ts.tables.nodes.time))
        self.assertEqual(list(range(ts.num_samples)), list(conv_tree.leaves()))
        self.assertEqual(list(range(ts.num_samples)), list(conv_tree.samples()))

    def test_msprime_binary(self):
        self.verify(msprime.simulate(10, random_seed=1))

    def test_msprime_non_binary(self):
        self.verify(get_nonbinary_example(8))
