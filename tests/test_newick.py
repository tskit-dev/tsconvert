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
import itertools

import dendropy
import msprime
import numpy as np
import pytest
import tskit

import tsconvert


def kc_distance(tree1, tree2):
    """
    Returns the Kendall-Colijn topological distance between the specified
    pair of trees. This is a very simple and direct implementation for testing.
    """
    samples = tree1.tree_sequence.samples()
    if not np.array_equal(samples, tree2.tree_sequence.samples()):
        raise ValueError("Trees must have the same samples")
    k = samples.shape[0]
    n = (k * (k - 1)) // 2
    M = [np.ones(n + k), np.ones(n + k)]
    for tree_index, tree in enumerate([tree1, tree2]):
        stack = [(tree.root, 0)]
        while len(stack) > 0:
            u, depth = stack.pop()
            children = tree.children(u)
            for v in children:
                stack.append((v, depth + 1))
            for c1, c2 in itertools.combinations(children, 2):
                for u in tree.samples(c1):
                    for v in tree.samples(c2):
                        if u < v:
                            a, b = u, v
                        else:
                            a, b = v, u
                        pair_index = a * (a - 2 * k + 1) // -2 + b - a - 1
                        assert M[tree_index][pair_index] == 1
                        M[tree_index][pair_index] = depth
    return np.linalg.norm(M[0] - M[1])


def get_nonbinary_example(sample_size=20, recombination_rate=0, random_seed=42):
    ts = msprime.simulate(
        sample_size=sample_size,
        recombination_rate=recombination_rate,
        random_seed=random_seed,
        demographic_events=[
            msprime.SimpleBottleneck(time=0.5, population=0, proportion=1)
        ],
    )
    # Make sure this really has some non-binary nodes
    found = False
    for e in ts.edgesets():
        if len(e.children) > 2:
            found = True
            break
    assert found
    return ts


class TestSingleTreeRoundTrip:
    """
    Tests that we can successfully roundtrip trees with various topologies
    from a tree sequence containing a single tree.
    """

    def verify(self, ts):
        assert ts.num_trees == 1
        source_tree = ts.first()
        newick = tsconvert.to_newick(source_tree)
        conv_ts = tsconvert.from_newick(newick)

        # Round tripping can introduce very slightly differing node times
        # due to floating point manipulations, so round the times
        conv_tables = conv_ts.dump_tables()
        nodes = conv_tables.nodes.copy()
        conv_tables.nodes.clear()
        for node in nodes:
            conv_tables.nodes.append(node.replace(time=round(node.time, 15)))
        conv_tables.sort()
        conv_ts = conv_tables.tree_sequence()

        assert conv_ts.num_trees == 1
        conv_tree = conv_ts.first()
        source_str = source_tree.draw(format="unicode", node_labels={}, order="tree")
        conv_str = conv_tree.draw(format="unicode", node_labels={}, order="tree")
        # The tree sequences are canonical, so the nodes are allocated in
        # time order. We should be identical other than the leaf labels.
        assert source_str == conv_str
        assert np.allclose(
            sorted(conv_ts.tables.nodes.time), sorted(ts.tables.nodes.time)
        )

    def test_msprime_binary(self):
        self.verify(msprime.simulate(10, random_seed=1))

    def test_msprime_non_binary(self):
        self.verify(get_nonbinary_example(8))


class TestMsRoundTrip:
    """
    Tests if we can round trip tree sequences through the ms format.
    """

    def verify(self, ts):
        msout = tsconvert.to_ms(ts.simplify())
        new_ts = tsconvert.from_ms(msout)
        assert ts.num_trees == new_ts.num_trees
        for t1, t2 in zip(ts.trees(), new_ts.trees()):
            assert t1.interval[0] == pytest.approx(t2.interval[0])
            assert t1.interval[1] == pytest.approx(t2.interval[1])
            assert kc_distance(t1, t2) == 0

    def test_msprime_single_tree(self):
        self.verify(msprime.simulate(10, random_seed=12))

    def test_msprime_binary(self):
        self.verify(msprime.simulate(10, recombination_rate=1, random_seed=1))

    def test_msprime_non_binary(self):
        ts = get_nonbinary_example(8, recombination_rate=1)
        assert ts.num_trees > 1
        self.verify(ts)

    # TODO more examples


class TestFromMs:
    """
    Tests for the from_ms function.
    """

    def test_empty_input(self):
        with pytest.raises(ValueError):
            tsconvert.from_ms("")

    def test_ms_without_trees_flag(self):
        msout = """ms 3 1 -t 4 -r 5 6 -seeds 1 2 3
        1 2 3

        //
        segsites: 3
        positions: 0.3054 0.3812 0.5338
        111
        000
        100
        """
        with pytest.raises(ValueError):
            tsconvert.from_ms(msout)

    def test_empty_trees(self):
        msout = """
        [];
        """
        with pytest.raises(ValueError):
            tsconvert.from_ms(msout)

        msout = """
        [1];
        """
        with pytest.raises(ValueError):
            tsconvert.from_ms(msout)

    def test_ms_without_recombination(self):
        msout = """
        ms 4 1 -t 5 -T -seeds 1 2 3
        1 2 3

        //
        ((2:0.0680,3:0.0680):0.1481,(1:0.2124,4:0.2124):0.0038);
        """
        with pytest.raises(ValueError):
            tsconvert.from_ms(msout)

    def test_malformed_length(self):
        msout = """
        5(1:0.27413282187548,2:0.27413282187548);
        """
        with pytest.raises(ValueError):
            tsconvert.from_ms(msout)

        msout = """
        [XXX](1:0.27413282187548,2:0.27413282187548);
        """
        with pytest.raises(ValueError):
            tsconvert.from_ms(msout)

        msout = """
        [](1:0.27413282187548,2:0.27413282187548);
        """
        with pytest.raises(ValueError):
            tsconvert.from_ms(msout)

        msout = """
        [5(1:0.27413282187548,2:0.27413282187548);
        """
        with pytest.raises(dendropy.utility.error.DataParseError):
            tsconvert.from_ms(msout)

    def test_nonmatching_tips(self):
        msout = """
        [2](1:0.2144,3:0.2144);
        [4](3:0.2144,(1:0.0768,2:0.0768):0.1376);
        """
        with pytest.raises(ValueError):
            tsconvert.from_ms(msout)

        msout = """
        [2](2:0.2930,(1:0.2144,4:0.2144):0.0786);
        [4](3:0.2144,(1:0.0768,2:0.0768):0.1376);
        """
        with pytest.raises(ValueError):
            tsconvert.from_ms(msout)

    def test_identical_node_times(self):
        msout = """
        [2](((1:1.5,2:1.5):1.7,3:3.2):1.1,(4:1.5,5:1.5):2.8);
        """
        with pytest.raises(ValueError):
            tsconvert.from_ms(msout)

    def test_zero_sequence_length(self):
        msout = """
        [0](1:0.27413282187548,2:0.27413282187548);
        """
        with pytest.raises(tskit.TskitException):
            tsconvert.from_ms(msout)

    def test_bad_edges(self):
        msout = """
        [1](1:0.27413282187548,2:0.27413282187548);
        [0](1:0.27413282187548,2:0.27413282187548);
        """
        with pytest.raises(tskit.TskitException):
            tsconvert.from_ms(msout)

    def test_single_tree(self):
        msout = """
        [1](1:0.27413282187548,2:0.27413282187548);
        """
        ts = tsconvert.from_ms(msout)
        assert ts.num_trees == 1

    def test_full_ms_output(self):
        msout = """
        ms 3 1 -t 4 -r 5 6 -seeds 1 2 3 -T
        1 2 3

        //
        [2](2:0.2930,(1:0.2144,3:0.2144):0.0786);
        [4](3:0.2144,(1:0.0768,2:0.0768):0.1376);
        segsites: 3
        positions: 0.3054 0.3812 0.5338
        111
        000
        100
        """
        ts = tsconvert.from_ms(msout)
        assert ts.num_samples == 3
        assert ts.sequence_length == 6
        assert ts.num_trees == 2
        assert ts.num_nodes == 6

        trees = ts.trees()
        tree = next(trees)
        assert tree.interval == (0, 2)
        internal_nodes = set(tree.nodes()) - set(ts.samples())
        assert tree.branch_length(0) == pytest.approx(0.2144)
        assert tree.branch_length(1) == pytest.approx(0.2930)
        assert tree.branch_length(2) == pytest.approx(0.2144)
        assert tree.branch_length(min(internal_nodes)) == pytest.approx(0.0786)

        tree = next(trees)
        assert tree.interval == (2, 6)
        internal_nodes = set(tree.nodes()) - set(ts.samples())
        assert tree.branch_length(0) == pytest.approx(0.0768)
        assert tree.branch_length(1) == pytest.approx(0.0768)
        assert tree.branch_length(2) == pytest.approx(0.2144)
        assert tree.branch_length(min(internal_nodes)) == pytest.approx(0.1376)

    def test_equal_internal_node_time(self):
        #     6
        #   ┏━┻━┓
        #   4   5
        #  ┏┻┓ ┏┻┓
        #  0 1 2 3
        tables = tskit.TableCollection(1)
        for _ in range(4):
            tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
        tables.nodes.add_row(0, time=1)
        tables.nodes.add_row(0, time=1)
        tables.nodes.add_row(0, time=2)

        tables.edges.add_row(0, 1, 4, 0)
        tables.edges.add_row(0, 1, 4, 1)
        tables.edges.add_row(0, 1, 5, 2)
        tables.edges.add_row(0, 1, 5, 3)
        tables.edges.add_row(0, 1, 6, 4)
        tables.edges.add_row(0, 1, 6, 5)
        tables.sort()
        ts = tables.tree_sequence()
        msout = tsconvert.to_ms(ts)
        # The current algorithm assumes node times are unique
        with pytest.raises(ValueError):
            tsconvert.from_ms(msout)

    def test_n2_example(self):
        # j$ mspms 2 1 -T -r 4 10 -p 14
        msout = """
        /home/jk/.local/bin/mspms 2 1 -T -r 4 10 -p 14
        1774173160 1383299789 1436475231

        //
        [5](1:0.27413282187548,2:0.27413282187548);
        [3](1:0.43103605328988,2:0.43103605328988);
        [1](1:1.96842212024363,2:1.96842212024363);
        [1](1:2.06027985820196,2:2.06027985820196);
        """
        ts = tsconvert.from_ms(msout)
        assert ts.num_samples == 2
        assert ts.sequence_length == 10
        assert ts.num_trees == 4

        trees = ts.trees()
        tree = next(trees)
        assert tree.interval == (0, 5)
        assert tree.branch_length(0) == pytest.approx(0.27413282187548)
        assert tree.branch_length(1) == pytest.approx(0.27413282187548)

        tree = next(trees)
        assert tree.interval == (5, 8)
        assert tree.branch_length(0) == pytest.approx(0.43103605328988)
        assert tree.branch_length(1) == pytest.approx(0.43103605328988)

        tree = next(trees)
        assert tree.interval == (8, 9)
        assert tree.branch_length(0) == pytest.approx(1.96842212024363)
        assert tree.branch_length(1) == pytest.approx(1.96842212024363)

        tree = next(trees)
        assert tree.interval == (9, 10)
        assert tree.branch_length(0) == pytest.approx(2.06027985820196)
        assert tree.branch_length(1) == pytest.approx(2.06027985820196)

    def test_n4_example(self):
        # $ mspms 4 1 -T -r 4 10 -p 8
        msout = """
        /home/jk/.local/bin/mspms 4 1 -T -r 4 10 -p 8
        961626313 1881970557 110898863
        
        //
        [5](1:0.70961771,(4:0.33536000,(2:0.12737966,3:0.12737966):0.20798034):0.37425772);
        [1]((2:0.12737966,3:0.12737966):0.20798034,(1:0.21249950,4:0.21249950):0.12286050);
        [2]((3:0.12737966,(2:0.02380236,4:0.02380236):0.10357730):0.20798034,1:0.33536000);
        [1](1:1.32624987,(3:0.12737966,(2:0.02380236,4:0.02380236):0.10357730):1.19887022);
        [2](1:1.80041212,(3:0.12737966,(2:0.02380236,4:0.02380236):0.10357730):1.67303246);
        """  # noqa: B950 W293
        ts = tsconvert.from_ms(msout)
        assert ts.num_samples == 4
        assert ts.sequence_length == 11
        assert ts.num_trees == 5

        trees = ts.trees()
        tree = next(trees)
        assert tree.interval == (0, 5)
        assert tree.branch_length(0) == pytest.approx(0.70961771)
        assert tree.branch_length(1) == pytest.approx(0.12737966)
        assert tree.branch_length(2) == pytest.approx(0.12737966)
        assert tree.branch_length(3) == pytest.approx(0.33536000)

        tree = next(trees)
        assert tree.interval == (5, 6)
        assert tree.branch_length(0) == pytest.approx(0.21249950)
        assert tree.branch_length(1) == pytest.approx(0.12737966)
        assert tree.branch_length(2) == pytest.approx(0.12737966)
        assert tree.branch_length(3) == pytest.approx(0.21249950)

        tree = next(trees)
        assert tree.interval == (6, 8)

        tree = next(trees)
        assert tree.interval == (8, 9)

        tree = next(trees)
        assert tree.interval == (9, 11)
        assert tree.branch_length(0) == pytest.approx(1.80041212)
        assert tree.branch_length(1) == pytest.approx(0.02380236)
        assert tree.branch_length(2) == pytest.approx(0.12737966)
        assert tree.branch_length(3) == pytest.approx(0.02380236)

        # Should get the same output if we strip off the header stuff.
        msout = """
        [5](1:0.70961771,(4:0.33536000,(2:0.12737966,3:0.12737966):0.20798034):0.37425772);
        [1]((2:0.12737966,3:0.12737966):0.20798034,(1:0.21249950,4:0.21249950):0.12286050);
        [2]((3:0.12737966,(2:0.02380236,4:0.02380236):0.10357730):0.20798034,1:0.33536000);
        [1](1:1.32624987,(3:0.12737966,(2:0.02380236,4:0.02380236):0.10357730):1.19887022);
        [2](1:1.80041212,(3:0.12737966,(2:0.02380236,4:0.02380236):0.10357730):1.67303246);
        """  # noqa: B950
        tables = tsconvert.from_ms(msout).tables
        tables.assert_equals(ts.tables, ignore_timestamps=True)


class TestNewicks:
    def test_non_ultrametric(self):
        newick = "((5:0.10,6:0.20):5.34,(7:1.00,8:0.01):0.02,4:42.42);"
        assert tsconvert.from_newick(newick).first().newick(precision=2) == newick

    def test_empty(self):
        with pytest.raises(ValueError, match="Newick string was empty"):
            tsconvert.from_newick("")

    def test_multitree(self):
        with pytest.raises(
            ValueError, match="Only one tree can be imported from a newick string"
        ):
            tsconvert.from_newick("(2:10.00,3:10.00);(2:10.00,3:10.00);")

    def test_annotation_ignored(self):
        assert (
            tsconvert.from_newick("(2:0.10,[IGNOREME]3:0.10);")
            .first()
            .newick(precision=2)
            == "(2:0.10,3:0.10);"
        )

    def test_root_time(self):
        ts = tsconvert.from_newick("(2:0.10,3:0.20):42.00;")
        # Root times are removed from tskit's newick output:
        assert ts.first().newick(precision=2) == "(2:0.10,3:0.20);"
        # But should still be in the ts:
        assert np.array_equal(ts.tables.nodes.time, [42.2, 42.1, 42.0])

    def test_zero_branch_length(self):
        with pytest.raises(
            ValueError,
            match="tskit tree sequences cannot contain edges" " with lengths <= 0",
        ):
            tsconvert.from_newick("(2:0.10,3:0.00);")

        ts = tsconvert.from_newick("(2:0.10,3:0.00);", 0.01)
        assert ts.first().newick(precision=2) == "(2:0.10,3:0.01);"

    def test_negative_branch_length(self):
        with pytest.raises(
            ValueError,
            match="tskit tree sequences cannot contain edges" " with lengths <= 0",
        ):
            tsconvert.from_newick("(2:0.10,3:0.00);")

        ts = tsconvert.from_newick("(2:0.10,3:-10.01);", 0.01)
        assert ts.first().newick(precision=2) == "(2:0.10,3:0.01);"

    def test_ids(self):
        ts = tsconvert.from_newick("(BILL:0.10,:0.30)LUCY;")
        assert ts.first().newick(precision=2) == "(2:0.10,3:0.30);"
        assert tuple(n.metadata for n in ts.nodes()) == (
            {"name": "LUCY"},
            {"name": "BILL"},
            {},
        )

    def test_comments(self):
        ts = tsconvert.from_newick(
            '(BILL[42]:0.10,[!"£$%^&*_+-={}<>,.?/~#|`¬]:0.30)LUCY;'
        )
        assert ts.first().newick(precision=2) == "(2:0.10,3:0.30);"
        assert tuple(n.metadata for n in ts.nodes()) == (
            {
                "name": "LUCY",
            },
            {"name": "BILL", "comment": "42"},
            {"comment": '!"£$%^&*_+-={}<>,.?/~#|`¬'},
        )
