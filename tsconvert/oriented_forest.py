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
"""
Converts trees specified in Oriented forest form to tskit. This is the output
produced by the discsim and ercs simulators.
"""
import tskit


def from_oriented_forest(n, pi, tau):
    """
    Returns a tree sequence from the specified oriented forest definition.
    Oriented forests are returned by the
    `discsim <https://github.com/jeromekelleher/discsim>`_ and
    `ercs <https://github.com/jeromekelleher/ercs>`_ simulators, which
    simulate models of evolution in a spatial continuum.

    An oriented forest is defined by two lists, pi and tau. These lists
    defined the trees at L loci, such that pi[l] and tau[l] define the
    topology and times at locus l. For a locus l, pi[l][j] gives the parent
    of node j and tau[l][j] is the time of node j. The samples are nodes
    1 to n (inclusive), and node 0 denotes the root (pi[l][0] and pi[l][1]
    are arbitrary). The converted tree sequence contains each of these
    trees with the same topology and times, but with sample nodes remapped
    to 0 to n - 1.

    :param int n: The number of samples.
    :param list(list) pi: The list of per-locus tree topologies as a parent
        array.
    :param list(list) tau: The list of per-locus node times.
    :rtype: tskit.TreeSquence
    :return: A TreeSequence object with the same topologies as the specified
        oriented forest.
    """
    L = len(pi)
    tables = tskit.TableCollection(L)
    # Allocate the samples
    for j in range(n):
        tables.nodes.add_row(time=tau[0][j + 1], flags=tskit.NODE_IS_SAMPLE)
    # We could do better here by mapping by node time, but we cannot assume that
    # node times are unique as discsim and ercs can have multiple nodes occuring
    # at the same time.
    for i in range(L):
        node_map = {j + 1: j for j in range(n)}
        for j in range(n + 1, len(pi[i])):
            node_map[j] = tables.nodes.add_row(time=tau[i][j], flags=0)
        for j in range(1, len(pi[i])):
            if pi[i][j] != 0:
                tables.edges.add_row(i, i + 1, node_map[pi[i][j]], node_map[j])
    tables.sort()
    tables.simplify()
    return tables.tree_sequence()
