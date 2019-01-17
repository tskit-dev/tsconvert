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
import tskit

import dendropy

def from_newick(string):
    """
    Returns a tree sequence representation of the specified newick string.
    """
    tree = dendropy.Tree.get(data=string, schema="newick")
    tables = tskit.TableCollection(1)

    id_map = {}
    for node in tree.ageorder_node_iter():
        children = list(node.child_nodes())
        if node not in id_map:
            flags = tskit.NODE_IS_SAMPLE if len(children) == 0 else 0
            # TODO derive information from the node and store it as JSON metadata.
            id_map[node] = tables.nodes.add_row(flags=flags, time=node.age)
        node_id = id_map[node]
        for child in children:
            tables.edges.add_row(0, 1, node_id, id_map[child])
    return tables.tree_sequence()


def to_newick(tree):
    # TODO implement a clean python version here.
    return tree.newick()
