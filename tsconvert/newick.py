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


def from_ms(string):
    """
    Returns a tree sequence represetation of the specified ms formatted tree output.
    """
    lines = [line.strip() for line in string.splitlines()]

    # Skip any leading non-tree lines.
    j = 0
    while j < len(lines) and not lines[j].startswith("["):
        j += 1
    if j == len(lines):
        raise ValueError("Malformed input: no lines starting with [")

    lines = lines[j:]
    sequence_length = 0
    trees = []
    for i, line in enumerate(lines):
        if len(line) == 0:
            break
        if not line.startswith("["):
            raise ValueError("Line {} not in ms format: missing [".format(i+j+1))
        index = line.index("]", 1)
        length = float(line[1: index])
        sequence_length += length
        tree = dendropy.Tree.get(data=line[index + 1:], schema="newick")
        node_ages = [node.age for node in tree.ageorder_node_iter(include_leaves=False)]
        if len(set(node_ages)) != len(node_ages):
            raise ValueError(
                "Line {}: cannot have two internal nodes with the same time"
                .format(i+j+1))
        trees.append((length, tree))

    tables = tskit.TableCollection(sequence_length)
    # Get the samples from the first tree
    for node in trees[0][1].leaf_node_iter():
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=node.age)

    age_id_map = {}
    left = 0
    for length, tree in trees:
        right = left + length
        for node in tree.ageorder_node_iter(include_leaves=False):
            children = list(node.child_nodes())
            if node.age not in age_id_map:
                age_id_map[node.age] = tables.nodes.add_row(flags=0, time=node.age)
            parent_id = age_id_map[node.age]
            for child in children:
                if child.is_leaf():
                    child_id = int(child.taxon.label) - 1
                else:
                    child_id = age_id_map[child.age]
                tables.edges.add_row(left, right, parent_id, child_id)
        left = right
    tables.sort()
    # Simplify will squash together any edges, removing redundancy.
    tables.simplify()
    return tables.tree_sequence()


def to_ms(ts):
    """
    Returns an ms-formatted version of the specified tree sequence.
    """
    output = ""
    for tree in ts.trees():
        length = tree.interval[1] - tree.interval[0]
        output += "[{}]{}\n".format(length, tree.newick())
    return output


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
