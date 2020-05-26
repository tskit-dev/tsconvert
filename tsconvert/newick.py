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
    Returns a tree sequence representation of the specified ms formatted tree output.
    """
    recording_trees = False
    tree_lines = ""
    for i, line in enumerate(string.splitlines()):
        line = line.strip()
        if line.startswith("["):
            recording_trees = True
            tree_lines += line
        else:
            if recording_trees:
                break  # stop after the first non-tree line

    if tree_lines == "":
        raise ValueError(
            "Malformed input: no lines starting with [."
            " Make sure you run ms with the -T and -r flags."
        )

    trees = dendropy.TreeList.get(
        data=tree_lines,
        schema="newick",
        extract_comment_metadata=True,
        rooting="force-rooted",
    )

    if len(trees) == 0:
        raise ValueError("No valid trees in ms file")

    spans = []
    for i, tree in enumerate(trees):
        try:
            spans.append(
                float(tree.comments[0])
            )  # the initial [X] is the first comment
        except (ValueError, IndexError):
            raise ValueError(
                "Problem reading integer # of positions spanned in tree {}".format(i)
                + " (in ms format this preceeds the tree, in square braces)"
            )
        if len(trees.taxon_namespace) != sum([1 for l in tree.leaf_node_iter()]):
            raise ValueError(
                "Tree {} does not have all {} expected tips".format(
                    i, len(trees.taxon_namespace)
                )
            )
        # below we might want to set is_force_max_age=True to allow ancient tips
        # and work around branch length precision errors in ms
        tree.calc_node_ages()
        node_ages = [n.age for n in tree.ageorder_node_iter(include_leaves=False)]
        if len(set(node_ages)) != len(node_ages):
            raise ValueError(
                "Tree {}: cannot have two internal nodes with the same time".format(i)
            )

    # NB: here we could check that the sequence_length == nsites, where nsites is given
    # in the ms_line, as the second number following the -r switch

    tables = tskit.TableCollection(sum(spans))
    # Get the samples from the first tree.
    for node in trees[0].leaf_node_iter():
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=node.age)

    age_id_map = {}
    left = 0
    for span, tree in zip(spans, trees):
        right = left + span
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


def from_beast(string, precision=6):
    """
    Reads in a BEAST nexus description.
    """
    trees = dendropy.TreeList.get(data=string, schema="nexus")

    print("Got", len(trees))
    positions = []
    for i, tree in enumerate(trees):
        label_toks = tree.label.split()
        assert label_toks[0] == "STATE"
        positions.append(int(label_toks[1]))
        # TODO likely some parameter tweaks would be needed to make this work reliably
        tree.calc_node_ages()
        # FIXME assuming for now what all samples are contemporary.
        node_ages = [
            round(n.age, precision)
            for n in tree.ageorder_node_iter(include_leaves=False)
        ]
        if len(set(node_ages)) != len(node_ages):
            raise ValueError(
                "Tree {}: cannot have two internal nodes with the same time".format(i)
            )

    assert len(positions) > 1
    assert positions[0] == 0
    diff = positions[-1] - positions[-2]
    positions.append(positions[-1] + diff)

    tables = tskit.TableCollection(positions[-1])
    sample_id_map = {}
    # Get the samples from the first tree.
    for node in trees[0].leaf_node_iter():
        node_id = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=node.age)
        sample_id_map[node.taxon.label] = node_id

    age_id_map = {}
    left = 0
    for right, tree in zip(positions[1:], trees):
        for node in tree.ageorder_node_iter(include_leaves=False):
            age = round(node.age, precision)
            # print("node", age, node.annotations)
            # print(node.description())
            children = list(node.child_nodes())
            if node.age not in age_id_map:
                age_id_map[age] = tables.nodes.add_row(flags=0, time=age)
            parent_id = age_id_map[age]
            for child in children:
                if child.is_leaf():
                    child_id = sample_id_map[child.taxon.label]
                else:
                    child_age = round(child.age, precision)
                    child_id = age_id_map[child_age]
                tables.edges.add_row(left, right, parent_id, child_id)
        left = right

    # import numpy as np
    # print(tables.nodes)
    # t = tables.nodes.time
    # np.sort(t)
    # print(t)
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
        span = tree.interval[1] - tree.interval[0]
        output += "[{}]{}\n".format(span, tree.newick())
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


def to_newick(tree, precision=16):
    # TODO implement a clean python version here.
    return tree.newick(precision=precision)
