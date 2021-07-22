#
# MIT License
#
# Copyright (c) 2021 Tskit Developers
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
import dendropy
import newick
import tskit


def from_ms(string):
    """
    Returns a tree sequence representation of the specified ms formatted tree output.
    """
    recording_trees = False
    tree_lines = ""
    for line in string.splitlines():
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
                f"Problem reading integer # of positions spanned in tree {i}"
                + " (in ms format this preceeds the tree, in square braces)"
            )
        if len(trees.taxon_namespace) != sum(1 for _ in tree.leaf_node_iter()):
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
                f"Tree {i}: cannot have two internal nodes with the same time"
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


def to_ms(ts):
    """
    Returns an ms-formatted version of the specified tree sequence.
    """
    output = ""
    for tree in ts.trees():
        span = tree.interval[1] - tree.interval[0]
        output += f"[{span}]{tree.newick()}\n"
    return output


def from_newick(string, min_edge_length=0):
    """
    Returns a tree sequence representation of the specified newick string.

    The tree sequence will contain a single tree, as specified by the newick. All
    leaf nodes will be marked as samples (``tskit.NODE_IS_SAMPLE``). Newick names and
    comments will be written to the node metadata.

    :param string string: Newick string
    :param float min_edge_length: Replace any edge length shorter than this value by this
        value. Unlike newick, tskit doesn't support zero or negative edge lengths, so
        setting this argument to a small value is necessary when importing trees with
        zero or negative lengths.
    """
    trees = newick.loads(string)
    if len(trees) > 1:
        raise ValueError("Only one tree can be imported from a newick string")
    if len(trees) == 0:
        raise ValueError("Newick string was empty")
    tree = trees[0]
    tables = tskit.TableCollection(1)
    nodes = tables.nodes
    nodes.metadata_schema = tskit.MetadataSchema(
        {
            "codec": "json",
            "type": "object",
            "properties": {
                "name": {
                    "type": ["string"],
                    "description": "Name from newick file",
                },
                "comment": {
                    "type": ["string"],
                    "description": "Comment from newick file",
                },
            },
        }
    )

    id_map = {}

    def get_or_add_node(newick_node, time):
        if newick_node not in id_map:
            flags = tskit.NODE_IS_SAMPLE if len(newick_node.descendants) == 0 else 0
            metadata = {}
            if newick_node.name:
                metadata["name"] = newick_node.name
            if newick_node.comment:
                metadata["comment"] = newick_node.comment
            id_map[newick_node] = tables.nodes.add_row(
                flags=flags, time=time, metadata=metadata
            )
        return id_map[newick_node]

    root = next(tree.walk())
    get_or_add_node(root, 0)
    for newick_node in tree.walk():
        node_id = id_map[newick_node]
        for child in newick_node.descendants:
            length = max(child.length, min_edge_length)
            if length <= 0:
                raise ValueError(
                    "tskit tree sequences cannot contain edges with lengths"
                    " <= 0. Set min_edge_length to force lengths to a"
                    " minimum size"
                )
            child_node_id = get_or_add_node(child, nodes[node_id].time - length)
            tables.edges.add_row(0, 1, node_id, child_node_id)
    # Rewrite node times to fit the tskit convention of zero at the youngest leaf
    nodes = tables.nodes.copy()
    youngest = min(tables.nodes.time)
    tables.nodes.clear()
    for node in nodes:
        tables.nodes.append(node.replace(time=node.time - youngest + root.length))
    tables.sort()
    return tables.tree_sequence()


def get_dendropy_node_id(node):
    """
    From a DendroPy Node object, returns the ID of that node.
    """
    if node._label is None:
        return int(node.taxon._label.split()[1])
    else:
        return int(node._label.split()[1])


def from_argon(string, sequence_length):
    """
    Returns a tree sequence representation of an ARGON .trees file.
    (Does not include mutations!)
    """
    string_by_lines = string.splitlines()
    tables = tskit.TableCollection(sequence_length)
    id_map = {}

    # Each line corresponds to a tree.
    for line in string_by_lines:
        argon_tree = line.split("\t")
        # unused tmrca = argon_tree[0]
        left = float(argon_tree[1]) - 1
        right = float(argon_tree[2])

        # Read in the tree.
        tree = dendropy.Tree.get(data=argon_tree[3], schema="newick")

        # Make the tables.
        for node in tree.ageorder_node_iter():
            time = node.age
            node_id_argon = get_dendropy_node_id(node)
            children = list(node.child_nodes())
            # The keys of id_map need to include the time of the node,
            # because, annoyingly, ARGON sometimes gives the same
            # id label to nodes in different times
            if (node_id_argon, time) not in id_map:
                flags = tskit.NODE_IS_SAMPLE if len(children) == 0 else 0
                # TODO derive information from the node and store it as JSON metadata.
                id_map[(node_id_argon, time)] = tables.nodes.add_row(
                    flags=flags, time=time
                )
            node_id = id_map[(node_id_argon, time)]
            for child in children:
                tables.edges.add_row(
                    left,
                    right,
                    node_id,
                    id_map[(get_dendropy_node_id(child), child.age)],
                )

    tables.sort()
    return tables.tree_sequence().simplify()


def to_newick(tree, precision=16):
    # TODO implement a clean python version here.
    return tree.newick(precision=precision)
