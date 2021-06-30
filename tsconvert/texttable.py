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
import csv

import tskit.util


def csv_to_node_metadata(ts, file_or_path, id_col, id_metadata_property):
    """
    Returns a new tree sequence with node metadata added from a csv/tsv.

    Rows in the tsv/csv are read and assigned as metadata on nodes, any property name
    collision existing metadata will result in the existing property being overwritten.
    Rows are matched by the values of "id_col" (in the CSV) and "id_metadata_property"
    (in the existing metadata) matching.

    :param TreeSequence ts: Tree sequence to add metadata to
    :param file_or_path file_or_path: CSV/TSV file to read
    :param string id_col: Name of the column in the file to match to
        "id_metadata_property" in the existing metadata.
    :param string id_metadata_property: Name of the metadata property in the tree
        sequence to match to the CSV.
    """
    csvfile, local_file = tskit.util.convert_file_like_to_open_file(file_or_path, "r")
    try:
        dialect = csv.Sniffer().sniff(csvfile.read(1048576))
        csvfile.seek(0)
        reader = csv.DictReader(csvfile, dialect=dialect)
        row_for_id = {row[id_col]: row for row in reader}
        tables = ts.dump_tables()
        nodes = tables.nodes.copy()
        tables.nodes.clear()
        for node in nodes:
            tables.nodes.append(
                node.replace(
                    metadata={
                        **node.metadata,
                        **row_for_id.get(node.metadata[id_metadata_property], {}),
                    }
                )
            )
        return tables.tree_sequence()
    finally:
        if local_file:
            csvfile.close()
