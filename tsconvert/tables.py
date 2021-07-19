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
import pandas as pd


def pandas_to_table_metadata(table, dataframe, index_metadata_property):
    """
    Adds metadata to a table, in place, from a pandas dataframe. The dataframe index
    is matched to the metadata property specified. The table must have a metadata schema
    with a codec that supports dicts, e.g. JSON.

    :param AnyTableClass table: The table to modify
    :param pandas.dataframe dataframe: The dataframe to insert as metadata
    :param string index_metadata_property: The metadata property to match to the
        dataframe index.
    """

    row_for_id = {
        index: row for index, row in zip(dataframe.index, dataframe.to_dict("records"))
    }
    table_copy = table.copy()
    table.clear()
    for row in table_copy:
        table.append(
            row.replace(
                metadata={
                    **row.metadata,
                    **row_for_id.get(row.metadata[index_metadata_property], {}),
                }
            )
        )


def csv_to_node_metadata(ts, file_or_path, id_col, id_metadata_property, **kwargs):
    """
    Returns a new tree sequence with node metadata added from a csv/tsv.

    Rows in the tsv/csv are read and assigned as metadata on nodes, any property name
    collision with existing metadata will result in the existing property being
    overwritten. Rows are matched by the values of "id_col" (in the CSV) and
    "id_metadata_property" (in the existing metadata). The node table must have
    a metadata schema with a codec that supports dicts, e.g. JSON.

    :param TreeSequence ts: Tree sequence to add metadata to
    :param file_or_path file_or_path: CSV/TSV file to read, passed to pandas so can
        be an open file-like or a string path/URL.
    :param string id_col: Name of the column in the file to match to
        "id_metadata_property" in the existing metadata.
    :param string id_metadata_property: Name of the metadata property in the tree
        sequence to match to the CSV.
    """

    df = pd.read_csv(file_or_path, **kwargs)
    df = df.set_index(id_col, drop=False)
    tables = ts.dump_tables()
    pandas_to_table_metadata(tables.nodes, df, id_metadata_property)
    return tables.tree_sequence()
