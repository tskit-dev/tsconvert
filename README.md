# tsconvert

Utilities for converting [tskit](https://github.com/tskit-dev/tskit) tree sequences to
and from other formats. 

`tsconvert` is in early development

## Roadmap for v0.1

 - Where conversion methods are widely used and do not require any dependencies,
   they will be implemented in `tskit`, for example writing VCF. `tsconvert` can then call
   the `tskit` method so that it has a complete list. Methods that require dependancies
   will be implemented here to avoid dependency bloat in `tskit`.


 - Provide either `to_`, `from_` (or both as appropriate) methods 
   for the following formats: 
 - [ ] vcf
 - [ ] ms
 - [ ] argon
 - [ ] newick
 - [ ] newick + csv/tsv metadata
 - [ ] nexus
 - [ ] sgkit
 

 - Where a format doesn't provide a full tree sequence it will
   have a `parse_` method: 
 - [ ] plink (pedigree)

 
 - Allow these methods to be used from the CLI, including streams.
