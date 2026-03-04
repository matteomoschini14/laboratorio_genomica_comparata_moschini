## Reference-based scaffolding

This stage of the assembly process involves scaffolding. A scaffold is a computational representation of contigs, reordered and oriented according to specific criteria to best reflect the original chromosomal structure. For this purpose, we employed RagTag, a suite of software tools designed to improve genome assemblies through homology-based correction, scaffolding, and merging.

### Reference-based error correction

 ```bash
ragtag.py correct -t 20 <REFERENCE_GENOME> <DRAFT_GENOME>
 ```

### Scaffolding 

```bash
ragtag.py scaffold -C -t 20 -o <OUTPUT_DIR> <REFERENCE_GENOME> <CORRECTED_DRAFTGENOME>
```

>The scaffolding process is highly time-consuming; therefore, due to logistical constraints, a link to a pre-computed result has been provided in this instance as well.
