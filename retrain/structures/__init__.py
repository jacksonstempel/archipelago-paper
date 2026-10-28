"""The structure bank used by DPFunc: for every protein an AlphaFold DB model of its exact sequence, else an AlphaFold2
structure predicted with ColabFold. Stages, in order: afdb (requests, search, fetch), colabfold (msa, fold), assemble."""
