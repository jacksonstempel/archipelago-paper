# Protein function review — batch 3 of 16

25 proteins; 25 assertions. Name your output `AGREEMENT_RESULTS_03.json`.

Review the supplied protein–GO molecular-function assertions. Research the evidence for each assertion and decide whether it is sufficiently supported to accept.

Work protein by protein, reusing relevant research across its assertions. The packet contains sequences, identifiers, GO definitions, and existing annotation provenance. Use primary literature and authoritative databases as appropriate. Investigate the underlying evidence rather than treating an existing annotation as confirmation.

## Decision rubric

Return **accept**, **do not accept**, or, exceptionally, **needs adjudication**.

Accept when either:

- **Experimental evidence** establishes the requested function for the protein. Biochemical, structural, interaction, and genetic experiments can qualify; a purified-protein assay is not required. The evidence must establish the particular molecular function, rather than merely a related phenotype or biological process.
- **A well-supported inference** establishes the function. Explain why the experimental evidence from another protein, or a curated functional rule, applies to the supplied protein. Relevant domain conservation, essential residues, orthology, and known functional differences should inform this judgment. Generic family resemblance or a repeated annotation alone is insufficient.

Apply these principles consistently:

- **Match the exact GO definition.** Evidence for a broader or related activity does not establish a more specific claim.
- **Weak activity is not automatically invalid.** Accept convincing, specific activity even if its physiological role is unknown. Do not accept a marginal signal that remains indistinguishable from nonspecific association or experimental artifact.
- **Complex-dependent activity can qualify.** Accept an established functional contribution to a complex; mere complex membership does not suffice. Record whether the evidence supports `enables` or `contributes_to`.
- **Relevant constructs and mature products can qualify.** Experiments on a functional domain or naturally processed product may establish a function of the supplied sequence. Strain, isoform, or construct differences should prevent acceptance only when they materially weaken the attribution.
- **Check the inference, not every historical pathway.** One adequately substantiated route to support is sufficient; you do not need to reconstruct every annotation’s provenance.
- **Failure to find support is not evidence of falsehood.** Return “do not accept” when sufficient support cannot be established. Reserve claims of contradiction for evidence genuinely opposing the assertion.

Be willing to accept justified functional inferences without demanding a new direct assay for every protein. Conversely, before rejecting an assertion, consider whether relevant evidence appears under another identifier, in a functional domain, or in a different study from the one originally cited.

Use **needs adjudication** only when a concrete conflict or ambiguity prevents a defensible decision under this rubric—for example, incompatible authoritative term definitions or strong conflicting evidence. Ordinary missing evidence should result in **do not accept**, not adjudication.

## Deliverable

Return a JSON file containing one entry per supplied assertion:

- `assertion_id`: preserve exactly.
- `decision`: `accept`, `do_not_accept`, or `needs_adjudication`.
- `support_type`: `experimental`, `inference`, or `none_established`.
- `supported_relation`: `enables`, `contributes_to`, or null.
- `rationale`: concise explanation connecting the evidence to the exact GO term.
- `sources`: URLs, relevant figure/section/record, what was established, and the studied protein or organism where relevant. State when only an abstract or limited material was accessible.
- `limitations`: any qualification material to the decision.
- `adjudication_question`: null unless adjudication is requested; explain precisely what remains undecidable.

Distinguish facts you verified from assumptions or inferences. Do not claim sequence comparisons, source access, or measurements you did not perform. Note evidence published after the packet’s annotation snapshot.

Review every supplied assertion. If you cannot finish, list unreviewed assertion IDs separately rather than assigning them an evidence verdict. Include a short summary of completion and decision counts.


## Protein and assertion data

```json
{
 "batch": 3,
 "source_release": "GOA 2026-07-29",
 "proteins": [
  {
   "sequence_sha256": "9ba5d2bcc7d021a803b02d1629ed3dfa05d68d14f4fa6bbc80346a09df0acd78",
   "sequence": "MAKYTIVDKDTCIACGACGAAAPDIYDYDDEGIAFVTLDENKGVVEVPEVLEEDMIDAFEGCPTDSIKVADEPFEGDPLKFE",
   "taxon_ids": [
    "224308"
   ],
   "accessions": [
    "A0A6H0WKK3",
    "A0A9Q4HNS4",
    "A0A9W5LEU7",
    "A0AAP3CKJ3",
    "A0AAP3FUX2",
    "A0ABT9DLF2",
    "A0ABU0V378",
    "A0ACD3ZU38",
    "G4NW02",
    "P50727"
   ],
   "assertions": [
    {
     "assertion_id": "9ba5d2bcc7d021a803b02d1629ed3dfa05d68d14f4fa6bbc80346a09df0acd78:GO:0005506",
     "go_id": "GO:0005506",
     "go_name": "iron ion binding",
     "go_definition": "\"Binding to an iron (Fe) ion.\" [GOC:ai]",
     "supporting_assertions": [
      {
       "accession": "A0ABU0V378",
       "symbol": "RAQ16_02540",
       "relation": "enables",
       "go_id": "GO:0005506",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001080"
       ],
       "object_name": "Ferredoxin",
       "taxon": [
        "taxon:2054641"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001080/"
       ]
      },
      {
       "accession": "A0ABU0V378",
       "symbol": "RAQ16_02540",
       "relation": "enables",
       "go_id": "GO:0005506",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR001348782"
       ],
       "object_name": "Ferredoxin",
       "taxon": [
        "taxon:2054641"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000104.html"
       ]
      },
      {
       "accession": "A0A6H0WKK3",
       "symbol": "G4P54_12070",
       "relation": "enables",
       "go_id": "GO:0005506",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001080"
       ],
       "object_name": "Ferredoxin",
       "taxon": [
        "taxon:227866"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001080/"
       ]
      },
      {
       "accession": "A0A6H0WKK3",
       "symbol": "G4P54_12070",
       "relation": "enables",
       "go_id": "GO:0005506",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR001348782"
       ],
       "object_name": "Ferredoxin",
       "taxon": [
        "taxon:227866"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000104.html"
       ]
      },
      {
       "accession": "G4NW02",
       "symbol": "GYO_2536",
       "relation": "enables",
       "go_id": "GO:0005506",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001080"
       ],
       "object_name": "Ferredoxin",
       "taxon": [
        "taxon:1052585"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001080/"
       ]
      },
      {
       "accession": "G4NW02",
       "symbol": "GYO_2536",
       "relation": "enables",
       "go_id": "GO:0005506",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR001348782"
       ],
       "object_name": "Ferredoxin",
       "taxon": [
        "taxon:1052585"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000104.html"
       ]
      },
      {
       "accession": "A0ABT9DLF2",
       "symbol": "KHP33_011835",
       "relation": "enables",
       "go_id": "GO:0005506",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001080"
       ],
       "object_name": "Ferredoxin",
       "taxon": [
        "taxon:2944916"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001080/"
       ]
      },
      {
       "accession": "A0ABT9DLF2",
       "symbol": "KHP33_011835",
       "relation": "enables",
       "go_id": "GO:0005506",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR001348782"
       ],
       "object_name": "Ferredoxin",
       "taxon": [
        "taxon:2944916"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000104.html"
       ]
      },
      {
       "accession": "A0ACD3ZU38",
       "symbol": "M0696_11390",
       "relation": "enables",
       "go_id": "GO:0005506",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001080"
       ],
       "object_name": "Ferredoxin",
       "taxon": [
        "taxon:2715209"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001080/"
       ]
      },
      {
       "accession": "A0ACD3ZU38",
       "symbol": "M0696_11390",
       "relation": "enables",
       "go_id": "GO:0005506",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR001348782"
       ],
       "object_name": "Ferredoxin",
       "taxon": [
        "taxon:2715209"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000104.html"
       ]
      },
      {
       "accession": "P50727",
       "symbol": "fer",
       "relation": "enables",
       "go_id": "GO:0005506",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001080"
       ],
       "object_name": "Ferredoxin",
       "taxon": [
        "taxon:224308"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001080/"
       ]
      },
      {
       "accession": "A0AAP3CKJ3",
       "symbol": "MKF32_11415",
       "relation": "enables",
       "go_id": "GO:0005506",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001080"
       ],
       "object_name": "Ferredoxin",
       "taxon": [
        "taxon:72361"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001080/"
       ]
      },
      {
       "accession": "A0AAP3CKJ3",
       "symbol": "MKF32_11415",
       "relation": "enables",
       "go_id": "GO:0005506",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR001348782"
       ],
       "object_name": "Ferredoxin",
       "taxon": [
        "taxon:72361"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000104.html"
       ]
      },
      {
       "accession": "A0A9W5LEU7",
       "symbol": "BSI_41720",
       "relation": "enables",
       "go_id": "GO:0005506",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001080"
       ],
       "object_name": "Ferredoxin",
       "taxon": [
        "taxon:1236548"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001080/"
       ]
      },
      {
       "accession": "A0A9W5LEU7",
       "symbol": "BSI_41720",
       "relation": "enables",
       "go_id": "GO:0005506",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR001348782"
       ],
       "object_name": "Ferredoxin",
       "taxon": [
        "taxon:1236548"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000104.html"
       ]
      },
      {
       "accession": "A0AAP3FUX2",
       "symbol": "fer",
       "relation": "enables",
       "go_id": "GO:0005506",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001080"
       ],
       "object_name": "Ferredoxin",
       "taxon": [
        "taxon:72360"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001080/"
       ]
      },
      {
       "accession": "A0AAP3FUX2",
       "symbol": "fer",
       "relation": "enables",
       "go_id": "GO:0005506",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR001348782"
       ],
       "object_name": "Ferredoxin",
       "taxon": [
        "taxon:72360"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000104.html"
       ]
      },
      {
       "accession": "A0A9Q4HNS4",
       "symbol": "O0R52_11805",
       "relation": "enables",
       "go_id": "GO:0005506",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001080"
       ],
       "object_name": "Ferredoxin",
       "taxon": [
        "taxon:260554"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001080/"
       ]
      },
      {
       "accession": "A0A9Q4HNS4",
       "symbol": "O0R52_11805",
       "relation": "enables",
       "go_id": "GO:0005506",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR001348782"
       ],
       "object_name": "Ferredoxin",
       "taxon": [
        "taxon:260554"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000104.html"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "0aa976750875ca9436456d55fa6a2d4c5da71a9d509b56735e213f87118c1a69",
   "sequence": "MGLNLNPILRQELANLDKDTESRKSAMKALKSYVKDLDSKAIPGFLAQVFETKETNSLSGEYTISLYEILARVHGPNIVPQIDTIMSTIVKTLASSAGSFPLQQACSKVIPAIARYGIDPTTTEDKKRVIIHSLCKPLTDSLLASQESLTSGAALCLKALVDSDNWRFASDEMVNRVCQNVVVALDSNSNQTHLQMGLVMSLAKHNPLIVEAYARLLIHTGLRILGFGVSEGNSQKRLSAVQMLNFLMKCLDPRSIYSEVELIIKEMERCQSDQMAYVRGAAYEAMMTSKRIAAELESKMEKGCRSVTGSNFSRRNCSSIVPDYSLSPESQTLGSFSGYDSPVESSPISHTSCNSEFDRRSVNRKLWRRDENGGVVDISLKDGLFSRVTKGSTTVSDSPLVPYDTCENGDEFEGFLMESLRNTTPSPQRQRSRRINAEDFNIFSTPRKLISSLQYPDDVDLDHSDIQSPILRGEREKTIGSRKNPKLRKQFPTMVETMSSTITVSEDTAQTQMITGKKKKKKMSYAKLVIAISFVVVALFATVILMVNQDDDVGYYTVPT",
   "taxon_ids": [
    "3702"
   ],
   "accessions": [
    "A0A8T2GR76",
    "A0A8T2HEX6",
    "Q5XVI1"
   ],
   "assertions": [
    {
     "assertion_id": "0aa976750875ca9436456d55fa6a2d4c5da71a9d509b56735e213f87118c1a69:GO:0008017",
     "go_id": "GO:0008017",
     "go_name": "microtubule binding",
     "go_definition": "\"Binding to a microtubule, a filament composed of tubulin monomers.\" [GOC:krc]",
     "supporting_assertions": [
      {
       "accession": "A0A8T2HEX6",
       "symbol": "ISN44_As01g045440",
       "relation": "enables",
       "go_id": "GO:0008017",
       "references": [
        "GO_REF:0000118"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "PANTHER:PTN002526378"
       ],
       "object_name": "TORTIFOLIA1/SINE1-2 N-terminal domain-containing protein",
       "taxon": [
        "taxon:45249"
       ],
       "annotation_date": "20260727",
       "assigned_by": "TreeGrafter",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000118.html"
       ]
      },
      {
       "accession": "A0A8T2HEX6",
       "symbol": "ISN44_As01g045440",
       "relation": "enables",
       "go_id": "GO:0008017",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR033337"
       ],
       "object_name": "TORTIFOLIA1/SINE1-2 N-terminal domain-containing protein",
       "taxon": [
        "taxon:45249"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR033337/"
       ]
      },
      {
       "accession": "A0A8T2GR76",
       "symbol": "ISN45_At01g046460",
       "relation": "enables",
       "go_id": "GO:0008017",
       "references": [
        "GO_REF:0000118"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "PANTHER:PTN002526378"
       ],
       "object_name": "TORTIFOLIA1/SINE1-2 N-terminal domain-containing protein",
       "taxon": [
        "taxon:1240361"
       ],
       "annotation_date": "20260727",
       "assigned_by": "TreeGrafter",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000118.html"
       ]
      },
      {
       "accession": "A0A8T2GR76",
       "symbol": "ISN45_At01g046460",
       "relation": "enables",
       "go_id": "GO:0008017",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR033337"
       ],
       "object_name": "TORTIFOLIA1/SINE1-2 N-terminal domain-containing protein",
       "taxon": [
        "taxon:1240361"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR033337/"
       ]
      },
      {
       "accession": "Q5XVI1",
       "symbol": "SINE1",
       "relation": "enables",
       "go_id": "GO:0008017",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR033337"
       ],
       "object_name": "Protein SINE1",
       "taxon": [
        "taxon:3702"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR033337/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "b17cb8007c43ef1c1c5769d4542e0af31ef615feb9fee0af3940522f7f7246b4",
   "sequence": "MPLHVSLANGNRDLDYDSVQPYFMCDDEEEDVHHQQPPQPPAPSEDIWKKFELLPTPRPSPGHAGLYSPPCEAVAVSFAPRDHDGDSFSIADLPELPGGDAVKQSFVCDPDDETFVKNIILQDCMWNGFSASAKLVSKLDPYQAVRKEGTGVSLAADVEPATPPDCTCNT",
   "taxon_ids": [
    "10090"
   ],
   "accessions": [
    "A0A8C6IEM5",
    "A0PJL6",
    "Q6P8Z1"
   ],
   "assertions": [
    {
     "assertion_id": "b17cb8007c43ef1c1c5769d4542e0af31ef615feb9fee0af3940522f7f7246b4:GO:0003700",
     "go_id": "GO:0003700",
     "go_name": "DNA-binding transcription factor activity",
     "go_definition": "\"A transcription regulator activity that modulates transcription of gene sets via selective and non-covalent binding to a specific double-stranded genomic DNA sequence (sometimes referred to as a motif) within a cis-regulatory region. Regulatory regions include promoters (proximal and distal) and enhancers. Genes are transcriptional units, and include bacterial operons.\" [GOC:txnOH-2018]",
     "supporting_assertions": [
      {
       "accession": "A0A8C6IEM5",
       "symbol": "A0A8C6IEM5",
       "relation": "enables",
       "go_id": "GO:0003700",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR002418",
        "InterPro:IPR012682"
       ],
       "object_name": "Transcription regulator Myc N-terminal domain-containing protein",
       "taxon": [
        "taxon:10103"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR002418/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR012682/"
       ]
      },
      {
       "accession": "A0PJL6",
       "symbol": "Bmyc",
       "relation": "enables",
       "go_id": "GO:0003700",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR002418",
        "InterPro:IPR012682"
       ],
       "object_name": "Transcription regulator Myc N-terminal domain-containing protein",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR002418/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR012682/"
       ]
      },
      {
       "accession": "Q6P8Z1",
       "symbol": "Mycb",
       "relation": "enables",
       "go_id": "GO:0003700",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR002418",
        "InterPro:IPR012682"
       ],
       "object_name": "Protein B-Myc",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR002418/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR012682/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "53f2211b9356924319046629198547d1b69bc0164d3a9929e2ef7e322db4a1d7",
   "sequence": "MAAVGRVGSFGSSPPGLASTYASGPLANELASGSGGPAAGDDEDGQNLWSRILREVSTRSRSKLPTGKNVLLLGEDGAGKTSLIRRIQGIEEYKKGRGLEYLYLNVHDEDRDDQTRCNVWILDGDLYHKGLLKFSLDALSLRDTLVMLVVDMSKPWTALDSLQKWASVVREHVDKLKIPPEEMKEMEQKLIRDFQEYVEPGEDFPASPQRRATAAQEDRDDSVVLPLGADTLTHNLGLPVLVVCTKCDAISVLEKEHDYRDEHFDFIQSHIRKFCLQYGAALIYTSVKENKNIDLVYKYIVQKLYGFPYKIPAVVVEKDAVFIPAGWDNDKKIGILHENFQTLKIEDNFEDIITKPPVRKFVHEKEIMAEDDQVFLMKLQSLLAKQPPTAAGRPVDASPRVPGGSPRTPNRSVSSNVASVSPIPAGSKKIDPNMKAGATSEGVLANFFNSLLSKKTGSPGGPGVGGSPGGGAAGASTSLPPSAKKSGQKPVLSDVHAELDRITRKPASVSPTTPPSPTEGEAS",
   "taxon_ids": [
    "10116"
   ],
   "accessions": [
    "Q9QXU8"
   ],
   "assertions": [
    {
     "assertion_id": "53f2211b9356924319046629198547d1b69bc0164d3a9929e2ef7e322db4a1d7:GO:0019003",
     "go_id": "GO:0019003",
     "go_name": "GDP binding",
     "go_definition": "\"Binding to GDP, guanosine 5'-diphosphate.\" [GOC:ai]",
     "supporting_assertions": [
      {
       "accession": "Q9QXU8",
       "symbol": "Dync1li1",
       "relation": "enables",
       "go_id": "GO:0019003",
       "references": [
        "GO_REF:0000121"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "UniProtKB:Q9Y6G9"
       ],
       "object_name": "Cytoplasmic dynein 1 light intermediate chain 1",
       "taxon": [
        "taxon:10116"
       ],
       "annotation_date": "20160121",
       "assigned_by": "RGD",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000121.html",
        "https://www.uniprot.org/uniprotkb/Q9Y6G9/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "10aca9a2d5ea0395fb4ac7d3ff219b906d0e0ee00dcb829d1a596ff74e69a79e",
   "sequence": "MQKEMKMIKDEDVHFDLAVKKTPSFPHCLQPVASRGKAPQRHPFPEALRGPFSQFRYEPPPGDLDGFPGVFEGAGSRKRKSMPTKMPYNHPAEEVTLALHSEENKNHGLPNLPLLFPQPPRPKYDSQMIDLCNVGFQFYRSLEHFGGKPVKQEPIKPSAVWPQPTPTPFLPTPYPYYPKVHPGLMFPFFVPSSSPFPFSRHTFLPKQPPEPLLPRKAEPQESEETKQKVERVDVNVQIDDSYYVDVGGSQKRWQCPTCEKSYTSKYNLVTHILGHSGIKPHACTHCGKLFKQLSHLHTHMLTHQGTRPHKCQVCHKAFTQTSHLKRHMMQHSEVKPHNCRVCGRGFAYPSELKAHEAKHASGRENICVECGLDFPTLAQLKRHLTTHRGPIQYNCSECDKTFQYPSQLQNHMMKHKDIRPYICSECGMEFVQPHHLKQHSLTHKGVKEHKCGICGREFTLLANMKRHVLIHTNIRAYQCHLCYKSFVQKQTLKAHMIVHSDVKPFKCKLCGKEFNRMHNLMGHMHLHSDSKPFKCLYCPSKFTLKGNLTRHMKVKHGVMERGLHSQGLGRGRIALAQTAGVLRSLEQEEPFDLSQKRRAKVPVFQSDGESAQGSHCHEEEEEDNCYEVEPYSPGLAPQSQQLCTPEDLSTKSEHAPEVLEEACKEEKEDASKGEWEKRSKGDLGAEGGQERDCAGRDECLSLRAFQSTRRGPSFSDYLYFKHRDESLKELLERKMEKQAVLLGI",
   "taxon_ids": [
    "9606"
   ],
   "accessions": [
    "Q8N895"
   ],
   "assertions": [
    {
     "assertion_id": "10aca9a2d5ea0395fb4ac7d3ff219b906d0e0ee00dcb829d1a596ff74e69a79e:GO:0003700",
     "go_id": "GO:0003700",
     "go_name": "DNA-binding transcription factor activity",
     "go_definition": "\"A transcription regulator activity that modulates transcription of gene sets via selective and non-covalent binding to a specific double-stranded genomic DNA sequence (sometimes referred to as a motif) within a cis-regulatory region. Regulatory regions include promoters (proximal and distal) and enhancers. Genes are transcriptional units, and include bacterial operons.\" [GOC:txnOH-2018]",
     "supporting_assertions": [
      {
       "accession": "Q8N895",
       "symbol": "ZNF366",
       "relation": "enables",
       "go_id": "GO:0003700",
       "references": [
        "GO_REF:0000033"
       ],
       "evidence_code": "IBA",
       "with_from": [
        "CGD:CAL0000176618",
        "CGD:CAL0000190493",
        "FB:FBgn0001325",
        "FB:FBgn0002573",
        "FB:FBgn0030012",
        "FB:FBgn0037446",
        "MGI:MGI:1920441",
        "MGI:MGI:2682609",
        "PANTHER:PTN002811104",
        "SGD:S000005639",
        "UniProtKB:P56270",
        "UniProtKB:Q66K89",
        "UniProtKB:Q9BU19",
        "UniProtKB:Q9UDV6",
        "WB:WBGene00003909"
       ],
       "object_name": "Zinc finger protein 366",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20250904",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000033.html",
        "https://www.uniprot.org/uniprotkb/P56270/entry",
        "https://www.uniprot.org/uniprotkb/Q66K89/entry",
        "https://www.uniprot.org/uniprotkb/Q9BU19/entry",
        "https://www.uniprot.org/uniprotkb/Q9UDV6/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "cb757784654e84d6045e2383dd0f39aa308801f69d39569eca3ef4d86f3ad789",
   "sequence": "MQTPVNIPVPVLRLPRGPDGFSRGFAPDGRRAPLRPEVPEIQECPIAQESLESQEQRARAALRERYLRSLLAMVGHQVSFTLHEGVRVAAHFGATDLDVANFYVSQLQTPIGVQAEALLRCSDIISYTFKP",
   "taxon_ids": [
    "9606"
   ],
   "accessions": [
    "Q9H840"
   ],
   "assertions": [
    {
     "assertion_id": "cb757784654e84d6045e2383dd0f39aa308801f69d39569eca3ef4d86f3ad789:GO:0003723",
     "go_id": "GO:0003723",
     "go_name": "RNA binding",
     "go_definition": "\"Binding to an RNA molecule or a portion thereof.\" [GOC:jl, GOC:mah]",
     "supporting_assertions": [
      {
       "accession": "Q9H840",
       "symbol": "GEMIN7",
       "relation": "enables",
       "go_id": "GO:0003723",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR047575"
       ],
       "object_name": "Gem-associated protein 7",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR047575/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "51268a552f903a9e5d04d911b075bf7508ddab4c20b16bc64236176f74085e6b",
   "sequence": "MDSILSKQLVDKTGFVRVPIKHFDCTMLTLALPTFDVSKMVDRITIDFNLDDIQGASEIGSTLLPSMSIDVEDMANFVHDFTFGHLADKTDRLLMREFPMMNDGFDHLSPDMIIKTTSGVYNIVEFTNFRGDERGAFQAAMIKLAKYEVPCENRSQGRTVVLYVVSAYRAWCMVYLELERTLKQREMVYRYRLALSVMDELRTLFPELSSTDEELGKTERELPAMVSSIQINWSVTESVFPPFSREMFDRFRSSPPDSEYITRIVSRCLINSQEKLINSSFFAEGNDKALRFSKNAEECSLAVERALNQYRAEDNLRDLNDHKSTIQLPPWLSYHDVDGKDLCPLQGLDVRGDHPMCNLWREVVTSANLEEIERMHDDAAAELEFAFGSKGQARERNRYHRVHLNMGSDVLVYIAALGVNGKKHKADTLVQQMRDRSKQPFSPDHDVITYLNFSLHALVTCGQQMRTCTALSLVIRDSVGSPEDSSAILVRKGFHEIITEHYKFMGSRIGHGCQMVSLIGAELSASVKQHVKPNYFVIKRLLGSGIFLLIKPTSSKSHIFVSFALSALAGPLISPLPGFSSPTKMLGILLVTDFVSYKLSKLTNLCKCVSLMESSFSFWAEAFGIQAGTLVGDFVPRSSDSAAMDASYMGKLSLLTLLEDKAATEELQTIARYIIMEGFVSPPEIPKPHKMTSKFPKVLRSELQVYLLNCLCRTIQRIAGEPFILKKKDGSISWGGMFNPFSGRPLLDMQPLISCCYNGYFKNKEEETEPSSLSGMYKKIIELEHLRPQSDAFLGYKDPELPRMHEFSVSYLKEACNHAKLVLRSLYGQNFMEQIDNQIIRELSGLTLERLATLKATSNFNENWYVYKDVADKNYTRDKLLVKMSKYASEGKSLAIQKFEDCMRQIESQGCMHICLFKKQQHGGLREIYVMGAEERIVQSVVETIARSIGKFFASDTLCNPPNKVKIPETHGIRARKQCKGPVWTCATSDDARKWNQGHFVTKFALMLCEFTSPKWWPLIIRGCSMFTKKRMMMNLNYLKILDGHRELDIRDDFVMDLFKAYHGEAEVPWAFKGKTYLETTTGMMQGILHYTSSLLHTIHQEYIRSLSFKIFNLKVAPEMSKSLVCDMMQGSDDSSMLISFPADDEKVLTRCKVAAAICFRMKKELGVYLAIYPSEKSTANTDFVMEYNSEFYFHTQHVRPTIRWIAACCSLPEVETLVARQEEASNLMTSVTEGGGSFSLAAIIQQAQCTLHYMLMGMGVSELFLEYKKAVLKWNDPGLGFFLLDNPYACGLGGFRFNLFKAITRTDLQKLYAFFMKKVKGSAARDWADEDVTIPETCSVSPGGALILSSSLKWGSRKKFQKLRDRLNIPENWIELINENPEVLYRAPRTGPEILLRIAEKVHSPGVVSSLSSGNAVCKVMASAVYFLSATIFEDTGRPEFNFLEDSKYSLLQKMAAYSGFHGFNDMEPEDILFLFPNIEELESLDSIVYNKGEIDIIPRVNIRDATQTRVTIFNEQKNLRTSPEKLVSDKWFGTQKSRIGKTTFLAEWEKLKKIVKWLEDAPEATLAHTPLNNHIQVRNFFARMESKPRTVRITGAPVKKRSGVSKIAMVIRDHFSRMGHLRGVEDLAGFTRSVSAEILKHFLFCILQGPYSESYKLQLIYRVLSSVSNVEIKESDGKTKTNLIGILQRFLDGDHVVPIIEEMGAGTVGGFIKRQQSKVVQNKVVYYGVGIWRGFMDGYQVHLEIENDIGQPPRLRNVTTNCQSSPWDLSIPIRQWAEDMGVTNNQDYSSKSSRGARYWMHSFRMQGPSKPFGCPVYIIKGDMSDVIRLRKEEVEMKVRGSTLNLYTKHHSHQDLHILSYTASDNDLSPGIFKSISDEGVAQALQLFEREPSNCWVRCESVAPKFISAILEICEGKRQIKGINRTRLSEIVRICSESSLRSKVGSMFSFVANVEEAHDVDYDALMDLMIEDAKNNAFSHVVDCIELDVSGPYEMESFHGRSTLTCTPSTILIRTYTFYLTLHQTMISVQAFQVILDEGVLLIALVNNYLRGSKANCWVRCESVAPKFISAILEICEGKRQIKGINRTRLSEIVEFVLNLPKIKSRIYVLICRQCHGANFPPISVRRLMLEDIASVARRLIIVASFGS",
   "taxon_ids": [
    "11589"
   ],
   "accessions": [
    "P27316"
   ],
   "assertions": [
    {
     "assertion_id": "51268a552f903a9e5d04d911b075bf7508ddab4c20b16bc64236176f74085e6b:GO:0046872",
     "go_id": "GO:0046872",
     "go_name": "metal ion binding",
     "go_definition": "\"Binding to a metal ion.\" [GOC:ai]",
     "supporting_assertions": [
      {
       "accession": "P27316",
       "symbol": "L",
       "relation": "enables",
       "go_id": "GO:0046872",
       "references": [
        "GO_REF:0000043"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB-KW:KW-0479"
       ],
       "object_name": "RNA-directed RNA polymerase L",
       "taxon": [
        "taxon:11589"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000043.html"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "0f57b974e0c77533efdb6aba44687ebdd5964014cd4ba6be30332305fa60970a",
   "sequence": "MKAAEIKRDLTNIQKSMSEINDLAKERITGGPGSISTTSASAITAPSTMSQTTTSRLAPKLTSAHPSIDDLRGLSRQDKITQLQKKIRASFENLVDHDDSNVIVTLPDDDDCPHNHFGSGLNLTHPTAAQLSASGLSGSSKTIDTIKFQEKSMKTESKTKVVTDGFSSEQATSNSAEMKRLQAGDIDYQESKGASAMRNRLEVDGVKTEENAAVIKEALSLRTGDITQQASNNVAASSITVQSENFSADKKAISQSQQSQTMTSNGIISQEKHVSSASQANYSMSHKGVSSTGSSMITSSSQMSAMNGQMLKLADLKLDDLKSLTAGSGQQEIEQTINKYSNMLTSIVSSLQEDERGGSAITVHDVGGKKSQYLEKINEVIRRAWAVPTHGHELGYSLCNSLRQSGGLDLLMKNCVKPDLQFSSAQLLEQCLTTENRKHVVDNGLDKVVNVACVCTKNSNMEHSRVGTGILEHLFKHSEGTCSDVIRLGGLDAVLFECRTSDLETLRHCASALANLSLYGGAENQEEMILRKVPMWLFPLAFHNDDNIKYYACLAIAVLVANKEIEAEVLKSGCLDLVEPFVTSHDPSAFARSNLAHAHGQSKHWLKRLVPVLSSNREEARNLAAFHFCMEAGIKREQGNTDIFREINAIEALKNVASCPNAIASKFAAQALRLIGETVPHKLSQQVPLWSVEDVQEWVKQIGFNDYIDKFNESQVDGDLLLKLNQDNLRADIGIGNGILLKRFERELQNLKRMADYSSKDTAKMHQFLSEIGTDYCTYTYAMLNAGIDKCALPHVNEDMLMTECGIHNSIHRLRILNAVKNLENSLPSSSEENMAKTLDVFVSYRRSNGSQLASLLKVHLQLRGFSVFIDVERLEAGKFDNGLLNSIRQAKNFVLVLTPDALHRCINDEDCKDWVHREIVAALNSNCNIIPIIDQQFDWPEVERLPEDMRSVAHFNGVNWIHDYQDACIDKLERFLRGEKNIDRIAAMVPGTPGSVSYQRMHSNDSDYQSGGAGAGSGAGTGGGGGGGVTGSVVDGLMVAANGSGQANHQANRYRQSPSPARQRGSTSQLSGYSRAPSKRSQILTPYRTQQAALLHKTGAGSASMQNMMPLAYLPPRRSSAAGLGHGSGSGMGSGYRSHSVDGLLDQAGSTPEQRIAAAAAKVTAGSTALTNASSTSTLQPEEEVTDAALNDSVTRRDKHTLSPPGNVQQHRKSRSLDHILSKQTLAELLPPSSELADGTQSMQNLAIPMTPQPQRRDTSSSSKSPTPERPPQPAMERVRERQSPEGVSATESEREDQPEECLRHGNQQRASASVHRGASLTSNKTSNSSLGSNFSAGGNNKTIFNRTMKKVRSLIKKP",
   "taxon_ids": [
    "7227"
   ],
   "accessions": [
    "Q6IDD9"
   ],
   "assertions": [
    {
     "assertion_id": "0f57b974e0c77533efdb6aba44687ebdd5964014cd4ba6be30332305fa60970a:GO:0035591",
     "go_id": "GO:0035591",
     "go_name": "signaling adaptor activity",
     "go_definition": "\"The binding activity of a molecule that brings together two or more molecules in a signaling pathway, permitting those molecules to function in a coordinated way. Adaptor molecules themselves do not have catalytic activity.\" [GOC:bf, PMID:19104498]",
     "supporting_assertions": [
      {
       "accession": "Q6IDD9",
       "symbol": "Sarm",
       "relation": "enables",
       "go_id": "GO:0035591",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR039184"
       ],
       "object_name": "NAD(+) hydrolase sarm1",
       "taxon": [
        "taxon:7227"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR039184/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "7d666f4d666cdc7d0c6f09d357923d325b0a34e930d4f8798ffc7dd9b974ed50",
   "sequence": "MMAKKPPKPAPRRIFQERLKITALPLYFEGFLLIKRSGYREYEHYWTELRGTTLFFYTDKKSIIYVDKLDIVDLTCLTEQNSTEKNCAKFTLVLPKEEVQLKTENTESGEEWRGFILTVTELSVPQNVSLLPGQVIKLHEVLEREKKRRIETEQSTSVEKEKEPTEDYVDVLNPMPACFYTVSRKEATEMLQKNPSLGNMILRPGSDSRNYSITIRQEIDIPRIKHYKVMSVGQNYTIELEKPVTLPNLFSVIDYFVKETRGNLRPFICSTDENTGQEPSMEGRSEKLKKNPHIA",
   "taxon_ids": [
    "9606"
   ],
   "accessions": [
    "Q9ULZ2"
   ],
   "assertions": [
    {
     "assertion_id": "7d666f4d666cdc7d0c6f09d357923d325b0a34e930d4f8798ffc7dd9b974ed50:GO:0005068",
     "go_id": "GO:0005068",
     "go_name": "transmembrane receptor protein tyrosine kinase adaptor activity",
     "go_definition": "\"The binding activity of a molecule that brings together a transmembrane receptor protein tyrosine kinase and one or more other molecules, permitting them to function in a coordinated way.\" [GOC:mtg_MIT_16mar07, PMID:10502414, PMID:20565848]",
     "supporting_assertions": [
      {
       "accession": "Q9ULZ2",
       "symbol": "STAP1",
       "relation": "enables",
       "go_id": "GO:0005068",
       "references": [
        "GO_REF:0000107"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB:Q9JM90",
        "ensembl:ENSMUSP00000143251"
       ],
       "object_name": "Signal-transducing adaptor protein 1",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20260724",
       "assigned_by": "Ensembl",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000107.html",
        "https://www.uniprot.org/uniprotkb/Q9JM90/entry"
       ]
      },
      {
       "accession": "Q9ULZ2",
       "symbol": "STAP1",
       "relation": "enables",
       "go_id": "GO:0005068",
       "references": [
        "GO_REF:0000024"
       ],
       "evidence_code": "ISS",
       "with_from": [
        "UniProtKB:Q9JM90"
       ],
       "object_name": "Signal-transducing adaptor protein 1",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20150306",
       "assigned_by": "BHF-UCL",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000024.html",
        "https://www.uniprot.org/uniprotkb/Q9JM90/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "d19e82125e9e71a4e6739fcda2ab13534c374566e2427299792919f589d35363",
   "sequence": "MFQIDLNCDLGESFGAYKIGLDQDILEYVTSANIACGFHAGDPSVMRKTVALAAERGVKMGAHPGLPDLLGFGRRNMAISPEEAYDLVVYQIGALSGFLKAEGLHMQHVKPHGALYNMAAVDQKLSDAIAKAVYKVDPGLILFGLAESELVKAGERIGLQTANEVFADRTYQSDGTLTPRSQPDALIESDDAAVTQVIKMVKEGAVKSQQGHDVSLKADTVCIHGDGAHALTFAQKIRKQLKAAGIEVTAISEQRST",
   "taxon_ids": [
    "224308"
   ],
   "accessions": [
    "P42963"
   ],
   "assertions": [
    {
     "assertion_id": "d19e82125e9e71a4e6739fcda2ab13534c374566e2427299792919f589d35363:GO:0005524",
     "go_id": "GO:0005524",
     "go_name": "ATP binding",
     "go_definition": "\"Binding to ATP, adenosine 5'-triphosphate, a universally important coenzyme and enzyme regulator.\" [ISBN:0198506732]",
     "supporting_assertions": [
      {
       "accession": "P42963",
       "symbol": "pxpA",
       "relation": "enables",
       "go_id": "GO:0005524",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR000031362"
       ],
       "object_name": "5-oxoprolinase subunit A",
       "taxon": [
        "taxon:224308"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000104.html"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "3ef951a363ba111562f14003951f3b98d50125e0f0f3aa182f102baec2f86397",
   "sequence": "MGNAQERPSETIDRERKRLVETLQADSGLLLDALLARGVLTGPEYEALDALPDAERRVRRLLLLVQGKGEAACQELLRCAQRTAGAPDPAWDWQHVGPGYRDRSYDPPCPGHWTPEAPGSGTTCPGLPRASDPDEAGGPEGSEAVQSGTPEEPEPELEAEASKEAEPEPEPEPELEPEAEAEPEPELEPEPDPEPEPDFEERDESEDS",
   "taxon_ids": [
    "9606"
   ],
   "accessions": [
    "O60936",
    "Q5TZN6"
   ],
   "assertions": [
    {
     "assertion_id": "3ef951a363ba111562f14003951f3b98d50125e0f0f3aa182f102baec2f86397:GO:0005102",
     "go_id": "GO:0005102",
     "go_name": "signaling receptor binding",
     "go_definition": "\"Binding to one or more specific sites on a receptor molecule, a macromolecule that undergoes combination with a hormone, neurotransmitter, drug or intracellular messenger to initiate a change in cell function.\" [GOC:bf, GOC:ceb, ISBN:0198506732]",
     "supporting_assertions": [
      {
       "accession": "Q5TZN6",
       "symbol": "NOL3",
       "relation": "enables",
       "go_id": "GO:0005102",
       "references": [
        "GO_REF:0000107"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB:Q9D1X0",
        "ensembl:ENSMUSP00000014920"
       ],
       "object_name": "Nucleolar protein 3",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20260724",
       "assigned_by": "Ensembl",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000107.html",
        "https://www.uniprot.org/uniprotkb/Q9D1X0/entry"
       ]
      },
      {
       "accession": "O60936",
       "symbol": "NOL3",
       "relation": "enables",
       "go_id": "GO:0005102",
       "references": [
        "GO_REF:0000107"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB:Q9D1X0",
        "ensembl:ENSMUSP00000014920"
       ],
       "object_name": "Nucleolar protein 3",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20260724",
       "assigned_by": "Ensembl",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000107.html",
        "https://www.uniprot.org/uniprotkb/Q9D1X0/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "d96234a1cd868e78c03079248ce90371e67ca60f66b7bc7e751dacebe439f5d6",
   "sequence": "MKTSEWIDISQPLNNNIATWPGDTPFSYEVSWSKEESGSVNVGKLTMSIHTGTHIDAPFHFDNDGKKVLDLDVQVYVGPARIIDVSNLESIGKKELESFHLEGVERLLLRTSSHGKAEEFPEVIPHLRADIASFLSEKGIRLIGVDVPSVDPLDDKELAAHHQLFKHGIHILENVVLDHVADGDYELIALPLALTDADGSPVRAVIRPI",
   "taxon_ids": [
    "226900"
   ],
   "accessions": [
    "A0A142GGI7",
    "Q81CK1"
   ],
   "assertions": [
    {
     "assertion_id": "d96234a1cd868e78c03079248ce90371e67ca60f66b7bc7e751dacebe439f5d6:GO:0004328",
     "go_id": "GO:0004328",
     "go_name": "formamidase activity",
     "go_definition": "\"Catalysis of the reaction: formamide + H2O = formate + NH4.\" [RHEA:21948]",
     "supporting_assertions": [
      {
       "accession": "A0A142GGI7",
       "symbol": "kynB",
       "relation": "enables",
       "go_id": "GO:0004328",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR017484"
       ],
       "object_name": "Kynurenine formamidase",
       "taxon": [
        "taxon:1428"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR017484/"
       ]
      },
      {
       "accession": "Q81CK1",
       "symbol": "kynB",
       "relation": "enables",
       "go_id": "GO:0004328",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR017484"
       ],
       "object_name": "Kynurenine formamidase",
       "taxon": [
        "taxon:226900"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR017484/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "0bddb1f8e86058d899eb4380098767c19f42245625f41985478ab33f65d74eaa",
   "sequence": "MSLLTEVETYVLSIIPSGPLKAEIAQRLEDVFAGKNTDLEVLMEWLKTRPILSPLTKGILGFVFTLTVPSERGLQRRRFVQNALNGNGDPNNMDKAVKLYRKLKREITFHGAKEISLSYSAGALASCMGLIYNRMGAVTTEVAFGLVCATCEQIADSQHRSHRQMVTTTNPLIRHENRMVLASTTAKAMEQMAGSSEQAAEAMEVASQARQMVQAMRTIGTHPSSSAGLKNDLLENLQAYQKRMGVQMQRFK",
   "taxon_ids": [
    "211044"
   ],
   "accessions": [
    "A0A0M3YZG2",
    "A0A2R4QPM6",
    "A4ZNR1",
    "B4UPA8",
    "B4URD8",
    "D5F1V4",
    "F0TTD6",
    "F5CGM6",
    "F5CGN7",
    "F5CGP8",
    "J7GV38",
    "P03485"
   ],
   "assertions": [
    {
     "assertion_id": "0bddb1f8e86058d899eb4380098767c19f42245625f41985478ab33f65d74eaa:GO:0005198",
     "go_id": "GO:0005198",
     "go_name": "structural molecule activity",
     "go_definition": "\"The action of a molecule that contributes to the structural integrity of a complex.\" [GOC:mah, GOC:vw]",
     "supporting_assertions": [
      {
       "accession": "F5CGP8",
       "symbol": "M1",
       "relation": "enables",
       "go_id": "GO:0005198",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001561",
        "InterPro:IPR013188",
        "InterPro:IPR015799",
        "InterPro:IPR036039"
       ],
       "object_name": "Matrix protein 1",
       "taxon": [
        "taxon:997555"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001561/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR013188/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR015799/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036039/"
       ]
      },
      {
       "accession": "A4ZNR1",
       "symbol": "M1",
       "relation": "enables",
       "go_id": "GO:0005198",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001561",
        "InterPro:IPR013188",
        "InterPro:IPR015799",
        "InterPro:IPR036039"
       ],
       "object_name": "Matrix protein 1",
       "taxon": [
        "taxon:437606"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001561/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR013188/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR015799/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036039/"
       ]
      },
      {
       "accession": "A0A0M3YZG2",
       "symbol": "M1",
       "relation": "enables",
       "go_id": "GO:0005198",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001561",
        "InterPro:IPR013188",
        "InterPro:IPR015799",
        "InterPro:IPR036039"
       ],
       "object_name": "Matrix protein 1",
       "taxon": [
        "taxon:1686408"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001561/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR013188/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR015799/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036039/"
       ]
      },
      {
       "accession": "B4UPA8",
       "symbol": "M1",
       "relation": "enables",
       "go_id": "GO:0005198",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001561",
        "InterPro:IPR013188",
        "InterPro:IPR015799",
        "InterPro:IPR036039"
       ],
       "object_name": "Matrix protein 1",
       "taxon": [
        "taxon:211044"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001561/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR013188/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR015799/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036039/"
       ]
      },
      {
       "accession": "J7GV38",
       "symbol": "M1",
       "relation": "enables",
       "go_id": "GO:0005198",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001561",
        "InterPro:IPR013188",
        "InterPro:IPR015799",
        "InterPro:IPR036039"
       ],
       "object_name": "Matrix protein 1",
       "taxon": [
        "taxon:1213881"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001561/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR013188/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR015799/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036039/"
       ]
      },
      {
       "accession": "B4URD8",
       "symbol": "M",
       "relation": "enables",
       "go_id": "GO:0005198",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001561",
        "InterPro:IPR013188",
        "InterPro:IPR015799",
        "InterPro:IPR036039"
       ],
       "object_name": "Matrix protein 1",
       "taxon": [
        "taxon:518998"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001561/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR013188/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR015799/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036039/"
       ]
      },
      {
       "accession": "F5CGM6",
       "symbol": "M1",
       "relation": "enables",
       "go_id": "GO:0005198",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001561",
        "InterPro:IPR013188",
        "InterPro:IPR015799",
        "InterPro:IPR036039"
       ],
       "object_name": "Matrix protein 1",
       "taxon": [
        "taxon:997563"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001561/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR013188/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR015799/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036039/"
       ]
      },
      {
       "accession": "A0A2R4QPM6",
       "symbol": "M1",
       "relation": "enables",
       "go_id": "GO:0005198",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001561",
        "InterPro:IPR013188",
        "InterPro:IPR015799",
        "InterPro:IPR036039"
       ],
       "object_name": "Matrix protein 1",
       "taxon": [
        "taxon:11320"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001561/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR013188/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR015799/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036039/"
       ]
      },
      {
       "accession": "F0TTD6",
       "symbol": "M1",
       "relation": "enables",
       "go_id": "GO:0005198",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001561",
        "InterPro:IPR013188",
        "InterPro:IPR015799",
        "InterPro:IPR036039"
       ],
       "object_name": "Matrix protein 1",
       "taxon": [
        "taxon:760066"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001561/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR013188/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR015799/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036039/"
       ]
      },
      {
       "accession": "D5F1V4",
       "symbol": "M1",
       "relation": "enables",
       "go_id": "GO:0005198",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001561",
        "InterPro:IPR013188",
        "InterPro:IPR015799",
        "InterPro:IPR036039"
       ],
       "object_name": "Matrix protein 1",
       "taxon": [
        "taxon:647923"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001561/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR013188/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR015799/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036039/"
       ]
      },
      {
       "accession": "F5CGN7",
       "symbol": "M1",
       "relation": "enables",
       "go_id": "GO:0005198",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001561",
        "InterPro:IPR013188",
        "InterPro:IPR015799",
        "InterPro:IPR036039"
       ],
       "object_name": "Matrix protein 1",
       "taxon": [
        "taxon:997561"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001561/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR013188/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR015799/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036039/"
       ]
      },
      {
       "accession": "P03485",
       "symbol": "M",
       "relation": "enables",
       "go_id": "GO:0005198",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001561",
        "InterPro:IPR013188",
        "InterPro:IPR015799",
        "InterPro:IPR036039"
       ],
       "object_name": "Matrix protein 1",
       "taxon": [
        "taxon:211044"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001561/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR013188/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR015799/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036039/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "e7e37baf98495be5c2c35cb67d85635c640dd97f48a142a7d6b4ad94a4dd9f1a",
   "sequence": "MTEMSFLNSEVLAGDLMSPFDQSGLGAEESLGLLDDYLEVAKHFKPHGFSSDKAGSSEWLAMDGLVSASDTGKEDAFSGTDWMLEKMDLKEFDFDALFRMDDLETMPDELLATLDDTCDLFAPLVQETNKEPPQTVNPIGHLPESVIKVDQAAPFTFLQPLPCSPGFLSSTPDHSFSLELGSEVDISEGDRKPDSAAYITLTPQCVKEEDTPSDSDSGICMSPESYLGSPQHSPSTSRAPPDSLPSPGVPRGSRPKPYDPPGVSVTAKVKTEKLDKKLKKMEQNKTAATRYRQKKRAEQEALTGECKELEKKNEALKEKADSLAKEIQYLKDLIEEVRKARGKKRVP",
   "taxon_ids": [
    "10116"
   ],
   "accessions": [
    "B0BMW3",
    "Q9ES19"
   ],
   "assertions": [
    {
     "assertion_id": "e7e37baf98495be5c2c35cb67d85635c640dd97f48a142a7d6b4ad94a4dd9f1a:GO:0008140",
     "go_id": "GO:0008140",
     "go_name": "cAMP response element binding protein binding",
     "go_definition": "\"Binding to a cAMP response element binding protein (a CREB protein).\" [GOC:mah]",
     "supporting_assertions": [
      {
       "accession": "Q9ES19",
       "symbol": "Atf4",
       "relation": "enables",
       "go_id": "GO:0008140",
       "references": [
        "GO_REF:0000024"
       ],
       "evidence_code": "ISS",
       "with_from": [
        "UniProtKB:Q06507"
       ],
       "object_name": "Cyclic AMP-dependent transcription factor ATF-4",
       "taxon": [
        "taxon:10116"
       ],
       "annotation_date": "20200422",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000024.html",
        "https://www.uniprot.org/uniprotkb/Q06507/entry"
       ]
      },
      {
       "accession": "Q9ES19",
       "symbol": "Atf4",
       "relation": "enables",
       "go_id": "GO:0008140",
       "references": [
        "GO_REF:0000121"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "MGI:MGI:88096"
       ],
       "object_name": "Cyclic AMP-dependent transcription factor ATF-4",
       "taxon": [
        "taxon:10116"
       ],
       "annotation_date": "20200422",
       "assigned_by": "RGD",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000121.html"
       ]
      },
      {
       "accession": "Q9ES19",
       "symbol": "Atf4",
       "relation": "enables",
       "go_id": "GO:0008140",
       "references": [
        "GO_REF:0000107"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB:Q06507",
        "ensembl:ENSMUSP00000105234"
       ],
       "object_name": "Cyclic AMP-dependent transcription factor ATF-4",
       "taxon": [
        "taxon:10116"
       ],
       "annotation_date": "20260724",
       "assigned_by": "Ensembl",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000107.html",
        "https://www.uniprot.org/uniprotkb/Q06507/entry"
       ]
      },
      {
       "accession": "B0BMW3",
       "symbol": "Atf4",
       "relation": "enables",
       "go_id": "GO:0008140",
       "references": [
        "GO_REF:0000107"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB:Q06507",
        "ensembl:ENSMUSP00000105234"
       ],
       "object_name": "Cyclic AMP-dependent transcription factor ATF-4",
       "taxon": [
        "taxon:10116"
       ],
       "annotation_date": "20260724",
       "assigned_by": "Ensembl",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000107.html",
        "https://www.uniprot.org/uniprotkb/Q06507/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "e5831534819e898b0151ede29f36914462d5c8864f4d56f45c5a964011a88c39",
   "sequence": "MSVDNWLLHPLWGQTFLLLLSVAVAQAHWPSEPSEAVRDWKNQLEASMHSVLSDFQEAVPTVVGIPDGTAVVGRSFRVSIPTDLIASSGEIIKVSAAGKEALPSWLHWDPHSHILEGLPLDTDKGVHYISVSAARLGANGSHVPQTSSVFSIEVYPEDHNEPQSVRAASSDPGEVVPSACAADEPVTVLTVILDADLTKMTPKQRIDLLNRMQSFSEVELHNMKLVPVVNNRLFDMSAFMAGPGNAKKVVENGALLSWKLGCSLNQNSVPDIRGVETPAREGAMSAQLGYPVVGWHIANKKPTLPKRLRRQIHATPTPVTAIGPPTTAIQEPPSRIVPTPTSPAIAPPTETMAPPVRDPVPGKPTVTIRTRGAIIQTPTLGPIQPTRVSEAGTTVPGQIRPTLTIPGYVEPTAVITPPTTTTKKPRVSTPKPATPSTDSSTTTTRRPTKKPRTPRPVPRVTTKAPITRLETASPPTRIRTTTSGVPRGGEPNQRPELKNHIDRVDAWVGTYFEVKIPSDTFYDNEDTTTDKLKLTLKLREQQLVGEKSWVQFNSNSQLMYGLPDSSHVGKHEYFMHATDKGGLSAVDAFEIHVHKRPQGDKAPARFKARLAGDPAPVVNDIHKKIALVKKLAFAFGDRNCSSITLQNITRGSIVVEWTNNTLPLEPCPKEQIIGLSRRIADENGKPRPAFSNALEPDFKALSIAVTGSGSCRHLQFIPVAPPSPGSSAAPATEVPDRDPEKSSEDDVYLHTVIPAVVVAAILLIAGIIAMICYRKKRKGKLTLEDQATFIKKGVPIIFADELDDSKPPPSSSMPLILQEEKAPLPPPEYPNQSMPETTPLNQDTVGEYTPLRDEDPNAPPYQPPPPFTAPMEGKGSRPKNMTPYRSPPPYVPP",
   "taxon_ids": [
    "10090"
   ],
   "accessions": [
    "Q544G5",
    "Q62165"
   ],
   "assertions": [
    {
     "assertion_id": "e5831534819e898b0151ede29f36914462d5c8864f4d56f45c5a964011a88c39:GO:0042802",
     "go_id": "GO:0042802",
     "go_name": "identical protein binding",
     "go_definition": "\"Binding to an identical protein or proteins.\" [GOC:jl]",
     "supporting_assertions": [
      {
       "accession": "Q62165",
       "symbol": "Dag1",
       "relation": "enables",
       "go_id": "GO:0042802",
       "references": [
        "GO_REF:0000119"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "UniProtKB:Q14118"
       ],
       "object_name": "Dystroglycan 1",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20251120",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000119.html",
        "https://www.uniprot.org/uniprotkb/Q14118/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "d7ba339757ca6cd2c873e6fbaebf884d910b060fb7dc3726313524a6a2e4d25b",
   "sequence": "MMNQSQRMAPVGSDKELSDLLDFSMMFPLPVANGKGRPASLAGTQFAGSGLEDRPSSESWGNSEQNSSSFDPSRAYSEGAHFSDSHSSLPPSTFLGAGLGGKGSERNAYATFGRDTSVGTLSQAGFLPGELGLSSPGPLSPSGVKSSSQYYTSFPSNPRRRAADGGLDTQPKKVRKVPPGLPSSVYPSSSGDNYSRDATAYPSAKTPSSAYPSPFYVADGSLHPSAELWSPPGQVGFGPMLGDGSAPLPLAPGSSSVSSGAFGGLQQQDRMGYQLHGSEVNGTLPAVSSFSAAPGTYSGTSGHTPPVSGADSLLGTRGTTASSSGDALGKALASIYSPDHSSNNFSPSPSTPVGSPQGLPGTSQWPRAGAPSALSPNYDAGLHGLSKMEDRLDEAIHVLRSHAVGTASELHGLLPGHSTLTTSFAGPMSLGGRHAGLVSGSHPEDGLTSGASLLHNHASLPSQPSSLPDLSQRPPDSFSGLGRAGVTAGASEIKREEKEDEEVTSVADAEEDKKDLKVPRTRTSPDEDEDDLLPPEQKAEREKERRVANNARERLRVRDINEAFKELGRMCQLHLSTEKPQTKLLILHQAVAVILSLEQQVRERNLNPKAACLKRREEEKVSGVVGDPQLALSAAHPGLGEAHNPAGHL",
   "taxon_ids": [
    "10116"
   ],
   "accessions": [
    "P21677"
   ],
   "assertions": [
    {
     "assertion_id": "d7ba339757ca6cd2c873e6fbaebf884d910b060fb7dc3726313524a6a2e4d25b:GO:0070644",
     "go_id": "GO:0070644",
     "go_name": "vitamin D response element binding",
     "go_definition": "\"Binding to a vitamin D response element (VDRE), a short sequence with dyad symmetry found in the promoters of some of the cellular immediate-early genes, regulated by serum.\" [GOC:BHF, GOC:vk, PMID:17426122]",
     "supporting_assertions": [
      {
       "accession": "P21677",
       "symbol": "Tcf3",
       "relation": "enables",
       "go_id": "GO:0070644",
       "references": [
        "GO_REF:0000121"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "UniProtKB:P15923"
       ],
       "object_name": "Transcription factor E2-alpha",
       "taxon": [
        "taxon:10116"
       ],
       "annotation_date": "20191018",
       "assigned_by": "RGD",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000121.html",
        "https://www.uniprot.org/uniprotkb/P15923/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "7134500addbfd6c216b05b1d51df8652a21da0b6ee3626a63c557359b91e1a44",
   "sequence": "MDVTVSELLELFLQSPLVTWVKTFGPFGSGSQDNLTMYMDLVDGIFLNQIMLQIDPRPTNQRINKHVNNDVNLRIQNLTILVRNIKTYYQEVLQQLIVMNLPNVLMIGRDPLSGKSMEEIKKVLLLVLGCAVQCERKEEFIERIKQLDIETQAGIVAHIQEVTHNQENVFDLQWLELPDVAPEELEALSRSMVLHLRRLIDQRDECTELIVDLTQERDYLQAQHPPSPIKSSSADSTPSPTSSLSSEDKQHLAVELADTKARLRRVRQELEDKTEQLVDTRHEVDQLVLELQKVKQENIQLAADARSARAYRDELDSLREKANRVERLELELTRCKEKLHDVDFYKARMEELREDNIILIETKAMLEEQLTAARARGDKVHELEKENLQLKSKLHDLELDRDTDKKRIEELLEENMVLEIAQKQSMNESAHLGWELEQLSKNADLSDASRKSFVFELNECASSRILKLEKENQSLQSTIQGLRDASLVLEESGLKCGELEKENHQLSKKIEKLQTQLEREKQSNQDLETLSEELIREKEQLQSDMETLKADKARQIKDLEQEKDHLNRAMWSLRERSQVSSEARMKDVEKENKALHQTVTEANGKLSQLEFEKRQLHRDLEQAKEKGERAEKLERELQRLQEENGRLARKVTSLETATEKVEALEHESQGLQLENRTLRKSLDTLQNVSLQLEGLERDNKQLDAENLELRRLVETMRFTSTKLAQMERENQQLEREKEELRKNVDLLKALGKKSERLELSYQSVSAENLRLQQSLESSSHKTQTLESELGELEAERQALRRDLEALRLANAQLEGAEKDRKALEQEVAQLEKDKKLLEKEAKRLWQQVELKDAVLDDSTAKLSAVEKESRALDKELARCRDAAGKLKELEKDNRDLTKQVTVHARTLTTLREDLVLEKLKSQQLSSELDKLSQELEKVGLNRELLLQEDDSGSDTKYKILEGRNESALKTTLAMKEEKIVLLEAQMEEKASLNRQLESELQMLKKECETLRQNQGEGQHLQNSFKHPAGKTAASHQGKEAWGPGHKEATMELLRVKDRAIELERNNAALQAEKQLLKEQLQHLETQNVTFSSQILTLQKQSAFLQEHNTTLQTQTAKLQVENSTLSSQSAALTAQYTLLQNHHTAKETENESLQRQQEQLTAAYEALLQDHEHLGTLHERQSAEYEALIRQHSCLKTLHRNLELEHKELGERHGDMLKRKAELEEREKVLTTEREALQQEQRTNALAMGENQRLRGELDRVNFLHHQLKGEYEELHAHTKELKTSLNNAQLELNRWQARFDELKEQHQTMDISLTKLDNHCELLSRLKGNLEEENHHLLSQIQLLSQQNQMLLEQNMENKEQYHEEQKQYIDKLNALRRHKEKLEEKIMDQYKFYDPPPKKKNHWIGAKALVKLIKPKKEGSRERLKSTVDSPPWQLESSDPASPAASQPLRSQAENPDTPALGSNCAEERDAHNGSVGKGPGDLKPKRGSPHRGSLDRTDASTDLAMRSWPSELGSRTCSTSATTTAPSNSTPIARHPGRTKGYNSDDNLCEPSLEFEVPNHRQYVSRPSSLESSRNTSSNSSPLNLKGSSEQLHGRSESFSSEDLIPSRDLATLPREASTPGRNALGRHEYPLPRNGPLPQEGAQKRGTAPPYVGVRPCSASPSSEMVTLEEFLEESNRSSPTHDTPSCRDDLLSDYFRKASDPPAIGGQPGPPAKKEGAKMPTNFVAPTVKMAAPTSEGRPLKPGQYVKPNFRLTEAEAPPSVAPRQAQPPQSLSLGRPRQAPVPPASHAPASRSASLSRAFSLASADLLRASGPEACKQESPQKLGAPEALGGRETGSHTLQSPAPPSSHSLARERTPLVGKAGSSCQGPGPRSRPLDTRRFSLAPPKEERLAPLHQSATAPAIATAGAGAAAAGSGSNSQLLHFSPAAAPAARTKPKAPPRSGEVATITPVRAGLSLSEGDGVPGQGCSEGLPAKSPGRSPDLAPHLGRALEDCSRGSVSKSSPASPEPGGDPQTVWYEYGCV",
   "taxon_ids": [
    "9606"
   ],
   "accessions": [
    "Q9P219"
   ],
   "assertions": [
    {
     "assertion_id": "7134500addbfd6c216b05b1d51df8652a21da0b6ee3626a63c557359b91e1a44:GO:0008017",
     "go_id": "GO:0008017",
     "go_name": "microtubule binding",
     "go_definition": "\"Binding to a microtubule, a filament composed of tubulin monomers.\" [GOC:krc]",
     "supporting_assertions": [
      {
       "accession": "Q9P219",
       "symbol": "CCDC88C",
       "relation": "enables",
       "go_id": "GO:0008017",
       "references": [
        "GO_REF:0000033"
       ],
       "evidence_code": "IBA",
       "with_from": [
        "FB:FBgn0001202",
        "MGI:MGI:1915589",
        "MGI:MGI:1925177",
        "MGI:MGI:1925213",
        "PANTHER:PTN000436622",
        "UniProtKB:Q86VS8"
       ],
       "object_name": "Protein Daple",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20240730",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000033.html",
        "https://www.uniprot.org/uniprotkb/Q86VS8/entry"
       ]
      },
      {
       "accession": "Q9P219",
       "symbol": "CCDC88C",
       "relation": "enables",
       "go_id": "GO:0008017",
       "references": [
        "GO_REF:0000107"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB:Q6VGS5",
        "ensembl:ENSMUSP00000068629"
       ],
       "object_name": "Protein Daple",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20260724",
       "assigned_by": "Ensembl",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000107.html",
        "https://www.uniprot.org/uniprotkb/Q6VGS5/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "4c31c455bf680b23d8a9d3bc621ff7596f1727d9c4aea5b1f1d0b5b4eb7f2b34",
   "sequence": "MDPALSAVRLTVQEAIHILSSSEDAGHILSTLGTLKRYLGGTEDPVLPEEKEEFATVHFSAVLRCLVSKLSPGWLELSPGGQLERLWESFFLDGPPDQAFLVLMEAIESTAGPSFRLMKMAQLLDTFLSTGRVAALMEEQCRPQTKPSFPLFQETLLSKVVGLPDLLGNCLQRDNLTQFFPQNYFPLLGQEVVQALKAVVNFLQDGLDCSVSFVSRVLGKVCIQGRKREILSVLVPQLTVLTQDSCLWQRVCWRLVEQVPDRAVEAVLTGLVEAAPRPEVLSRLLGNLVVKNKKARFVVTRKLLLLQYQHTTPMVQSLLGYLALDSQRRPLLIQVLKELLETWGCSSAVRHTPLEQQCYISKAILVCLAHLGEPELQDIRDELLASMMAGVKCRLDSSLPPVRRLGMIVAEVISSRIHPEGPLLKFQYEDDEMSRELLALATPEPAGDCSSVSRGPSPAPVDTESPVEMPEKAVESDVPPTQPQGSDSELDSDDEFIPYDMSGDRELKSSKEPLYIRDCVEALTTSEDMERWEASLKGLEGLVYRSPTATREVSVELAKVLLHLEEKTCVAEFEQLRQSALVAVTVTDPEQVAKYLTSQFYGLNYSLRQRMDILDVLVLAAQALSRPKSLQRRSQHGPPVPGTMCSPALAVSQTGNVAAPDWQVVVEERIRSKTRRFSKGCPQRELSGVPNEFSSVAGYFFFPLLQHFDRPLVTFDLLGDDQLVLGRLTHTLASLMYLAVNTTVAVPMGKALLEFVWALRFHVDIYVRRGLLSAVSSVLLSVPTERLLGDLPDELLEARSWLADVAEKDVDEDCRELAVRALLLLERLKDKLLSSSSPQP",
   "taxon_ids": [
    "10090"
   ],
   "accessions": [
    "Q9DC40"
   ],
   "assertions": [
    {
     "assertion_id": "4c31c455bf680b23d8a9d3bc621ff7596f1727d9c4aea5b1f1d0b5b4eb7f2b34:GO:0140777",
     "go_id": "GO:0140777",
     "go_name": "protein-containing complex stabilizing activity",
     "go_definition": "\"A molecular function that involves direct binding to one of the subunits of a protein-containing complex, thus preventing an interaction with a factor that would promote dissociation of the complex.\" [PMID:34040253]",
     "supporting_assertions": [
      {
       "accession": "Q9DC40",
       "symbol": "Telo2",
       "relation": "enables",
       "go_id": "GO:0140777",
       "references": [
        "PMID:24657168"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "UniProtKB:Q9Y4R8"
       ],
       "object_name": "Telomere length regulation protein TEL2 homolog",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20230106",
       "assigned_by": "MGI",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://pubmed.ncbi.nlm.nih.gov/24657168/",
        "https://www.uniprot.org/uniprotkb/Q9Y4R8/entry"
       ]
      },
      {
       "accession": "Q9DC40",
       "symbol": "Telo2",
       "relation": "enables",
       "go_id": "GO:0140777",
       "references": [
        "GO_REF:0000107"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB:Q9Y4R8",
        "ensembl:ENSP00000262319"
       ],
       "object_name": "Telomere length regulation protein TEL2 homolog",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20260724",
       "assigned_by": "Ensembl",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000107.html",
        "https://www.uniprot.org/uniprotkb/Q9Y4R8/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "4fe5a2a78853c457b011250812a422409cdbd80140b051ca12d05e39c7c61a08",
   "sequence": "MTLLSGANSVGYPSQMPAKRISWRRLVALLPNSIGLSLKPVTPHLRDLHLVQRRLQDDFLGCEALWTIFMAAAWSYRYRTRLRPFPSHWNTIDLVFNTLGDAPRLEVLQQQLMHCDYQACSPNVVRLLTDILVDQADRVSLSSLRPCEFQELYAHLGMSPPKQPPTQIFEVRTGKGNEKGEAYASLRQENKESVRLGFYGCKLEKVYALLNSQNSLDNGKYLELTCDINEALARSKPQAGVGGSRCGSILRCVAVVEFVFQDNETSRDKKQVIIKDANTMQVSYLLLYGQSCNEYAMERQIKLMAEPARELLNWLERYHKKAISLGVGLLIVSSMAHFGSTFFRLLARTGFYVFKRGIL",
   "taxon_ids": [
    "7227"
   ],
   "accessions": [
    "A1ZAG4"
   ],
   "assertions": [
    {
     "assertion_id": "4fe5a2a78853c457b011250812a422409cdbd80140b051ca12d05e39c7c61a08:GO:0003950",
     "go_id": "GO:0003950",
     "go_name": "NAD+ poly-ADP-ribosyltransferase activity",
     "go_definition": "\"Catalysis of the reaction: NAD+ + (ADP-D-ribosyl)(n)-acceptor = nicotinamide + (ADP-D-ribosyl)(n+1)-acceptor.\" [EC:2.4.2.30]",
     "supporting_assertions": [
      {
       "accession": "A1ZAG4",
       "symbol": "Parp16",
       "relation": "enables",
       "go_id": "GO:0003950",
       "references": [
        "GO_REF:0000024"
       ],
       "evidence_code": "ISS",
       "with_from": [
        "UniProtKB:Q8N5Y8"
       ],
       "object_name": "Protein mono-ADP-ribosyltransferase Parp16",
       "taxon": [
        "taxon:7227"
       ],
       "annotation_date": "20220912",
       "assigned_by": "FlyBase",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000024.html",
        "https://www.uniprot.org/uniprotkb/Q8N5Y8/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "6f04b325eaafec58879bf2a89c083ed547a72940138491b6ad151527b390ab4f",
   "sequence": "MSGSFWLLLSLVAVTAAQSTTEDEAKKFLDKFNSKAEDLSYESSLASWDYNTNISDENVQKMDEAGAKWSAFYEEQSKLAKNYPLEEIQNDTVKRQLQILQQSGSPVLSEDKSKRLNSILNAMSTIYSTGKVCKPNNPQECLLLEPGLDNIMGTSKDYNERLWAWEGWRAEVGKQLRPLYEEYVVLKNEMARGYHYEDYGDYWRRDYETEESSGPGYSRDQLMKDVDRIFTEIKPLYEHLHAYVRAKLMDTYPLHISPTGCLPAHLLGDMWGRFWTNLYPLTVPFGQKPNIDVTDAMLNQGWDANRIFKEAEKFFVSVSLPKMTEGFWNKSMLTEPGDGRKVVCHPTAWDLGKGDFRIKMCTKVTMEDFLTAHHEMGHIQYDMAYASQPYLLRNGANEGFHEAVGEVMSLSVATPKHLKTMGLLSPDFREDDETEINFLLKQALNIVGTLPFTYMLEKWRWMVFKGEIPKEEWMKKWWEMKREIVGVVEPVPHDETYCDPASLFHVANDYSFIRYYTRTIFEFQFHEALCRIAQHNGPLHKCDISNSTDAGKKLHQMLSVGKSQAWTKTLEDIVDSRNMDVGPLLRYFKPLYTWLQEQNRKSYVGWNTDWSPYADQSIKVWISLKSALGEKAYEWNDNEMYLFRSSVAYAMREYFLKTKNQTILFGDENVWVSNLKPRISFNFHVTSPGNVSDIIPRPEVEGAIRMSRSRINDAFRLDDNSLEFLGIQPTLGPPYQPPVTIWLIVFGVVMAVVVVGIVVLIITGIRDRRKTDQARSEENPYPSVDLSKGENNPGFQNGDDVQTSF",
   "taxon_ids": [
    "196889"
   ],
   "accessions": [
    "E2DHI3"
   ],
   "assertions": [
    {
     "assertion_id": "6f04b325eaafec58879bf2a89c083ed547a72940138491b6ad151527b390ab4f:GO:0004180",
     "go_id": "GO:0004180",
     "go_name": "carboxypeptidase activity",
     "go_definition": "\"Catalysis of the hydrolysis of a single C-terminal amino acid residue from a polypeptide chain.\" [https://www.ebi.ac.uk/merops/about/glossary.shtml#CARBOXYPEPTIDASE]",
     "supporting_assertions": [
      {
       "accession": "E2DHI3",
       "symbol": "ACE2",
       "relation": "enables",
       "go_id": "GO:0004180",
       "references": [
        "GO_REF:0000117"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "ARBA:ARBA00028164"
       ],
       "object_name": "Angiotensin-converting enzyme",
       "taxon": [
        "taxon:196889"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000117.html"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "2be1462ce9456241aaaa2421ee07e9d854f76c3abe07648af8558bfef543c8a4",
   "sequence": "MNCARLHQRIPLRAMALTTTSYPALAPSHSFANASTSVMTASAMAVAHRAAATTGASDATTTAQVVHYSKLGMDGFTDIKFNTAADELLELVETKVDALDSAAVEDVSCNGGVLTLETTERGTFILNKQAPNVQLWLSSPISGPHHYDMITVTQDGHEKVSWKSDHDGHDLVKKLEKELTEVLGTAFKL",
   "taxon_ids": [
    "5661"
   ],
   "accessions": [
    "A4I1E7",
    "I7HD01"
   ],
   "assertions": [
    {
     "assertion_id": "2be1462ce9456241aaaa2421ee07e9d854f76c3abe07648af8558bfef543c8a4:GO:0051537",
     "go_id": "GO:0051537",
     "go_name": "2 iron, 2 sulfur cluster binding",
     "go_definition": "\"Binding to a 2 iron, 2 sulfur (2Fe-2S) cluster; this cluster consists of two iron atoms, with two inorganic sulfur atoms found between the irons and acting as bridging ligands.\" [GOC:ai, PMID:15952888, Wikipedia:Iron-sulfur_cluster]",
     "supporting_assertions": [
      {
       "accession": "I7HD01",
       "symbol": "LDBPK_251090",
       "relation": "enables",
       "go_id": "GO:0051537",
       "references": [
        "GO_REF:0000118"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "PANTHER:PTN000423756"
       ],
       "object_name": "Frataxin-like protein, mitochondrial",
       "taxon": [
        "taxon:5661"
       ],
       "annotation_date": "20260727",
       "assigned_by": "TreeGrafter",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000118.html"
       ]
      },
      {
       "accession": "A4I1E7",
       "symbol": "LINJ_25_1090",
       "relation": "enables",
       "go_id": "GO:0051537",
       "references": [
        "GO_REF:0000118"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "PANTHER:PTN000423756"
       ],
       "object_name": "ferroxidase",
       "taxon": [
        "taxon:5671"
       ],
       "annotation_date": "20260727",
       "assigned_by": "TreeGrafter",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000118.html"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "f47b889573e76ad26980f54cd6c6980dabd22cc9f2017babfe76f63851ef4e73",
   "sequence": "MSANDILDVLNIQRDESNQPPKKKQKSSSTPTLPDGKQLTGMARELYNLVGPNTPPINLNSNSYTANKEKMKKFKPSPWTRMPFTPKQGIELNHWVKGSKELIEQQEFEEDGTPKPYFFEKYNVQLEIPEFVDEDTYDLYMIEIKEYESKMKEERARREKERKEREKRDLEEKKKKQQQQQQKSQQNPQNQIKDDEKNQDTRNNTDKKDSEQKSEDKPTVEAKKETDEKKDDVVLKDNTNETKPVTETTKSETETTEQNNSEKTNENETNKTNDKDGEGNLTKSKDSATEDQSNNKKENDEDTESEWTYKETKHLFELCQAFELKWPIIHDRFPNPNRTAEDLKEQFYRICIKILENQKNKNQALIDSLKAYCKPRELERKQYLENLLKRTPAEIAEEESLVIEARRFEIAAKKMLMERSNLLTLLDSPQTTQNVSQYQSSQGITNLYNNLLIYDKHQKKKQMANKSNPQQEPVPPPIPLAASSSVKRDRGFQTQLQQYLSSFLKQNHHTNPAVKQEINSIQQLLMKRLTQKEEEAYGLYFHGTEKLNPGVMLRSQQKLPGLNQRQSILKSVNILLQEMDIPTGGGTSWKPIMPTRKTMAKYDELIRSVVTLLDVKKAKDKLESEIKLIKSQRGL",
   "taxon_ids": [
    "237561"
   ],
   "accessions": [
    "Q5AAJ7"
   ],
   "assertions": [
    {
     "assertion_id": "f47b889573e76ad26980f54cd6c6980dabd22cc9f2017babfe76f63851ef4e73:GO:0003714",
     "go_id": "GO:0003714",
     "go_name": "transcription corepressor activity",
     "go_definition": "\"A transcription coregulator activity that represses or decreases the transcription of specific gene sets via binding to a DNA-binding transcription factor at a specific genomic locus, either on its own or as part of a complex. Corepressors often act by altering chromatin structure and modifications. For example, one class of transcription corepressors modifies chromatin structure through covalent modification of histones. A second class remodels the conformation of chromatin in an ATP-dependent fashion. A third class modulates interactions of DNA-bound DNA-binding transcription factors with other transcription coregulators.\" [GOC:txnOH-2018, PMID:10213677, PMID:16858867]",
     "supporting_assertions": [
      {
       "accession": "Q5AAJ7",
       "symbol": "SWC4",
       "relation": "enables",
       "go_id": "GO:0003714",
       "references": [
        "GO_REF:0000033"
       ],
       "evidence_code": "IBA",
       "with_from": [
        "MGI:MGI:1913483",
        "PANTHER:PTN000304415"
       ],
       "object_name": "SWR1-complex protein 4",
       "taxon": [
        "taxon:237561"
       ],
       "annotation_date": "20170228",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000033.html"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "ef7ea7231d382877816d0dfb3804c51321936a6c0eb61eeb2a660c9c9da95c90",
   "sequence": "MAFQLALDALAPTTHRDPSLHPILESTVDSIRSSIQTYPWSIPKELLPLLNSYGIPTSGLGTSHHPHAAHKTIETFLLCTHWSFQATTPSSVMFMKPSKFNKLAQVNSNFRELKNYRLHPNDSTRYPFTSPDLPVFPTIFMHDALMYYHPSQIMDLFLRKPNLERLYASLVVPPEAHLSDQSFYPKLYTYTTTRHTLHYVPEGHEAGSYNQPSDAHSWLRINSIRLGNHHLSVTILESWGPVHSLLIQRGTPPPDPSLQAPPTLMTSDLFRSYQEPRLDVVSFRIPDAIELPQATFLQQPLRDRLVPRAVYNALFTYTRAVRTLRTSDPAAFVRMHSSKPDHDWVTSNAWDNLQTFALLNVPLRPNVVYHVLQSPIASLSLYLRQHWRRLTATAVPILSFLTLLQRFLPLPIPLAEVKSITAFRRELYRKKEPHHPLDVFHLQHRVRNYHSAISAVRPASPPHQKLPHALQKAALLLLRPISPLLTATPFFRSEQKSMLPNAELSWTLKRFALPWQASLVLLALSESSILLHKLFSPPTLQAQHDTYHRHLHPGSYSLQWERTPLSIPRTTAFLPFTPTTSTAPPDRSEASLPPAFASTFVPRPPPAASSPGAQPPTTTAAPPTPIEPTQRTHQNSDLALESSTSTEPPPPPIRSPDMTPSAPVLFPEINSPRRFPPQLPATPDLEPAHTPPPLSIPHQDPTDSADPLMGSHLLHHSLPAPPTHPLPSSQLLPAPLTNDPTAIGPVLPFEELHPRRYPENTATFLTRLRSLPSNHLPQPTLNCLLSAVSDQTKVSEEHLWESLQTILPDSQLSNEETNTLGLSTEHLTALAHLYNFQATVYSDRGPILFGPSDTIKRIDITHTTGPPSHFSPGKRLLGSQPSAKGHPSDPLIRAMKSFKVSGNYLPFSEAHNHPTSISHAKNLISNMKNGFDGVLSLLDVSTGQRTGPTPKERIIQIDHYLDTNPGKTTPVVHFAGFAGCGKTYPIQQLLKTKLFKDFRVSCPTTELRTEWKTAMELHGSQSWRFNTWESSILKSSRILVIDEIYKMPRGYLDLSILADPALELVIILGDPLQGEYHSQSKDSSNHRLPSETLRLLPYIDMYCWWSYRIPQCIARLFQIHSFNAWQGVIGSVSTPHDQSPVLTNSHASSLTFNSLGYRSCTISSSQGLTFCDPAIIVLDNYTKWLSSANGLVALTRSRSGVQFMGPSSYVGGTNGSSAMFSDAFNNSLIIMDRYFPSLFPQLKLITSPLTTRGPKLNGATPSASPTHRSPNFHLPPHIPLSYDRDFVTVNPTLPDQGPETRLDTHFLPPSRLPLHFDLPPAITPPPVSTSVDPPQAKASPVYPGEFFDSLAAFFLPAHDPSTREILHKDQSSNQFPWFDRPFSLSCQPSSLISAKHAPNHDPTLLPASINKRLRFRPSDSPHQITADDVVLGLQLFHSLCRAYSRQPNSTVPFNPELFAECISLNEYAQLSSKTQSTIVANASRSDPDWRHTTVKIFAKAQHKVNDGSIFGSWKACQTLALMHDYVILVLGPVKKYQRIFDNADRPPNIYSHCGKTPNQLRDWCQEHLTHSTPKIANDYTAFDQSQHGESVVLEALKMKRLNIPSHLIQLHVHLKTNVSTQFGPLTCMRLTGEPGTYDDNTDYNLAVIYSQYDVGSCPIMVSGDDSLIDHPLPTRHDWPSVLKRLHLRFKLELTSHPLFCGYYVGPAGCIRNPLALFCKLMIAVDDDALDDRRLSYLTEFTTGHLLGESLWHLLPETHVQYQSACFDFFCRRCPRHEKMLLDDSTPALSLLERITSSPRWLTKNAMYLLPAKLRLAITSLSQTQSFPESIEVSHAESELLHYVQ",
   "taxon_ids": [
    "12154"
   ],
   "accessions": [
    "P10358"
   ],
   "assertions": [
    {
     "assertion_id": "ef7ea7231d382877816d0dfb3804c51321936a6c0eb61eeb2a660c9c9da95c90:GO:0034062",
     "go_id": "GO:0034062",
     "go_name": "5'-3' RNA polymerase activity",
     "go_definition": "\"Catalysis of the reaction: nucleoside triphosphate + RNA(n) = diphosphate + RNA(n+1); the synthesis of RNA from ribonucleotide triphosphates in the presence of a nucleic acid template, via extension of the 3'-end.\" [GOC:mah, GOC:pf, RHEA:21248]",
     "supporting_assertions": [
      {
       "accession": "P10358",
       "symbol": "P10358",
       "relation": "enables",
       "go_id": "GO:0034062",
       "references": [
        "GO_REF:0000116"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "RHEA:21248"
       ],
       "object_name": "Non-structural replication polyprotein",
       "taxon": [
        "taxon:12154"
       ],
       "annotation_date": "20260727",
       "assigned_by": "RHEA",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000116.html"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "2936953b9659b27b2b76e0c08793623202770017858bc6c194a2e9bf016b4b71",
   "sequence": "MLPLRRSLLKPLQDVLRHSHRQMSQKNLLELTDRGFFHGIFPDTAAPRMKQLFTRGQQSIYAGFDPTADSLHVGNLLVIMGLIHCQRAGHRPIALVGGATGLIGDPSGRKTERNQLGETVIETNLKAIEQQLRRVFENHENCLWDSKKQKLPLAPLIIVNNADWYADLQLIDFVANMGRHFRMGSMLSRSSVQSRLESEDGMSFTEFTYQIFQAYDWLHLLRRHNCCFQMGGSDQTGNLMTGHELISRVERKREVFGLTLPLVTTEEGDKFGKSAGNAVWLDGNKTSPFALYQFFLRMPDSEVEKLLKLFTFIPLPQVEQLMREHTKEPEKRKAQTLLAEDVTLLVHGESGLKQAERVTNALYKGNVEGLAELNLSEIQQTFQGATMVNLLTEPGMSILELAMKAKCFPTETDAVRIINAGGFYVNQKRVQNIAEVLTTGVHILRNGISLLRVGKRNFYIVRWQ",
   "taxon_ids": [
    "7227"
   ],
   "accessions": [
    "Q9W107"
   ],
   "assertions": [
    {
     "assertion_id": "2936953b9659b27b2b76e0c08793623202770017858bc6c194a2e9bf016b4b71:GO:0000166",
     "go_id": "GO:0000166",
     "go_name": "nucleotide binding",
     "go_definition": "\"Binding to a nucleotide, any compound consisting of a nucleoside that is esterified with (ortho)phosphate or an oligophosphate at any hydroxyl group on the ribose or deoxyribose.\" [GOC:mah, ISBN:0198547684]",
     "supporting_assertions": [
      {
       "accession": "Q9W107",
       "symbol": "TyrRS-m",
       "relation": "enables",
       "go_id": "GO:0000166",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001412",
        "InterPro:IPR002305",
        "InterPro:IPR002307"
       ],
       "object_name": "Tyrosine--tRNA ligase, mitochondrial",
       "taxon": [
        "taxon:7227"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001412/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR002305/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR002307/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "69706cd52d6f574832f2271bd3f5f0f142bbfa0b0078f9e7134f9efd8aea867b",
   "sequence": "MEELVVEVRGSNGAFYKAFVKDVHEDSITVAFENNWQPDRQIPFHDVRFPPPVGYNKDINESDEVEVYSRANEKEPCCWWLAKVRMIKGEFYVIEYAACDATYNEIVTIERLRSVNPNKPATKDTFHKIKLDVPEDLRQMCAKEAAHKDFKKAVGAFSVTYDPENYQLVILSINEVTSKRAHMLIDMHFRSLRTKLSLIMRNEEASKQLESSRQLASRFHEQFIVREDLMGLAIGTHGANIQQARKVPGVTAIDLDEDTCTFHIYGEDQDAVKKARSFLEFAEDVIQVPRNLVGKVIGKNGKLIQEIVDKSGVVRVRIEAENEKNVPQEEEIMPPNSLPSNNSRVGPNAPEEKKHLDIKENSTHFSQPNSTKVQRVLVASSVVAGESQKPELKAWQGMVPFVFVGTKDSIANATVLLDYHLNYLKEVDQLRLERLQIDEQLRQIGASSRPPPNRTDKEKSYVTDDGQGMGRGSRPYRNRGHGRRGPGYTSGTNSEASNASETESDHRDELSDWSLAPTEEERESFLRRGDGRRRGGGGRGQGGRGRGGGFKGNDDHSRTDNRPRNPREAKGRTTDGSLQIRVDCNNERSVHTKTLQNTSSEGSRLRTGKDRNQKKEKPDSVDGQQPLVNGVP",
   "taxon_ids": [
    "9606"
   ],
   "accessions": [
    "Q06787"
   ],
   "assertions": [
    {
     "assertion_id": "69706cd52d6f574832f2271bd3f5f0f142bbfa0b0078f9e7134f9efd8aea867b:GO:0003682",
     "go_id": "GO:0003682",
     "go_name": "chromatin binding",
     "go_definition": "\"Binding to chromatin, the network of fibers of DNA, protein, and sometimes RNA, that make up the chromosomes of the eukaryotic nucleus during interphase.\" [GOC:jl, ISBN:0198506732, PMID:20404130]",
     "supporting_assertions": [
      {
       "accession": "Q06787",
       "symbol": "FMR1",
       "relation": "enables",
       "go_id": "GO:0003682",
       "references": [
        "GO_REF:0000024"
       ],
       "evidence_code": "ISS",
       "with_from": [
        "UniProtKB:P35922"
       ],
       "object_name": "Fragile X messenger ribonucleoprotein 1",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20160418",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000024.html",
        "https://www.uniprot.org/uniprotkb/P35922/entry"
       ]
      }
     ]
    }
   ]
  }
 ]
}
```
