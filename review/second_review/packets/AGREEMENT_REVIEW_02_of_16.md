# Protein function review — batch 2 of 16

25 proteins; 25 assertions. Name your output `AGREEMENT_RESULTS_02.json`.

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
 "batch": 2,
 "source_release": "GOA 2026-07-29",
 "proteins": [
  {
   "sequence_sha256": "0b0b625f889ca27c25c2ea760f56492fa6b1b20e5dc5be156fb4715281fd4f43",
   "sequence": "MEDSGKTFSSEEEEANYWKDLAMTYKQRAENTQEELREFQEGSREYEAELETQLQQIETRNRDLLSENNRLRMELETIKEKFEVQHSEGYRQISALEDDLAQTKAIKDQLQKYIRELEQANDDLERAKRATIMSLEDFEQRLNQAIERNAFLESELDEKENLLESVQRLKDEARDLRQELAVQQKQEKPRTPMPSSVEAERTDTAVQATGSVPSTPIAHRGPSSSLNTPGSFRRGLDDSTGGTPLTPAARISALNIVGDLLRKVGALESKLASCRNLVYDQSPNRTGGPASGRSSKNRDGGERRPSSTSVPLGDKGLDTSCRWLSKSTTRSSSSC",
   "taxon_ids": [
    "9606"
   ],
   "accessions": [
    "A0A2J8R0D9",
    "Q9NXR1",
    "X5DR54"
   ],
   "assertions": [
    {
     "assertion_id": "0b0b625f889ca27c25c2ea760f56492fa6b1b20e5dc5be156fb4715281fd4f43:GO:0008017",
     "go_id": "GO:0008017",
     "go_name": "microtubule binding",
     "go_definition": "\"Binding to a microtubule, a filament composed of tubulin monomers.\" [GOC:krc]",
     "supporting_assertions": [
      {
       "accession": "Q9NXR1",
       "symbol": "NDE1",
       "relation": "enables",
       "go_id": "GO:0008017",
       "references": [
        "GO_REF:0000033"
       ],
       "evidence_code": "IBA",
       "with_from": [
        "MGI:MGI:1914453",
        "MGI:MGI:1932915",
        "PANTHER:PTN000100519"
       ],
       "object_name": "Nuclear distribution protein nudE homolog 1",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20170228",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000033.html"
       ]
      },
      {
       "accession": "Q9NXR1",
       "symbol": "NDE1",
       "relation": "enables",
       "go_id": "GO:0008017",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR033494"
       ],
       "object_name": "Nuclear distribution protein nudE homolog 1",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR033494/"
       ]
      },
      {
       "accession": "Q9NXR1",
       "symbol": "NDE1",
       "relation": "enables",
       "go_id": "GO:0008017",
       "references": [
        "GO_REF:0000107"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB:Q9CZA6",
        "ensembl:ENSMUSP00000023359"
       ],
       "object_name": "Nuclear distribution protein nudE homolog 1",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20260724",
       "assigned_by": "Ensembl",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000107.html",
        "https://www.uniprot.org/uniprotkb/Q9CZA6/entry"
       ]
      },
      {
       "accession": "Q9NXR1",
       "symbol": "NDE1",
       "relation": "enables",
       "go_id": "GO:0008017",
       "references": [
        "GO_REF:0000024"
       ],
       "evidence_code": "ISS",
       "with_from": [
        "UniProtKB:Q9CZA6"
       ],
       "object_name": "Nuclear distribution protein nudE homolog 1",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20060616",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000024.html",
        "https://www.uniprot.org/uniprotkb/Q9CZA6/entry"
       ]
      },
      {
       "accession": "A0A2J8R0D9",
       "symbol": "NDE1",
       "relation": "enables",
       "go_id": "GO:0008017",
       "references": [
        "GO_REF:0000118"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "PANTHER:PTN002589138"
       ],
       "object_name": "NUDE domain-containing protein",
       "taxon": [
        "taxon:9601"
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
       "accession": "A0A2J8R0D9",
       "symbol": "NDE1",
       "relation": "enables",
       "go_id": "GO:0008017",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR033494"
       ],
       "object_name": "NUDE domain-containing protein",
       "taxon": [
        "taxon:9601"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR033494/"
       ]
      },
      {
       "accession": "X5DR54",
       "symbol": "NDE1",
       "relation": "enables",
       "go_id": "GO:0008017",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR033494"
       ],
       "object_name": "NUDE domain-containing protein",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR033494/"
       ]
      },
      {
       "accession": "X5DR54",
       "symbol": "NDE1",
       "relation": "enables",
       "go_id": "GO:0008017",
       "references": [
        "GO_REF:0000107"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB:Q9CZA6",
        "ensembl:ENSMUSP00000023359"
       ],
       "object_name": "NUDE domain-containing protein",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20260724",
       "assigned_by": "Ensembl",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000107.html",
        "https://www.uniprot.org/uniprotkb/Q9CZA6/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "400fced191c9a75c735f5a60a73c29548cb733d002b21393f4f744613539f1a7",
   "sequence": "MMLPRTWAVETEQLLTPNPGYGTQAGTSPAPTTPPEEEDLRRRLKYFFMSPCDKFRAKGRKPCKLMLQVVKILVVTVQLILFGLSNQLVVTFREENTIAFRHLFLLGYSDGSDDTFAAYTREQLYQAIFYAVDQYLTLPEISLGRYAYVRGGGGPWANGSALALCQRYYHRGHVDPANDTFDIDPRVVTDCIQVDPPDRPPDIPSEDLDFLDGSTSYKNITLKFHKLINVTIHFQLKTINLQSLINNEIPDCYTFSILITFDNKAHSGRIPIRLETQTHIQECKHPSVSRHGDNSFRLLFDVVVILTCSLSFLLCARSLLRGFLLQNEFVVFMWRRRGREISLWERLEFVNGWYILLVTSDVLTISGTIMKIGIEAKNLASYDVCSILLGTSTLLVWVGVIRYLTFFHKYNILIATLRVALPSVMRFCCCVAVIYLGYCFCGWIVLGPYHVKFRSLSMVSECLFSLINGDDMFVTFAAMQAQQGRSSLVWLFSQLYLYSFISLFIYMVLSLFIALITGAYDTIKHPGGTGTEKSELQAYIEQCQDSPTSGKFRRGSGSTCSLLCCCGRNSPEDHSLLVN",
   "taxon_ids": [
    "10116"
   ],
   "accessions": [
    "A0A8I6G3L4"
   ],
   "assertions": [
    {
     "assertion_id": "400fced191c9a75c735f5a60a73c29548cb733d002b21393f4f744613539f1a7:GO:0005381",
     "go_id": "GO:0005381",
     "go_name": "iron ion transmembrane transporter activity",
     "go_definition": "\"Enables the transfer of iron (Fe) ions from one side of a membrane to the other.\" [GOC:ai]",
     "supporting_assertions": [
      {
       "accession": "A0A8I6G3L4",
       "symbol": "Mcoln1",
       "relation": "enables",
       "go_id": "GO:0005381",
       "references": [
        "GO_REF:0000121"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "MGI:MGI:1890498"
       ],
       "object_name": "Mucolipin-1",
       "taxon": [
        "taxon:10116"
       ],
       "annotation_date": "20240319",
       "assigned_by": "RGD",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000121.html"
       ]
      },
      {
       "accession": "A0A8I6G3L4",
       "symbol": "Mcoln1",
       "relation": "enables",
       "go_id": "GO:0005381",
       "references": [
        "GO_REF:0000121"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "UniProtKB:Q9GZU1"
       ],
       "object_name": "Mucolipin-1",
       "taxon": [
        "taxon:10116"
       ],
       "annotation_date": "20240319",
       "assigned_by": "RGD",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000121.html",
        "https://www.uniprot.org/uniprotkb/Q9GZU1/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "05cc820efa81a7b3a62924577ca010a38119277947a0aa23832f77eca482fc4c",
   "sequence": "MTESNRRDISWIFALLAAYFVLQVGVRLATSHSLDLDEAEQAFRSQWLAAGYGPQPPFYNWLQYTVFQFAGVSLTALSIVKNLLLFISYLLYGLTARLVLRDKALVAIATLGLLTIPQMAFEMQRDLTHTVAVFFSASIFFYGFIRSLKQPSLASYLIAGIGIGFGLLAKYNFAILPAAALIAALSDARLRPRIFDWRLGLTAAVALVITLPHLFWLKDNLDFATARTLEKMTASGDASYLTQVAMGVSSLALAIISFAALTVAVFAIVFGKSLRPALGSGSEWTRLLERMMLVFLAGILLLIVFGGAAGIKDRWLVPMLFILPLYFCLKIEAAGVETGKALRRFIPVVAVIMIGVPAALYGSVAAARFTGHYERLNRPYAGMLEILRKQAEPAAILAGDSLLAGNLRQDIPGVPILSADYPGFNPDLTSRRPLLLVWLLPKGGSEALPPDMAEWLQANLGTSAPEASVIDVPYFYGRGDDRYRFGYAWVNQPG",
   "taxon_ids": [
    "216596"
   ],
   "accessions": [
    "Q1MJ97"
   ],
   "assertions": [
    {
     "assertion_id": "05cc820efa81a7b3a62924577ca010a38119277947a0aa23832f77eca482fc4c:GO:0016763",
     "go_id": "GO:0016763",
     "go_name": "pentosyltransferase activity",
     "go_definition": "\"Catalysis of the transfer of a pentosyl group from one compound (donor) to another (acceptor).\" [GOC:jl]",
     "supporting_assertions": [
      {
       "accession": "Q1MJ97",
       "symbol": "rgtB",
       "relation": "enables",
       "go_id": "GO:0016763",
       "references": [
        "GO_REF:0000118"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "PANTHER:PTN008018715"
       ],
       "object_name": "Lipopolysaccharide core galacturonosyltransferase RgtB",
       "taxon": [
        "taxon:216596"
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
   "sequence_sha256": "1a58512a9c22a39acc9ac9b2ae748fe59e57e45487d775d7274874bec776e1ae",
   "sequence": "MKPSHSSCEAAPLLPNMAETHYAPLSSAFPFVTSYQTGSSRLPEVSRSTERALREGKLLELVYGIKETVATLSQIPVSIFVTGDSGNGMSSFINALRVIGHDEDASAPTGVVRTTKTRTEYSSSHFPNVVLWDLPGLGATAQTVEDYVEEMKFSTCDLFIIIASEQFSSNHVKLSKIIQSMGKRFYIVWTKLDRDLSTSVLSEVRLLQNIQENIRENLQKEKVKYPPVFLVSSLDPLLYDFPKLRDTLHKDLSNIRCCEPLKTLYGTYEKIVGDKVAVWKQRIANESLKNSLGVRDDDNMGECLKVYRLIFGVDDESVQQVAQSMGTVVMEYKDNMKSQNFYTLRREDWKLRLMTCAIVNAFFRLLRFLPCVCCCLRRLRHKRMLFLVAQDTKNILEKILRDSIFPPQI",
   "taxon_ids": [
    "10090"
   ],
   "accessions": [
    "A0A8C6G445",
    "Q60766"
   ],
   "assertions": [
    {
     "assertion_id": "1a58512a9c22a39acc9ac9b2ae748fe59e57e45487d775d7274874bec776e1ae:GO:0019901",
     "go_id": "GO:0019901",
     "go_name": "protein kinase binding",
     "go_definition": "\"Binding to a protein kinase, any enzyme that catalyzes the transfer of a phosphate group, usually from ATP, to a protein substrate.\" [GOC:jl]",
     "supporting_assertions": [
      {
       "accession": "Q60766",
       "symbol": "Irgm1",
       "relation": "enables",
       "go_id": "GO:0019901",
       "references": [
        "GO_REF:0000119"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "UniProtKB:A1A4Y4"
       ],
       "object_name": "Immunity-related GTPase family M protein 1",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20260702",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000119.html",
        "https://www.uniprot.org/uniprotkb/A1A4Y4/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "f9e5e2678617ec09bbb5cc8c5662525b1e370b18f0fca7b264a8ed8819733bef",
   "sequence": "MAVAGAVSGEPLVHWCTQQLRKTFGLDVSEEIIQYVLSIESAEEIREYVTDLLQGNEGKKGQFIEELITKWQKNDQELISDPLQQCFKKDEILDGQKSGDHLKRGRKKGRNRQEVPAFTEPDTTAEVKTPFDLAKAQENSNSVKKKTKFVNLYTREGQDRLAVLLPGRHPCDCLGQKHKLINNCLICGRIVCEQEGSGPCLFCGTLVCTHEEQDILQRDSNKSQKLLKKLMSGVENSGKVDISTKDLLPHQELRIKSGLEKAIKHKDKLLEFDRTSIRRTQVIDDESDYFASDSNQWLSKLERETLQKREEELRELRHASRLSKKVTIDFAGRKILEEENSLAEYHSRLDETIQAIANGTLNQPLTKLDRSSEEPLGVLVNPNMYQSPPQWVDHTGAASQKKAFRSSGFGLEFNSFQHQLRIQDQEFQEGFDGGWCLSVHQPWASLLVRGIKRVEGRSWYTPHRGRLWIAATAKKPSPQEVSELQATYRLLRGKDVEFPNDYPSGCLLGCVDLIDCLSQKQFKEQFPDISQESDSPFVFICKNPQEMVVKFPIKGNPKIWKLDSKIHQGAKKGLMKQNKAV",
   "taxon_ids": [
    "9606"
   ],
   "accessions": [
    "Q15650"
   ],
   "assertions": [
    {
     "assertion_id": "f9e5e2678617ec09bbb5cc8c5662525b1e370b18f0fca7b264a8ed8819733bef:GO:0008270",
     "go_id": "GO:0008270",
     "go_name": "zinc ion binding",
     "go_definition": "\"Binding to a zinc ion (Zn).\" [GOC:ai]",
     "supporting_assertions": [
      {
       "accession": "Q15650",
       "symbol": "TRIP4",
       "relation": "enables",
       "go_id": "GO:0008270",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR009349"
       ],
       "object_name": "Activating signal cointegrator 1",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR009349/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "160379d844fb72682f97b572d2344bb12457e6d0f4dd744e981ebe92a3b65397",
   "sequence": "MVTLDKTATKIIKLNSGAIPAIQSVLEQEAMVLLGLATTYVLVVNGASPRAVAQLRAFKQFDKDQPLGILTRGDRAMDVAKIPPAAKKLMGLFPLPITLIVEALPKLDPGITQGFRNIFIACPDAFVFDLVGAMPFPLVCATAKVGGEPVTSFEAASRYFEGQVPLIGDGGRCRYARRGTLIDCTLPLPTIMNFGPISFDDLRPLIPEIILPSHLMK",
   "taxon_ids": [],
   "accessions": [
    "Q55667"
   ],
   "assertions": [
    {
     "assertion_id": "160379d844fb72682f97b572d2344bb12457e6d0f4dd744e981ebe92a3b65397:GO:0000049",
     "go_id": "GO:0000049",
     "go_name": "tRNA binding",
     "go_definition": "\"Binding to a transfer RNA.\" [GOC:ai]",
     "supporting_assertions": [
      {
       "accession": "Q55667",
       "symbol": "slr0006",
       "relation": "enables",
       "go_id": "GO:0000049",
       "references": [
        "GO_REF:0000033"
       ],
       "evidence_code": "IBA",
       "with_from": [
        "PANTHER:PTN000426057",
        "UniProtKB:P45748"
       ],
       "object_name": "L-threonylcarbamoyladenylate synthase",
       "taxon": [
        "taxon:1111708"
       ],
       "annotation_date": "20170228",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000033.html",
        "https://www.uniprot.org/uniprotkb/P45748/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "36271c1ed723dc07d3fbf66c1fddf5147a0be57e372ae451cacc79d01c849197",
   "sequence": "MSEKTFLVEIGTEELPPKALRSLAESFAANFTAELDNAGLAHGTVQWFAAPRRLALKVANLAEAQPDREIEKRGPAIAQAFDAEGKPSKAAEGWARGCGITVDQAERLTTDKGEWLLYRAHVKGESTEALLPNMVATSLAKLPIPKLMRWGASDVHFVRPVHTVTLLLGDKVIPATILGIQSDRVIRGHRFMGEPEFTIDNADQYPEILRERGKVIADYEERKAKIKADAEEAARKIGGNADLSESLLEEVASLVEWPVVLTAKFEEKFLAVPAEALVYTMKGDQKYFPVYANDGKLLPNFIFVANIESKDPQQIISGNEKVVRPRLADAEFFFNTDRKKRLEDNLPRLQTVLFQQQLGTLRDKTDRIQALAGWIAEQIGADVNHATRAGLLSKCDLMTNMVFEFTDTQGVMGMHYARHDGEAEDVAVALNEQYQPRFAGDDLPSNPVACALAIADKMDTLAGIFGIGQHPKGDKDPFALRRAALGVLRIIVEKNLNLDLQTLTEEAVRLYGDKLTNANVVDDVIDFMLGRFRAWYQDEGYTVDTIQAVLARRPTRPADFDARMKAVSHFRTLDAAAALAAANKRVSNILAKSDEVLSDRVNASTLKEPEEIKLAMQVVVLRDKLEPYFTEGRYQDALVELAELREPVDAFFDKVMVMVDDKELRINRLTMLEKLRELFLRVADISLLQ",
   "taxon_ids": [
    "83333"
   ],
   "accessions": [
    "A8A616",
    "B1IZN4",
    "B1X8H5",
    "C4ZXE9",
    "P00961"
   ],
   "assertions": [
    {
     "assertion_id": "36271c1ed723dc07d3fbf66c1fddf5147a0be57e372ae451cacc79d01c849197:GO:0005524",
     "go_id": "GO:0005524",
     "go_name": "ATP binding",
     "go_definition": "\"Binding to ATP, adenosine 5'-triphosphate, a universally important coenzyme and enzyme regulator.\" [ISBN:0198506732]",
     "supporting_assertions": [
      {
       "accession": "B1IZN4",
       "symbol": "glyS",
       "relation": "enables",
       "go_id": "GO:0005524",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR006194",
        "InterPro:IPR008909",
        "InterPro:IPR015944"
       ],
       "object_name": "Glycine--tRNA ligase beta subunit",
       "taxon": [
        "taxon:481805"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR006194/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR008909/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR015944/"
       ]
      },
      {
       "accession": "B1IZN4",
       "symbol": "glyS",
       "relation": "enables",
       "go_id": "GO:0005524",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR000034380"
       ],
       "object_name": "Glycine--tRNA ligase beta subunit",
       "taxon": [
        "taxon:481805"
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
       "accession": "A8A616",
       "symbol": "glyS",
       "relation": "enables",
       "go_id": "GO:0005524",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR006194",
        "InterPro:IPR008909",
        "InterPro:IPR015944"
       ],
       "object_name": "Glycine--tRNA ligase beta subunit",
       "taxon": [
        "taxon:331112"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR006194/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR008909/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR015944/"
       ]
      },
      {
       "accession": "A8A616",
       "symbol": "glyS",
       "relation": "enables",
       "go_id": "GO:0005524",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR000034380"
       ],
       "object_name": "Glycine--tRNA ligase beta subunit",
       "taxon": [
        "taxon:331112"
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
       "accession": "P00961",
       "symbol": "glyS",
       "relation": "enables",
       "go_id": "GO:0005524",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR006194",
        "InterPro:IPR008909",
        "InterPro:IPR015944"
       ],
       "object_name": "Glycine--tRNA ligase beta subunit",
       "taxon": [
        "taxon:83333"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR006194/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR008909/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR015944/"
       ]
      },
      {
       "accession": "P00961",
       "symbol": "glyS",
       "relation": "enables",
       "go_id": "GO:0005524",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR000034380"
       ],
       "object_name": "Glycine--tRNA ligase beta subunit",
       "taxon": [
        "taxon:83333"
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
       "accession": "C4ZXE9",
       "symbol": "glyS",
       "relation": "enables",
       "go_id": "GO:0005524",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR006194",
        "InterPro:IPR008909",
        "InterPro:IPR015944"
       ],
       "object_name": "Glycine--tRNA ligase beta subunit",
       "taxon": [
        "taxon:595496"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR006194/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR008909/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR015944/"
       ]
      },
      {
       "accession": "C4ZXE9",
       "symbol": "glyS",
       "relation": "enables",
       "go_id": "GO:0005524",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR000034380"
       ],
       "object_name": "Glycine--tRNA ligase beta subunit",
       "taxon": [
        "taxon:595496"
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
       "accession": "B1X8H5",
       "symbol": "glyS",
       "relation": "enables",
       "go_id": "GO:0005524",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR006194",
        "InterPro:IPR008909",
        "InterPro:IPR015944"
       ],
       "object_name": "Glycine--tRNA ligase beta subunit",
       "taxon": [
        "taxon:316385"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR006194/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR008909/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR015944/"
       ]
      },
      {
       "accession": "B1X8H5",
       "symbol": "glyS",
       "relation": "enables",
       "go_id": "GO:0005524",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR000034380"
       ],
       "object_name": "Glycine--tRNA ligase beta subunit",
       "taxon": [
        "taxon:316385"
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
   "sequence_sha256": "40d3c32c3a19c7ff0d58987f1624271cf6d358c783e4214c4a15552afc9f436b",
   "sequence": "MLRLVLVALVVAVSSADFTVLKSPQYQIFQEGKWPVPGDRIPDVISLSMGFSVEEGLDWPGLGVGNLFQRPRSTVLVTVTGVNKLPLSGNGMSYPVENAVPYSLDSVVNSVHSVFSEEMPVVLQLAPIEERVYMVGKANTVFEDLSVTLRQLRTRLEQDNSIIQSLPVSSLYKNDETDRLFLSELQVLQDIITLLSGHKHLAKDNVPDIYSLELTGLEEIKKRYGEESAKFKDAVKILSDFLQKFADDMYSIYGGNAIVEVVAVESFEIPLVRRSRSILASEAISNPGSPYNLAYQYNFDYAVIFNIILWIMIGLALAVIAISYSLWNMDPGYDSIIYRMTNQKIRMD",
   "taxon_ids": [
    "8355"
   ],
   "accessions": [
    "A0A8J0PRR3",
    "Q7T0S3"
   ],
   "assertions": [
    {
     "assertion_id": "40d3c32c3a19c7ff0d58987f1624271cf6d358c783e4214c4a15552afc9f436b:GO:0038023",
     "go_id": "GO:0038023",
     "go_name": "signaling receptor activity",
     "go_definition": "\"Receiving a signal and transmitting it in the cell to initiate a change in cell activity. A signal is a physical entity or change in state that is used to transfer information in order to trigger a response.\" [GOC:bf, GOC:signaling]",
     "supporting_assertions": [
      {
       "accession": "A0A8J0PRR3",
       "symbol": "atp6ap2.S",
       "relation": "enables",
       "go_id": "GO:0038023",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR012493"
       ],
       "object_name": "Renin receptor",
       "taxon": [
        "taxon:8355"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR012493/"
       ]
      },
      {
       "accession": "Q7T0S3",
       "symbol": "atp6ap2.S",
       "relation": "enables",
       "go_id": "GO:0038023",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR012493"
       ],
       "object_name": "Renin receptor",
       "taxon": [
        "taxon:8355"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR012493/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "f1f44589619ff2613cf22d7cd06d7cd9cde096160cf6ab3e783b657ecd74299d",
   "sequence": "RAEPGRDRVGGSMNSIKNVPARVLSRRPGHSLEAEREQFDKTQAISISKAINSQEAPVKEKHARRIILGTHHEKGAFTFWSYAIGLPLPSSSILSWKFCHVLHKVLRDGHPNVLHDCQRYRSNIREIGDLWGHLRDQYGHLVNIYTKLLLTKISFHLKHPQFPAGLEVTDEVLEKAAGTDVNNIFQLTVEMFDYMDCELKLSESVFRQLNTAIAVSQMSSGQCRLAPLIQVIQDCSHLYHYTVKLMFKLHSCLPADTLQGHRDRFHEQFHSLRNFFRRASDMLYFKRLIQIPRLPEGPPNFLRASALAEHIKPVVVIPEEAPEEEEPENLIEINSGPPAGEPVQVVADLFDQTFGPPNGSMKDDRDLQIENLKREVETLRAELEKIKMEAQRYISQLKGQVNSLEAELEEQRKQKQKALVDNEQLRHELAQLKALQLEGARNQGLREEAERKASATEARYSKLKEKHNELINTHAELLRKNADTAKQLTVTAEQEEVARVKEQLAFQMEQVKRESEMKMEEQSDQLEKLKRELVAKAGELAHAQEALSRTEQSGSELSSRLDTLNAEKEALSGAIRQREAELLAAQSLVREKEEALSQEQQRAQEKGELQGRLAEKESQEQGLQQKLLDEQFAVLRGAAAEAEAILQDAVSKLDDPLHLRCTSSPDYLVSRAQAALDSVSGLEKGHTQYLASSEDASALVAALTRFSHLAADTIVNGGATSHLAPTDPADRLIDTCRECGARALELMGQLQDQTKLPRAQPSLMRAPLQGILQLGQDLKPKSLDVRQEELGAMVDKEMAATSAAIEDAVRRIEDMMNQARHESSGVKLEVNERILNSCTDLMKAIRPLVMTSTSLQKEIVESGRGAATQQEFYAKNSRWTEGLISASKAVGWGATQLVESADKVVLHMGKYEELIVCSHEIAASTAQLVAASKVKANKNSPHLSRLQECSRTVNERAANVVASTKSGQEQIEDRDTMDFSGLSLIKLKKQEMETQVRVLELEKTLEAERVRLGELRKQHYVLAGGMGTPGEEEPSRPSPAPRSGATKKPPLAQKPSIAPRPDNQLDKKDGVYPAQLVNY",
   "taxon_ids": [
    "10116"
   ],
   "accessions": [
    "Q99PW9"
   ],
   "assertions": [
    {
     "assertion_id": "f1f44589619ff2613cf22d7cd06d7cd9cde096160cf6ab3e783b657ecd74299d:GO:0030276",
     "go_id": "GO:0030276",
     "go_name": "clathrin binding",
     "go_definition": "\"Binding to a clathrin heavy or light chain, the main components of the coat of coated vesicles and coated pits, and which also occurs in synaptic vesicles.\" [GOC:jl, GOC:mah, ISBN:0198506732]",
     "supporting_assertions": [
      {
       "accession": "Q99PW9",
       "symbol": "Hip1r",
       "relation": "enables",
       "go_id": "GO:0030276",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR030224"
       ],
       "object_name": "Huntingtin-interacting protein 1-related protein",
       "taxon": [
        "taxon:10116"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR030224/"
       ]
      },
      {
       "accession": "Q99PW9",
       "symbol": "Hip1r",
       "relation": "enables",
       "go_id": "GO:0030276",
       "references": [
        "GO_REF:0000117"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "ARBA:ARBA00085651"
       ],
       "object_name": "Huntingtin-interacting protein 1-related protein",
       "taxon": [
        "taxon:10116"
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
   "sequence_sha256": "0d36b6183fa580e29232bee841a76ddc989e92f769673f7742a24ac6f35f808b",
   "sequence": "METSVSEIQVETKDEKGPVAASPQKERQERKTATLCFKRRKKANKTKPKAGSRTAEETKKHTPEAGGSGQRQPAGAWASIKGLVTHRKRSEPAKKQKPPEAEVQPEDGALPKKKAKSRLKFPCLRFSRGAKRSRHSKLTEDSGYVRVQGEADDLEIKAQTQPDDQAIQAGSTQGLQEGVLVRDGKKSQESHISNSVTSGENVIAIELELENKSSAIQMGTPELEKETKVITEKPSVQTQRASLLESSAAGSPRSVTSAAPPSPATTHQHSLEEPSNGIRESAPSGKDDRRKTAAEEKKSGETALGQAEEAAVGQADKRALSQAGEATAGHPEEATVIQAESQAKEGKLSQAEETTVAQAKETVLSQAKEGELSQAKKATVGQAEEATIDHTEKVTVDQAEETTVGQAEEATVGQAGEAILSQAKEATVVGQAEEATVDRAEEATVGQAEEATVGHTEKVTVDQAEEATVGQAEEATVGQAEEATVDWAEKPTVGQAEEATVGQAEEATVGHTEKVTVDQAEEATVGQAEEATVGHTEKVTVDHAEEATVGQAEEATVGQAEKVTVDHAEEATVGQAEEATVGQAEKVTVDHAEEATVGQAEEATVGQAEKVTVDQAEEPTVDQAEEAISSHAPDLKENGIDTEKPRSEESKRMEPIAIIITDTEISEFDVKKSKNVPKQFLISMENEQVGVFANDSDFEGRTSEQYETLLIETASSLVKNAIELSVEQLVNEMVSEDNQINTLFQ",
   "taxon_ids": [
    "10090"
   ],
   "accessions": [
    "D3YVF0"
   ],
   "assertions": [
    {
     "assertion_id": "0d36b6183fa580e29232bee841a76ddc989e92f769673f7742a24ac6f35f808b:GO:0044877",
     "go_id": "GO:0044877",
     "go_name": "protein-containing complex binding",
     "go_definition": "\"Binding to a macromolecular complex.\" [GOC:jl]",
     "supporting_assertions": [
      {
       "accession": "D3YVF0",
       "symbol": "Akap5",
       "relation": "enables",
       "go_id": "GO:0044877",
       "references": [
        "GO_REF:0000096"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "RGD:620829"
       ],
       "object_name": "A-kinase anchor protein 5",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20250828",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000096.html"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "16c31ae9f730d5497e2102a3945750dcb9f2e8273f4771dbc85f5a03b9bdeb0f",
   "sequence": "MLLCKAFWEMHASRPGAPPHAGLPSRRTAGDQDHRADPKVTRIMSASTLEQPAAAHVDELVARMRGRLLDPLAIAVLAAVISGAWASRPSLWFDEGATISASASRTLPELWSLLGHIDAVHGLYYLLMHGWFAIFPPTELWSRLPSCLAIGAAAAGVVVFAKQFSGRTTAVCAGAVFAILPRVTWAGIEARSSALSVAAAVWLTVLLVAAVRCNTQRRWLLYALVLMLSILVSINLALLVPAYATMVPLLASGKSRKSPVIWWTVVTAAALGAMTPFILFAHGQVWQVGWIAGLNRNIILDVIHRQYFDHSVPFAILAGLIVAAGIAAHLAGARGPGGDTHRLVLVSAAWIVVPTAVVLIYSATVEPIYYPRYLILTAPAAAVILAVCVVTIARKPWLIAGVVFLLAAAAFPNYFFTQRGPYAKEGWDYSQVADVISAHAKPGDCLLVDNTAGWRPGPIRALLATRPAAFRSLIDVERGTYGPKVGTLWDGHVAVWLTTAKIDKCPTLWTIANRDKSLPDHQVGEMLSPGTGFGRTPVYRFPSYLGFRIVERWQFHYSQVVKSTR",
   "taxon_ids": [
    "83331"
   ],
   "accessions": [
    "Q7D885"
   ],
   "assertions": [
    {
     "assertion_id": "16c31ae9f730d5497e2102a3945750dcb9f2e8273f4771dbc85f5a03b9bdeb0f:GO:0016763",
     "go_id": "GO:0016763",
     "go_name": "pentosyltransferase activity",
     "go_definition": "\"Catalysis of the transfer of a pentosyl group from one compound (donor) to another (acceptor).\" [GOC:jl]",
     "supporting_assertions": [
      {
       "accession": "Q7D885",
       "symbol": "MT1671",
       "relation": "enables",
       "go_id": "GO:0016763",
       "references": [
        "GO_REF:0000118"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "PANTHER:PTN005029416"
       ],
       "object_name": "Uncharacterized protein",
       "taxon": [
        "taxon:83331"
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
   "sequence_sha256": "5cb755fefb6323f1b454ec92533fa7e1035f6b2f266b9f4720a023ed4be47471",
   "sequence": "MTMESGADNQQSGDAAVTEAENQQMTVQAQPQIATLAQVSMPAAHATSSAPTVTLVQLPNGQTVQVHGVIQAAQPSVIQSPQVQTVQISTIAESEDSQESVDSVTDSQKRREILSRRPSYRKILNDLSSDAPGVPRIEEEKSEEETSAPAITTVTVPTPIYQTSSGQYIAITQGGAIQLANNGTDGVQGLQTLTMTNAAATQPGTTILQYAQTTDGQQILVPSNQVVVQAASGDVQTYQIRTAPTSTIAPGVVMASSPALPTQPAEEAARKREVRLMKNREAARECRRKKKEYVKCLENRVAVLENQNKTLIEELKALKDLYCHKSD",
   "taxon_ids": [
    "10090"
   ],
   "accessions": [
    "Q01147",
    "Q543W0"
   ],
   "assertions": [
    {
     "assertion_id": "5cb755fefb6323f1b454ec92533fa7e1035f6b2f266b9f4720a023ed4be47471:GO:0042802",
     "go_id": "GO:0042802",
     "go_name": "identical protein binding",
     "go_definition": "\"Binding to an identical protein or proteins.\" [GOC:jl]",
     "supporting_assertions": [
      {
       "accession": "Q543W0",
       "symbol": "Creb1",
       "relation": "enables",
       "go_id": "GO:0042802",
       "references": [
        "GO_REF:0000107"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB:P16220",
        "ensembl:ENSP00000236995"
       ],
       "object_name": "Cyclic AMP-responsive element-binding protein 1",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20260724",
       "assigned_by": "Ensembl",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000107.html",
        "https://www.uniprot.org/uniprotkb/P16220/entry"
       ]
      },
      {
       "accession": "Q01147",
       "symbol": "Creb1",
       "relation": "enables",
       "go_id": "GO:0042802",
       "references": [
        "GO_REF:0000119"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "UniProtKB:P16220"
       ],
       "object_name": "Cyclic AMP-responsive element-binding protein 1",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20260723",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000119.html",
        "https://www.uniprot.org/uniprotkb/P16220/entry"
       ]
      },
      {
       "accession": "Q01147",
       "symbol": "Creb1",
       "relation": "enables",
       "go_id": "GO:0042802",
       "references": [
        "GO_REF:0000107"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB:P16220",
        "ensembl:ENSP00000236995"
       ],
       "object_name": "Cyclic AMP-responsive element-binding protein 1",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20260724",
       "assigned_by": "Ensembl",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000107.html",
        "https://www.uniprot.org/uniprotkb/P16220/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "90403d58e3a0bd9abd9be8379acae41076e0a6c93bb133d7366ca88085204929",
   "sequence": "MHQGCAAVLCVLAVLGLEVAALGECELARLLQDKLRYEMRLQYMKHNFPIDYTLRVQHEEVLRTANVTRLRDGKVSEASLRYLWFHACSQAVLHILEVLPEKHPSRGYTQELSQLLDALGVEYSGYRQSDVDAVVADLVKQLHSGDSRQKAVRPKALLDNCLKVLRMLFGAHCRWDSA",
   "taxon_ids": [
    "9031"
   ],
   "accessions": [
    "R4GIU3"
   ],
   "assertions": [
    {
     "assertion_id": "90403d58e3a0bd9abd9be8379acae41076e0a6c93bb133d7366ca88085204929:GO:0005157",
     "go_id": "GO:0005157",
     "go_name": "macrophage colony-stimulating factor receptor binding",
     "go_definition": "\"Binding to a macrophage colony-stimulating factor receptor.\" [GOC:ai]",
     "supporting_assertions": [
      {
       "accession": "R4GIU3",
       "symbol": "IL34",
       "relation": "enables",
       "go_id": "GO:0005157",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR020415"
       ],
       "object_name": "Interleukin 34",
       "taxon": [
        "taxon:9031"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR020415/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "73c056b6308e57e3a0d3d3a8558ad0859f2f58cbe233d933ec531a027ae48293",
   "sequence": "MVMSQGTYTFLTCFAGFWLIWGLIVLLCCFCSFLRRRLKRRQEERLREQNLRALELEPLELEGSLAGSPPGLAPPQPPPHRSRLEAPAHAHSHPHVHVHPLLHHGPAQPHAHAHPHPHHHALPHPPPTHLSVPPRPWSYPRQAESDMSKPPCYEEAVLMAEPPPPYSEVLTDTRGLYRKIVTPFLSRRDSAEKQEQPPPSYKPLFLDRGYTSALHLPSAPRPAPPCPALCLQADRGRRVFPSWTDSELSSREPLEHGAWRLPVSIPLFGRTTAV",
   "taxon_ids": [
    "9606"
   ],
   "accessions": [
    "Q8TB68"
   ],
   "assertions": [
    {
     "assertion_id": "73c056b6308e57e3a0d3d3a8558ad0859f2f58cbe233d933ec531a027ae48293:GO:0044877",
     "go_id": "GO:0044877",
     "go_name": "protein-containing complex binding",
     "go_definition": "\"Binding to a macromolecular complex.\" [GOC:jl]",
     "supporting_assertions": [
      {
       "accession": "Q8TB68",
       "symbol": "PRR7",
       "relation": "enables",
       "go_id": "GO:0044877",
       "references": [
        "GO_REF:0000107"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB:P0C6T3",
        "ensembl:ENSRNOP00000038274"
       ],
       "object_name": "Proline-rich protein 7",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20260724",
       "assigned_by": "Ensembl",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000107.html",
        "https://www.uniprot.org/uniprotkb/P0C6T3/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "bb884f2939cb91546f53f9471fb098e3d926cea41120dba7b7798efc23f0991d",
   "sequence": "MSGSVWNSNSWHWEERNYNKWAESYIKYNLSNLKIEKEDLTIYFDNLQVSGNACVSIRKGKQINSFEYIIKFEWLYSKKKEGKDYFGGSVEIPDFSTFSLEENDYAINIERTDESENLRFIYDSILKKEGKEKIKECLKNFQEDLLKHDKNESNKELKIKEEEKIKLENIKTCNEKKTEEENKSNQNINNNINDEKKEGSVWNINNYHWEEKCLTKWAIEELQNIFNKSIIELSNNIFLEFFSCDVEGEASSSLRKKKKILMYDLKITSEWKAYQKNKNQQIEIESKGHVSINDILSDFSSDDNTKYSYYFIFDNKTDEYNQINDVIKLEGPNKINQIIDDFILKMREK",
   "taxon_ids": [
    "36329"
   ],
   "accessions": [
    "O97244",
    "Q689C6"
   ],
   "assertions": [
    {
     "assertion_id": "bb884f2939cb91546f53f9471fb098e3d926cea41120dba7b7798efc23f0991d:GO:0051087",
     "go_id": "GO:0051087",
     "go_name": "protein-folding chaperone binding",
     "go_definition": "\"Binding to a chaperone protein, a class of proteins that bind to nascent or unfolded polypeptides and ensure correct folding or transport.\" [PMID:10585443]",
     "supporting_assertions": [
      {
       "accession": "O97244",
       "symbol": "PF3D7_0306200",
       "relation": "enables",
       "go_id": "GO:0051087",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR015310"
       ],
       "object_name": "Activator of Hsp90 ATPase AHSA1-like N-terminal domain-containing protein",
       "taxon": [
        "taxon:36329"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR015310/"
       ]
      },
      {
       "accession": "Q689C6",
       "symbol": "PFC0270w",
       "relation": "enables",
       "go_id": "GO:0051087",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR015310"
       ],
       "object_name": "Activator of Hsp90 ATPase AHSA1-like N-terminal domain-containing protein",
       "taxon": [
        "taxon:5833"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR015310/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "f193c39fb35c917f19479c729f06f6a62a890f04ab9d07aa372882f3705210a1",
   "sequence": "MSLLSAIDTSAASVYQPAQLLNWVYLSLQDTHQASAFDAFRPEPTAGAAPPELAFGKGRPEQLGSPLHSSYLNSFFQLQRGEALSNSVYKGASPYGSLNNIADGLSSLTEHFSDLTLTSEARKPSKRPPPNYLCHLCFNKGHYIKDCPQARPKGEGLTPYQGKKRCFGEYKCPKCKRKWMSGNSWANMGQECIKCHINVYPHKQRPLEKPDGLDVSDQSKEHPQHLCEKCKVLGYYCRRVQ",
   "taxon_ids": [
    "9606"
   ],
   "accessions": [
    "A0A0D9RAE2",
    "A0A2I2YK31",
    "A0A2K5CAQ0",
    "A0A2K5N0Z8",
    "A0A2K5RHT2",
    "A0A2K6DQ18",
    "A0A2K6LJF0",
    "A0A2K6THT2",
    "A0A5F8AIR8",
    "A0A663D9F1",
    "A0A6J3GAS2",
    "A0A7N9CUK3",
    "A0A8C9IT94",
    "A0A8D2FFR9",
    "A0A8I5MW21",
    "A0ACF8HS03",
    "A0ACF8KDD8",
    "A0ACF8PUK9",
    "F7D1P2",
    "H2Q254",
    "K7CLX7",
    "Q8N2G6"
   ],
   "assertions": [
    {
     "assertion_id": "f193c39fb35c917f19479c729f06f6a62a890f04ab9d07aa372882f3705210a1:GO:0008270",
     "go_id": "GO:0008270",
     "go_name": "zinc ion binding",
     "go_definition": "\"Binding to a zinc ion (Zn).\" [GOC:ai]",
     "supporting_assertions": [
      {
       "accession": "A0A2K5CAQ0",
       "symbol": "A0A2K5CAQ0",
       "relation": "enables",
       "go_id": "GO:0008270",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001878",
        "InterPro:IPR036875"
       ],
       "object_name": "CCHC-type domain-containing protein",
       "taxon": [
        "taxon:37293"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001878/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036875/"
       ]
      },
      {
       "accession": "A0A663D9F1",
       "symbol": "ZCCHC24",
       "relation": "enables",
       "go_id": "GO:0008270",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001878",
        "InterPro:IPR036875"
       ],
       "object_name": "CCHC-type domain-containing protein",
       "taxon": [
        "taxon:9601"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001878/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036875/"
       ]
      },
      {
       "accession": "A0A8C9IT94",
       "symbol": "ZCCHC24",
       "relation": "enables",
       "go_id": "GO:0008270",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001878",
        "InterPro:IPR036875"
       ],
       "object_name": "CCHC-type domain-containing protein",
       "taxon": [
        "taxon:591936"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001878/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036875/"
       ]
      },
      {
       "accession": "A0A2K6DQ18",
       "symbol": "A0A2K6DQ18",
       "relation": "enables",
       "go_id": "GO:0008270",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001878",
        "InterPro:IPR036875"
       ],
       "object_name": "CCHC-type domain-containing protein",
       "taxon": [
        "taxon:9545"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001878/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036875/"
       ]
      },
      {
       "accession": "A0ACF8KDD8",
       "symbol": "H8959_011081",
       "relation": "enables",
       "go_id": "GO:0008270",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001878",
        "InterPro:IPR036875"
       ],
       "object_name": "CCHC-type domain-containing protein",
       "taxon": [
        "taxon:310352"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001878/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036875/"
       ]
      },
      {
       "accession": "A0A0D9RAE2",
       "symbol": "A0A0D9RAE2",
       "relation": "enables",
       "go_id": "GO:0008270",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001878",
        "InterPro:IPR036875"
       ],
       "object_name": "CCHC-type domain-containing protein",
       "taxon": [
        "taxon:60711"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001878/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036875/"
       ]
      },
      {
       "accession": "A0A2K6THT2",
       "symbol": "A0A2K6THT2",
       "relation": "enables",
       "go_id": "GO:0008270",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001878",
        "InterPro:IPR036875"
       ],
       "object_name": "CCHC-type domain-containing protein",
       "taxon": [
        "taxon:39432"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001878/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036875/"
       ]
      },
      {
       "accession": "H2Q254",
       "symbol": "ZCCHC24",
       "relation": "enables",
       "go_id": "GO:0008270",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001878",
        "InterPro:IPR036875"
       ],
       "object_name": "CCHC-type domain-containing protein",
       "taxon": [
        "taxon:9598"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001878/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036875/"
       ]
      },
      {
       "accession": "K7CLX7",
       "symbol": "ZCCHC24",
       "relation": "enables",
       "go_id": "GO:0008270",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001878",
        "InterPro:IPR036875"
       ],
       "object_name": "CCHC-type domain-containing protein",
       "taxon": [
        "taxon:9598"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001878/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036875/"
       ]
      },
      {
       "accession": "A0A2I2YK31",
       "symbol": "ZCCHC24",
       "relation": "enables",
       "go_id": "GO:0008270",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001878",
        "InterPro:IPR036875"
       ],
       "object_name": "CCHC-type domain-containing protein",
       "taxon": [
        "taxon:9595"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001878/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036875/"
       ]
      },
      {
       "accession": "A0A6J3GAS2",
       "symbol": "ZCCHC24",
       "relation": "enables",
       "go_id": "GO:0008270",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001878",
        "InterPro:IPR036875"
       ],
       "object_name": "CCHC-type domain-containing protein",
       "taxon": [
        "taxon:9515"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001878/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036875/"
       ]
      },
      {
       "accession": "A0A2K6LJF0",
       "symbol": "A0A2K6LJF0",
       "relation": "enables",
       "go_id": "GO:0008270",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001878",
        "InterPro:IPR036875"
       ],
       "object_name": "CCHC-type domain-containing protein",
       "taxon": [
        "taxon:61621"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001878/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036875/"
       ]
      },
      {
       "accession": "A0A8D2FFR9",
       "symbol": "ZCCHC24",
       "relation": "enables",
       "go_id": "GO:0008270",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001878",
        "InterPro:IPR036875"
       ],
       "object_name": "CCHC-type domain-containing protein",
       "taxon": [
        "taxon:9565"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001878/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036875/"
       ]
      },
      {
       "accession": "F7D1P2",
       "symbol": "ZCCHC24",
       "relation": "enables",
       "go_id": "GO:0008270",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001878",
        "InterPro:IPR036875"
       ],
       "object_name": "CCHC-type domain-containing protein",
       "taxon": [
        "taxon:9544"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001878/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036875/"
       ]
      },
      {
       "accession": "A0ACF8PUK9",
       "symbol": "H8958_022462",
       "relation": "enables",
       "go_id": "GO:0008270",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001878",
        "InterPro:IPR036875"
       ],
       "object_name": "CCHC-type domain-containing protein",
       "taxon": [
        "taxon:43780"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001878/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036875/"
       ]
      },
      {
       "accession": "A0A7N9CUK3",
       "symbol": "A0A7N9CUK3",
       "relation": "enables",
       "go_id": "GO:0008270",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001878",
        "InterPro:IPR036875"
       ],
       "object_name": "CCHC-type domain-containing protein",
       "taxon": [
        "taxon:9541"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001878/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036875/"
       ]
      },
      {
       "accession": "A0A5F8AIR8",
       "symbol": "ZCCHC24",
       "relation": "enables",
       "go_id": "GO:0008270",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001878",
        "InterPro:IPR036875"
       ],
       "object_name": "CCHC-type domain-containing protein",
       "taxon": [
        "taxon:9544"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001878/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036875/"
       ]
      },
      {
       "accession": "A0A2K5RHT2",
       "symbol": "ZCCHC24",
       "relation": "enables",
       "go_id": "GO:0008270",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001878",
        "InterPro:IPR036875"
       ],
       "object_name": "CCHC-type domain-containing protein",
       "taxon": [
        "taxon:2715852"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001878/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036875/"
       ]
      },
      {
       "accession": "Q8N2G6",
       "symbol": "ZCCHC24",
       "relation": "enables",
       "go_id": "GO:0008270",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001878",
        "InterPro:IPR036875"
       ],
       "object_name": "Zinc finger CCHC domain-containing protein 24",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001878/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036875/"
       ]
      },
      {
       "accession": "A0A2K5N0Z8",
       "symbol": "A0A2K5N0Z8",
       "relation": "enables",
       "go_id": "GO:0008270",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001878",
        "InterPro:IPR036875"
       ],
       "object_name": "CCHC-type domain-containing protein",
       "taxon": [
        "taxon:9531"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001878/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036875/"
       ]
      },
      {
       "accession": "A0A8I5MW21",
       "symbol": "A0A8I5MW21",
       "relation": "enables",
       "go_id": "GO:0008270",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001878",
        "InterPro:IPR036875"
       ],
       "object_name": "CCHC-type domain-containing protein",
       "taxon": [
        "taxon:9555"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001878/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036875/"
       ]
      },
      {
       "accession": "A0ACF8HS03",
       "symbol": "H8957_011638",
       "relation": "enables",
       "go_id": "GO:0008270",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001878",
        "InterPro:IPR036875"
       ],
       "object_name": "CCHC-type domain-containing protein",
       "taxon": [
        "taxon:88029"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001878/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036875/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "dd921fb39ec53cbc5d6b357b28437456c35b1a2f7abf386d32083c78095d6277",
   "sequence": "MYNGIGLPTPRGSGTNGYVQRNLSLVRGRRGERPDYKGEEELRRLEAALVKRPNPDILDHERKRRVELRCLELEEMMEEQGYEEQQIQEKVATFRLMLLEKDVNPGGKEETPGQRPAVTETHQLAELNEKKNERLRAAFGISDSYVDGSSFDPQRRAREAKQPAPEPPKPYSLVRESSSSRSPTPKQKKKKKKKDRGRRSESSSPRRERKKSSKKKKHRSESESKKRKHRSPTPKSKRKSKDKKRKRSRSTTPAPKSRRAHRSTSADSASSSDTSRSRSRSAAAKTHTTALAGRSPSPASGRRGEGDAPFSEPGTTSTQRPSSPETATKQPSSPYEDKDKDKKEKSATRPSPSPERSSTGPEPPAPTPLLAERHGGSPQPLATTPLSQEPVNPPSEASPTRDRSPPKSPEKLPQSSSSESSPPSPQPTKVSRHASSSPESPKPAPAPGSHREISSSPTSKNRSHGRAKRDKSHSHTPSRRMGRSRSPATAKRGRSRSRTPTKRGHSRSRSPQWRRSRSAQRWGRSRSPQRRGRSRSPQRPGWSRSRNTQRRGRSRSARRGRSHSRSPATRGRSRSRTPARRGRSRSRTPARRRSRSRTPTRRRSRSRTPARRGRSRSRTPARRRSRTRSPVRRRSRSRSPARRSGRSRSRTPARRGRSRSRTPARRGRSRSRTPARRSGRSRSRTPARRGRSRSRTPRRGRSRSRSLVRRGRSHSRTPQRRGRSGSSSERKNKSRTSQRRSRSNSSPEMKKSRISSRRSRSLSSPRSKAKSRLSLRRSLSGSSPCPKQKSQTPPRRSRSGSSQPKAKSRTPPRRSRSSSSPPPKQKSKTPSRQSHSSSSPHPKVKSGTPPRQGSITSPQANEQSVTPQRRSCFESSPDPELKSRTPSRHSCSGSSPPRVKSSTPPRQSPSRSSSPQPKVKAIISPRQRSHSGSSSPSPSRVTSRTTPRRSRSVSPCSNVESRLLPRYSHSGSSSPDTKVKPETPPRQSHSGSISPYPKVKAQTPPGPSLSGSKSPCPQEKSKDSLVQSCPGSLSLCAGVKSSTPPGESYFGVSSLQLKGQSQTSPDHRSDTSSPEVRQSHSESPSLQSKSQTSPKGGRSRSSSPVTELASRSPIRQDRGEFSASPMLKSGMSPEQSRFQSDSSSYPTVDSNSLLGQSRLETAESKEKMALPPQEDATASPPRQKDKFSPFPVQDRPESSLVFKDTLRTPPRERSGAGSSPETKEQNSALPTSSQDEELMEVVEKSEEPAGQILSHLSSELKEMSTSNFESSPEVEERPAVSLTLDQSQSQASLEAVEVPSMASSWGGPHFSPEHKELSNSPLRENSFGSPLEFRNSGPLGTEMNTGFSSEVKEDLNGPFLNQLETDPSLDMKEQSTRSSGHSSSELSPDAVEKAGMSSNQSISSPVLDAVPRTPSRERSSSASSPEMKDGLPRTPSRRSRSGSSPGLRDGSGTPSRHSLSGSSPGMKDIPRTPSRGRSECDSSPEPKALPQTPRPRSRSPSSPELNNKCLTPQRERSGSESSVDQKTVARTPLGQRSRSGSSQELDVKPSASPQERSESDSSPDSKAKTRTPLRQRSRSGSSPEVDSKSRLSPRRSRSGSSPEVKDKPRAAPRAQSGSDSSPEPKAPAPRALPRRSRSGSSSKGRGPSPEGSSSTESSPEHPPKSRTARRGSRSSPEPKTKSRTPPRRRSSRSSPELTRKARLSRRSRSASSSPETRSRTPPRHRRSPSVSSPEPAEKSRSSRRRRSASSPRTKTTSRRGRSPSPKPRGLQRSRSRSRREKTRTTRRRDRSGSSQSTSRRRQRSRSRSRVTRRRRGGSGYHSRSPARQESSRTSSRRRRGRSRTPPTSRKRSRSRTSPAPWKRSRSRASPATHRRSRSRTPLISRRRSRSRTSPVSRRRSRSRTSVTRRRSRSRASPVSRRRSRSRTPPVTRRRSRSRTPTTRRRSRSRTPPVTRRRSRSRTPPVTRRRSRSRTSPITRRRSRSRTSPVTRRRSRSRTSPVTRRRSRSRTSPVTRRRSRSRTPPAIRRRSRSRTPLLPRKRSRSRSPLAIRRRSRSRTPRTARGKRSLTRSPPAIRRRSASGSSSDRSRSATPPATRNHSGSRTPPVALNSSRMSCFSRPSMSPTPLDRCRSPGMLEPLGSSRTPMSVLQQAGGSMMDGPGPRIPDHQRTSVPENHAQSRIALALTAISLGTARPPPSMSAAGLAARMSQVPAPVPLMSLRTAPAANLASRIPAASAAAMNLASARTPAIPTAVNLADSRTPAAAAAMNLASPRTAVAPSAVNLADPRTPTAPAVNLAGARTPAALAALSLTGSGTPPTAANYPSSSRTPQAPASANLVGPRSAHATAPVNIAGSRTAAALAPASLTSARMAPALSGANLTSPRVPLSAYERVSGRTSPPLLDRARSRTPPSAPSQSRMTSERAPSPSSRMGQAPSQSLLPPAQDQPRSPVPSAFSDQSRCLIAQTTPVAGSQSLSSGAVATTTSSAGDHNGMLSVPAPGVPHSDVGEPPASTGAQQPSALAALQPAKERRSSSSSSSSSSSSSSSSSSSSSSSSSGSSSSDSEGSSLPVQPEVALKRVPSPTPAPKEAVREGRPPEPTPAKRKRRSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSPSPAKPGPQALPKPASPKKPPPGERRSRSPRKPIDSLRDSRSLSYSPVERRRPSPQPSPRDQQSSSSERGSRRGQRGDSRSPSHKRRRETPSPRPMRHRSSRSP",
   "taxon_ids": [
    "9606"
   ],
   "accessions": [
    "Q9UQ35"
   ],
   "assertions": [
    {
     "assertion_id": "dd921fb39ec53cbc5d6b357b28437456c35b1a2f7abf386d32083c78095d6277:GO:0003729",
     "go_id": "GO:0003729",
     "go_name": "mRNA binding",
     "go_definition": "\"Binding to messenger RNA (mRNA), an intermediate molecule between DNA and protein. mRNA includes UTR and coding sequences, but does not contain introns.\" [GOC:kmv, GOC:pr, SO:0000234]",
     "supporting_assertions": [
      {
       "accession": "Q9UQ35",
       "symbol": "SRRM2",
       "relation": "enables",
       "go_id": "GO:0003729",
       "references": [
        "GO_REF:0000033"
       ],
       "evidence_code": "IBA",
       "with_from": [
        "MGI:MGI:1916205",
        "PANTHER:PTN002857912"
       ],
       "object_name": "Serine/arginine repetitive matrix protein 2",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20180315",
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
   "sequence_sha256": "44aaa3b89f059a76f4ff72aa851a62903a233a4c0a069e6c3f106e549632fb51",
   "sequence": "MTEMSFLSSEVLVGDLMSPFDQSGLGAEESLGLLDDYLEVAKHFKPHGFSSDKAKAGSSEWLAVDGLVSPSNNSKEDAFSGTDWMLEKMDLKEFDLDALLGIDDLETMPDDLLTTLDDTCDLFAPLVQETNKQPPQTVNPIGHLPESLTKPDQVAPFTFLQPLPLSPGVLSSTPDHSFSLELGSEVDITEGDRKPDYTAYVAMIPQCIKEEDTPSDNDSGICMSPESYLGSPQHSPSTRGSPNRSLPSPGVLCGSARPKPYDPPGEKMVAAKVKGEKLDKKLKKMEQNKTAATRYRQKKRAEQEALTGECKELEKKNEALKERADSLAKEIQYLKDLIEEVRKARGKKRVP",
   "taxon_ids": [
    "9606"
   ],
   "accessions": [
    "P18848"
   ],
   "assertions": [
    {
     "assertion_id": "44aaa3b89f059a76f4ff72aa851a62903a233a4c0a069e6c3f106e549632fb51:GO:0008140",
     "go_id": "GO:0008140",
     "go_name": "cAMP response element binding protein binding",
     "go_definition": "\"Binding to a cAMP response element binding protein (a CREB protein).\" [GOC:mah]",
     "supporting_assertions": [
      {
       "accession": "P18848",
       "symbol": "ATF4",
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
        "taxon:9606"
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
       "accession": "P18848",
       "symbol": "ATF4",
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
        "taxon:9606"
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
   "sequence_sha256": "256aa3183a862c10acf3c7d6ebe5afb53edcb1bf2d4b768308e144788d89d54a",
   "sequence": "MDKVWWHTAWGLLCLLQLSLAQQQIDLNITCRYAGVFHVEKNGRYSISRTEAADLCEAFNTTLPTMAQMELALRKGFETCRYGFIEGHVVIPRIHPNAICAANNTGVYILLASNTSHYDTYCFNASAPLEEDCTSVTDLPNSFDGPVTITIVNRDGTRYSKKGEYRTHQEDIDASNIIDEDVSSGSTIEKSTPEGYILHTDLPTSQPTGDRDDAFFIGSTLATIATTPWVSAHTKQNQERTQWNPIHSNPEVLLQTTTRMTDIDRNSTSAHGENWTQEPQPPFNNHEYQDEEETPHATSTTWADPNSTTEEAATQKEKWFENEWQGKNPPTPSEDSHVTEGTTASAHNNHPSQRMTTQSQEDVSWTDFFDPISHPMGQGHQTESKGHSSGNQDSGVTTTSGPARRPQIPEWLIILASLLALALILAVCIAVNSRRRCGQKKKLVINSGNGTVEDRKPSELNGEASKSQEMVHLVNKEPTETPDQFMTADETRNLQSVDMKIGV",
   "taxon_ids": [
    "10116"
   ],
   "accessions": [
    "P26051"
   ],
   "assertions": [
    {
     "assertion_id": "256aa3183a862c10acf3c7d6ebe5afb53edcb1bf2d4b768308e144788d89d54a:GO:0038023",
     "go_id": "GO:0038023",
     "go_name": "signaling receptor activity",
     "go_definition": "\"Receiving a signal and transmitting it in the cell to initiate a change in cell activity. A signal is a physical entity or change in state that is used to transfer information in order to trigger a response.\" [GOC:bf, GOC:signaling]",
     "supporting_assertions": [
      {
       "accession": "P26051",
       "symbol": "Cd44",
       "relation": "enables",
       "go_id": "GO:0038023",
       "references": [
        "GO_REF:0000121"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "UniProtKB:P16070"
       ],
       "object_name": "CD44 antigen",
       "taxon": [
        "taxon:10116"
       ],
       "annotation_date": "20250902",
       "assigned_by": "RGD",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000121.html",
        "https://www.uniprot.org/uniprotkb/P16070/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "200e28ef2dbd5fe21278f2c9f8ab2eddfa7f5b5d9b2660fb2331d7c601829b9d",
   "sequence": "MSAAILLAPGDVIKRSSEELKQRQIQINLVDWMESEGGKEDKTEPKEESKAEGSKDGEGTQSESGQKEEGGKETKDADVDRRIHTAVGSGSGTKGSGERANENANRGDGKVGGGGGDADAGVGATGTNGGRWVVLTEEIARAIESKYGTKIDVYRDEVPAQIIEVERSLQKELGISREGVAEQTERLRDLRRKEKNGTHAKAVERGGRKQRKKAHGDAQREGVEEEKTSEEPARIGITIEGVMSQKKLLSMIGGVERKMAPIGARESAVMLVSNSIKDVVRATAYFTAPTGDPHWKEVAREASKKKNILAYTSTGGDVKTEFLHLIDHL",
   "taxon_ids": [
    "10900"
   ],
   "accessions": [
    "A0A5B8YW99",
    "Q98829"
   ],
   "assertions": [
    {
     "assertion_id": "200e28ef2dbd5fe21278f2c9f8ab2eddfa7f5b5d9b2660fb2331d7c601829b9d:GO:0005198",
     "go_id": "GO:0005198",
     "go_name": "structural molecule activity",
     "go_definition": "\"The action of a molecule that contributes to the structural integrity of a complex.\" [GOC:mah, GOC:vw]",
     "supporting_assertions": [
      {
       "accession": "A0A5B8YW99",
       "symbol": "VP6",
       "relation": "enables",
       "go_id": "GO:0005198",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001399"
       ],
       "object_name": "VP6 protein",
       "taxon": [
        "taxon:40051"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001399/"
       ]
      },
      {
       "accession": "Q98829",
       "symbol": "Segment-9",
       "relation": "enables",
       "go_id": "GO:0005198",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001399"
       ],
       "object_name": "Helicase VP6-A",
       "taxon": [
        "taxon:10900"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001399/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "87d30f2de5a1e726c6af6db080ec6c6deb3b7672b4c82352120081b94b07e956",
   "sequence": "MNIRQSGKYYEYKTLEILEKNGFKALRIPVSGTGKQALPDLIATKNNTIYPIEVKSTSKDVVTVRNFQIEKLFKFCEIFNFCECHPLVTVYYKKYKIVIVYELSQDVRTKEKIKFKYGINS",
   "taxon_ids": [
    "157898"
   ],
   "accessions": [
    "Q98VP9"
   ],
   "assertions": [
    {
     "assertion_id": "87d30f2de5a1e726c6af6db080ec6c6deb3b7672b4c82352120081b94b07e956:GO:0003677",
     "go_id": "GO:0003677",
     "go_name": "DNA binding",
     "go_definition": "\"Any molecular function by which a gene product interacts selectively and non-covalently with DNA (deoxyribonucleic acid).\" [GOC:dph, GOC:jl, GOC:tb, GOC:vw]",
     "supporting_assertions": [
      {
       "accession": "Q98VP9",
       "symbol": "hjc",
       "relation": "enables",
       "go_id": "GO:0003677",
       "references": [
        "GO_REF:0000043"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB-KW:KW-0238"
       ],
       "object_name": "Crossover junction endodeoxyribonuclease",
       "taxon": [
        "taxon:157898"
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
   "sequence_sha256": "4b8bf0d5886abc28883f3d605eddc25c790a3548750a1b4ebe9b8a6807684025",
   "sequence": "MAPSPRPQHVLHWRDAHNFYLLSPLMGLLSRAWSRLRGPEVPEAWLAKTVTGADQIEAAALLTPTPVSGNLLPHGETEESGSPEQSQAAQRLCLVEAESSPPETWGLSNVDEYNAKPGQDDLREKEMERTAGKATLQPAGLQGADKRLGEVVAREEGVAEPAYPTSQLEGGPAENEEDGETVKTYQASAASIAPGYKPSTPVPFLGEAEHQATEEKGTENKADPSNSPSSGSHSRAWEYYSREKPKQEGEAKVEAHRAGQGHPCRNAEAEEGGPETTFVCTGNAFLKAWVYRPGEDTEEEDNSDSDSAEEDTAQTGATPHTSAFLKAWVYRPGEDTEEEDSDSDSAEEDTAQTGATPHTSAFLKAWVYRPGEDTEEENSDLDSAEEDTAQTGATPHTSAFLKAWVYRPGEDTEEENSDLDSAEEDTAQTGATPHTSPFLKAWVYRPGEDTEDDTEEEEDSENVAPGDSETADSSQSPCLQPQRCLPGEKTKGRGEEPPLFQVAFYLPGEKPESPWAAPKLPLRLQRRLRLFKAPTRDQDPEIPLKARKVHFAEKVTVHFLAVWAGPAQAARRGPWEQFARDRSRFARRIAQAEEKLGPYLTPDSRARAWARLRNPSLPQSEPRSSSEATPLTQDVTTPSPLPSETPSPSLYLGGRRG",
   "taxon_ids": [
    "10090"
   ],
   "accessions": [
    "P17564"
   ],
   "assertions": [
    {
     "assertion_id": "4b8bf0d5886abc28883f3d605eddc25c790a3548750a1b4ebe9b8a6807684025:GO:0008157",
     "go_id": "GO:0008157",
     "go_name": "protein phosphatase 1 binding",
     "go_definition": "\"Binding to a protein phosphatase 1.\" [GOC:jl]",
     "supporting_assertions": [
      {
       "accession": "P17564",
       "symbol": "Ppp1r15a",
       "relation": "enables",
       "go_id": "GO:0008157",
       "references": [
        "GO_REF:0000096"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "RGD:621526"
       ],
       "object_name": "Protein phosphatase 1 regulatory subunit 15A",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20250828",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000096.html"
       ]
      },
      {
       "accession": "P17564",
       "symbol": "Ppp1r15a",
       "relation": "enables",
       "go_id": "GO:0008157",
       "references": [
        "GO_REF:0000119"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "UniProtKB:O75807"
       ],
       "object_name": "Protein phosphatase 1 regulatory subunit 15A",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20260122",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000119.html",
        "https://www.uniprot.org/uniprotkb/O75807/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "15a1ba52355f6ad7c0ab48869591df3901d66b535d5ba0e9fa231939fb45c49e",
   "sequence": "MRVSGVLRLLALIFAIVTTWMFIRSYMSFSMKTIRLPRWLAASPTKEIQVKKYKCGLIKPCPANYFAFKICSGAANVVGPTMCFEDRMIMSPVKNNVGRGLNIALVNGTTGAVLGQKAFDMYSGDVMHLVKFLKEIPGGALVLVASYDDPGTKMNDESRKLFSDLGSSYAKQLGFRDSWVFIGAKDLRGKSPFEQFLKNSPDTNKYEGWPELLEMEGCMPPKPF",
   "taxon_ids": [
    "9606"
   ],
   "accessions": [
    "A0A0A8K9B4",
    "A0A2J8P7L1",
    "Q96BQ1"
   ],
   "assertions": [
    {
     "assertion_id": "15a1ba52355f6ad7c0ab48869591df3901d66b535d5ba0e9fa231939fb45c49e:GO:0030246",
     "go_id": "GO:0030246",
     "go_name": "carbohydrate binding",
     "go_definition": "\"Binding to a carbohydrate, which includes monosaccharides, oligosaccharides and polysaccharides as well as substances derived from monosaccharides by reduction of the carbonyl group (alditols), by oxidation of one or more hydroxy groups to afford the corresponding aldehydes, ketones, or carboxylic acids, or by replacement of one or more hydroxy group(s) by a hydrogen atom. Cyclitols are generally not regarded as carbohydrates.\" [GOC:mah]",
     "supporting_assertions": [
      {
       "accession": "A0A0A8K9B4",
       "symbol": "FAM3D",
       "relation": "enables",
       "go_id": "GO:0030246",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR001922981"
       ],
       "object_name": "ILEI/PANDER domain-containing protein",
       "taxon": [
        "taxon:9606"
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
       "accession": "A0A2J8P7L1",
       "symbol": "FAM3D",
       "relation": "enables",
       "go_id": "GO:0030246",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR001922981"
       ],
       "object_name": "ILEI/PANDER domain-containing protein",
       "taxon": [
        "taxon:9598"
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
   "sequence_sha256": "592767da9d744e81fbe40d1b8483f97c4d9d92fb3ca160cd559b035e530d03dc",
   "sequence": "MYAAGSRGHSPAFLQPQNGNGHRSPGYVPGKVVPLRPAPPPKNHASAKLTSRSQDAPATFAFSPEEQRTPSESRKRKRHKNTFICFAITSFSFFVALAVILGISSKYAPDENCPDQNPRLRNWDPGQDSAKHIVIKEGDLFRLTSDATVDSIVIQDGGLLVFGDDKDGSKNITLRTRYILIQDGGALHIGAEKCRYRSKATITLYGKSDERESMPIFGKKFIGVEAGGTLELHGAQRTSWTMLARTLHSSGLPFGSYAFEKDFSRGLNVRVIDQDTARVLENEKFDTHEYHNESRRLQEFLRAQEPGRIVAIAVGDSAVKSLLQGTIQMIQDRLGSKLIQGLGYRQAWALVGVIDGGSSSCNESVRNYENHSTGGKALAQGEFYTLDGQKFSVTAYSEWSQGISLSGFRVDIADGVKLHLLDDVSTWEAGDRIVVASTDYSMYQAEELTLLRCPECSRSQVKVKEIPQYLHVGEIIDGIDMRAEVGLLTRNIVIQGEMEDSCYAENHCQFFDYDTFGGHVMIEKNFTSVHLSYVELKHMGQQHMGRYPVHFHLCGDVDSKGGYSQPASVDGLSVHHSFSRCITVHGTSGLLIKDTIGFDTLGHCFFLEDGVEQRNILYHNLGLLTKPGTLLPTDRNSSMCTVMRDGVFGNYVPVPTTDCMAVSTFWIAHPNNHLINNAAAGSQDAGIWYLFHKEPTGESSGLQLLEKPELTPLGIFYNNRVHSNFKAGLFVDKGVKTTNASASDPREYLCLDNSARFRPHQDADPEKPRVAAIIDRLIAFKNNDNGAWVRGGDIIVQNSAFADNGKGLTFASDGSFPSDEGSSQEVTESLFVGESRNYGFQGGQNKYMGTGGIDQKPRTLPRNRTFPIRGFQIYDGPIHLTKSTFKKYVPTPDRYSSAIGFLMKNSWQTTPRNNVSLVKFGPQVSLNVFFGKPGPWFEDCELDGDKNSIFHDIDGSVTGYKDTYVGRMDNYLIRHPNCVNVTKWNAVICSGTYAQVYVQTWNTPNLSMIITRDEYPSHPMVLRGINQRAISPQYQPVVMLEKGYTIHWNGPAPRTTFLYLVNFNKDDWIRVGLCYPANTSFQVTVGFLQRQNGSLSRIEDYEPARSMEELQKKPSERKFYFDSGTGLLFLYLRAHSHRDGHSYCSSQGCERVKIQAATDSKDISNCMAKAYPQYYKKPSAVKRMPAMLTGLCQGCGTHQMVFTSDPHKSYLPVRFQSPGKAEIQRGDPSIISVNGTDFTFRSAGALLLIVDACSVPFRVKEKRMFLSADVSHMEEYFKASIPPRSIVLLSTRGEIKQLNISDSLAVLGLAKPAHLYSKGSVVFLGFSGNFAPSWTKLFTSPDEQGLGVLEQFLPLQMEEYGCSRTGSVHRRDLDLLQQALKVL",
   "taxon_ids": [
    "10090"
   ],
   "accessions": [
    "Q5FWI3"
   ],
   "assertions": [
    {
     "assertion_id": "592767da9d744e81fbe40d1b8483f97c4d9d92fb3ca160cd559b035e530d03dc:GO:0005509",
     "go_id": "GO:0005509",
     "go_name": "calcium ion binding",
     "go_definition": "\"Binding to a calcium ion (Ca2+).\" [GOC:ai]",
     "supporting_assertions": [
      {
       "accession": "Q5FWI3",
       "symbol": "Cemip2",
       "relation": "enables",
       "go_id": "GO:0005509",
       "references": [
        "GO_REF:0000119"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "UniProtKB:Q9UHN6"
       ],
       "object_name": "Cell surface hyaluronidase",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20260702",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000119.html",
        "https://www.uniprot.org/uniprotkb/Q9UHN6/entry"
       ]
      },
      {
       "accession": "Q5FWI3",
       "symbol": "Cemip2",
       "relation": "enables",
       "go_id": "GO:0005509",
       "references": [
        "GO_REF:0000024"
       ],
       "evidence_code": "ISS",
       "with_from": [
        "UniProtKB:Q9UHN6"
       ],
       "object_name": "Cell surface hyaluronidase",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20170428",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000024.html",
        "https://www.uniprot.org/uniprotkb/Q9UHN6/entry"
       ]
      },
      {
       "accession": "Q5FWI3",
       "symbol": "Cemip2",
       "relation": "enables",
       "go_id": "GO:0005509",
       "references": [
        "GO_REF:0000107"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB:Q9UHN6",
        "ensembl:ENSP00000366243"
       ],
       "object_name": "Cell surface hyaluronidase",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20260724",
       "assigned_by": "Ensembl",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000107.html",
        "https://www.uniprot.org/uniprotkb/Q9UHN6/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "16a6fe970745bf4d3c92f8d4ada3438b4dc4ed04ad0c9ad903075f301eabee82",
   "sequence": "MGKRKKSTRKPTKRLVQKLDTKFNCLFCNHEKSVSCTLDKKNSIGTLSCKICGQSFQTRINSLSQPVDVYSDWFDAVEEVNSGRGSDTDDGDEGSDSDYESDSEQDAKTQNDGEIDSDEEEVDSDEERIGQVKRGRGALVDSDDE",
   "taxon_ids": [
    "559292"
   ],
   "accessions": [
    "P36053"
   ],
   "assertions": [
    {
     "assertion_id": "16a6fe970745bf4d3c92f8d4ada3438b4dc4ed04ad0c9ad903075f301eabee82:GO:0008270",
     "go_id": "GO:0008270",
     "go_name": "zinc ion binding",
     "go_definition": "\"Binding to a zinc ion (Zn).\" [GOC:ai]",
     "supporting_assertions": [
      {
       "accession": "P36053",
       "symbol": "ELF1",
       "relation": "enables",
       "go_id": "GO:0008270",
       "references": [
        "PMID:30358795"
       ],
       "evidence_code": "RCA",
       "with_from": [],
       "object_name": "Transcription elongation factor 1",
       "taxon": [
        "taxon:559292"
       ],
       "annotation_date": "20220420",
       "assigned_by": "SGD",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://pubmed.ncbi.nlm.nih.gov/30358795/"
       ]
      }
     ]
    }
   ]
  }
 ]
}
```
