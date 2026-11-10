# Protein function review — batch 1 of 16

25 proteins; 25 assertions. Name your output `AGREEMENT_RESULTS_01.json`.

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
 "batch": 1,
 "source_release": "GOA 2026-07-29",
 "proteins": [
  {
   "sequence_sha256": "7564ae12de848e03bee8c8542ffebcfbfb00e42928723026a4347856fad4a272",
   "sequence": "MATTVSTQRGPVYIGELPQDFLRITPTQQQQQIQLDAQAAQQLQYGGTVGTVGRLSITVVQAKLAKNYGMTRMDPYCRLRLGYAVYETPTAHNGAKNPRWNKVIQCTVPPGVDSFYLEIFDERAFSMDDRIAWTHITIPESLKQGQVEDEWYSLSGRQGDDKEGMINLVMSYTSLPAAMMMPPQPVVLMPTVYQQGVGYVPIAGMPAVCSPGMVPMAMPPPAVAPQPRCNEEDLKAIQDMFPNMDREVIRSVLEAQRGNKDAAINSLLQMGEES",
   "taxon_ids": [
    "10090"
   ],
   "accessions": [
    "A0A6P5Q2D5",
    "A0A8C6G8T4",
    "Q9QZ06"
   ],
   "assertions": [
    {
     "assertion_id": "7564ae12de848e03bee8c8542ffebcfbfb00e42928723026a4347856fad4a272:GO:0032183",
     "go_id": "GO:0032183",
     "go_name": "SUMO binding",
     "go_definition": "\"Binding to the small ubiquitin-like protein SUMO.\" [GOC:mah]",
     "supporting_assertions": [
      {
       "accession": "Q9QZ06",
       "symbol": "Tollip",
       "relation": "enables",
       "go_id": "GO:0032183",
       "references": [
        "GO_REF:0000096"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "RGD:1304676"
       ],
       "object_name": "Toll-interacting protein",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20260723",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000096.html"
       ]
      },
      {
       "accession": "Q9QZ06",
       "symbol": "Tollip",
       "relation": "enables",
       "go_id": "GO:0032183",
       "references": [
        "GO_REF:0000107"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB:A2RUW1",
        "ensembl:ENSRNOP00000049498"
       ],
       "object_name": "Toll-interacting protein",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20260724",
       "assigned_by": "Ensembl",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000107.html",
        "https://www.uniprot.org/uniprotkb/A2RUW1/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "676c344a7056b61c94f638b3577b18df6a3ad959e029a6f31db58a294b50ae46",
   "sequence": "MIVVVAVTKHCTSPAGMKLFIEDVPVLFPYEYIYPEQLEYMTELKKGLDQGGHMVLEMPSGTGKTTTLLSLLIAYLHHHADEKRKVVYCTRTVGEVDKTASELRKILGLWETEEVGSRPLRGVCLTAKKNLCIEPSVVSRRHLDDVDAGCRSLTAPWQTDQRCVYFDTLENQGFDLRPGAYSLSDLKRLGEENHVCPYYLARKALKVADIVIHSFLYIVDPIVAEVTKEYMDENTIVVMDEAHNVDDVCIEAMSLIVTKDDAFQAKDNVKKLNEAVDRMKATNRQKLQEEYDRLVNGLAMAELARSGESTAGVSVQALNIPQDIAEEAIPGSLRQANHFLAFMQRLVDFTHRVVCRISRTYVADPLTFLTKVKEECSVDVRHLRYLTERLKVLLNTLQITSAHSFHNVSLIAQMFMTLSTHYTDDRYEKPGFVVVCEASDPTRPTIPDPVIRLVCVDASLALREVFAKYRSVILTSGTLSPLDIYPKMLGFSPVIAKSFQMTLSRKCIAPVVVTRSSESVNITAEEVTSSFVVRKNPEAQATVSAAYEDLLRGLAKTVPDGIVCFFTGYQYMSEVLLGWHRSGFLKELSRHKLIFVETQSVDDTATALESYRRACDIGRGAIFMSIARGKIAEGIDFDSHYGRAVVMFGVPFLPPIDEPLRQRIHWMEICLGIPGSEYRNFDAMRQASQCIGRVLRNKTDYGMMLLVDRRFALNDKRKKLPIWIQQCLKENTNLSVDAAVAVARGFFKEMAQPWEYKRDLGTTLFSVETLARKGFLKPSKSSVSRKTVDESRATTTVPPPAEEIIRPVREPDYRKRSRE",
   "taxon_ids": [
    "185431"
   ],
   "accessions": [
    "Q57YJ9"
   ],
   "assertions": [
    {
     "assertion_id": "676c344a7056b61c94f638b3577b18df6a3ad959e029a6f31db58a294b50ae46:GO:0016787",
     "go_id": "GO:0016787",
     "go_name": "hydrolase activity",
     "go_definition": "\"Catalysis of the hydrolysis of various bonds, e.g. C-O, C-N, C-C, phosphoric anhydride bonds, etc.\" [ISBN:0198506732]",
     "supporting_assertions": [
      {
       "accession": "Q57YJ9",
       "symbol": "Tb08.11J15.890",
       "relation": "enables",
       "go_id": "GO:0016787",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR000400010"
       ],
       "object_name": "DNA 5'-3' helicase",
       "taxon": [
        "taxon:185431"
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
   "sequence_sha256": "ab00c9f05d6e0f42bb7f78dfe68d68b098b7d030b7ac572e5ac49ad1f3b738a7",
   "sequence": "MKFFIDDLPVLFPYPRIYPEQYSYMCDIKKTLDVGGNCILEMPSGTGKTVSLLSLTVAYQMHYPEHRKIVYCSRTMSEIEKALIELHKLMEFRANELGEVEDFRGLGLTSRKNLCLHPTISKERKGVVVDEKCRRITNGQLKEKIEKGVPTEDSLCSFHEKLYDLEPHNLVPPGVYSFDALIKYCKEEGTCPYFTVRRMIPFCNIIIYSYHYLLDPKIAERVSRELSKDSIVIFDEAHNIDNVCIESLSLDLTDETLKKASRGANKLGEAIEDMKAQDSEKLQNEYEQLVEGLRQAEIARDQELFMANPVLPQDLLDEAIPGNIRKGEHFIAFLKRFIEYLKTRMKVLHVISETPTSFLQHLKELTYIERKPLRFCSERLSLLVRTLEVTEIDDFNALKDIATFATLVSTYDSGFQLILEPYETEGSTVPNPILHFTCLDASIAIKPVFDRFSSVIITSGTISPLDMYPKMLKFETVIQESYAMTLARRSFLPMIVTKGSDQVSISSRFEIRNDPSVVRNYGSLLIEFAKITPDGMVVFFPSYLYMESIISMWQSMGVLDEVWKHKLILVETPDAQETSLALETYRKACSNGRGAVLLSVARGKVSEGIDFDHHYGRTVLMIGIPFQYTESRILKARLEFMRDHFQIKENDFLSFDAMRHAAQCLGRVLRGKDDYGIMVLADRRFARKKNQLPKWIAQALNDSDTNLSTDMALATAKKFLRSLAQPTNPKDSEGVSIWNIEQLQEYQQQFQKDEMDIDEEDLDLL",
   "taxon_ids": [
    "237561"
   ],
   "accessions": [
    "A0A1D8PQH0"
   ],
   "assertions": [
    {
     "assertion_id": "ab00c9f05d6e0f42bb7f78dfe68d68b098b7d030b7ac572e5ac49ad1f3b738a7:GO:0003677",
     "go_id": "GO:0003677",
     "go_name": "DNA binding",
     "go_definition": "\"Any molecular function by which a gene product interacts selectively and non-covalently with DNA (deoxyribonucleic acid).\" [GOC:dph, GOC:jl, GOC:tb, GOC:vw]",
     "supporting_assertions": [
      {
       "accession": "A0A1D8PQH0",
       "symbol": "RAD3",
       "relation": "enables",
       "go_id": "GO:0003677",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001945",
        "InterPro:IPR010614"
       ],
       "object_name": "DNA 5'-3' helicase",
       "taxon": [
        "taxon:237561"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001945/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR010614/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "7692819e223e3c2bca2f5f3fdb0cce8a00265cfffc25f434afdfe9b081742424",
   "sequence": "MKKKRAFLKWAGGKYPLIEEIRRHLPAGERLVEPFVGAGSVFLNTDYDEYLLADINADLVGLYNIVKDETDRFIGDARQMFCPELNVSEEYYRIREAFNLCDDRYQRALYFLYLNRHCYNGLCRYNAKGEFNVPFGRYKKPYFPQEELEFFARKAERAVFVCEHYQQTMLRAQPGAVVYCDPPYAPLSATANFTAYHVNGFGFAEQMRLAELAGELAQRQVPVLISNHETPLTREWYAQADALHVVKVRRTISRKTDGRSKVNELLALYGKR",
   "taxon_ids": [],
   "accessions": [
    "A0A0H3DWW9",
    "A0AAU8PIT3"
   ],
   "assertions": [
    {
     "assertion_id": "7692819e223e3c2bca2f5f3fdb0cce8a00265cfffc25f434afdfe9b081742424:GO:1904047",
     "go_id": "GO:1904047",
     "go_name": "S-adenosyl-L-methionine binding",
     "go_definition": "\"Binding to S-adenosyl-L-methionine.\" [GO_REF:0000067, GOC:BHF, GOC:hal, GOC:TermGenie, PMID:22985361]",
     "supporting_assertions": [
      {
       "accession": "A0A0H3DWW9",
       "symbol": "ETAF_2948",
       "relation": "enables",
       "go_id": "GO:1904047",
       "references": [
        "GO_REF:0000118"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "PANTHER:PTN007825871"
       ],
       "object_name": "Site-specific DNA-methyltransferase (adenine-specific)",
       "taxon": [
        "taxon:718251"
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
       "accession": "A0AAU8PIT3",
       "symbol": "dam",
       "relation": "enables",
       "go_id": "GO:1904047",
       "references": [
        "GO_REF:0000118"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "PANTHER:PTN007825871"
       ],
       "object_name": "Site-specific DNA-methyltransferase (adenine-specific)",
       "taxon": [
        "taxon:1263550"
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
   "sequence_sha256": "05eb8f1a36ba790f88e4aa24996d86f186f71f1a54d3cf3bee948e9304d0ef93",
   "sequence": "MQVTLKTLQQQTFKIDIDPEETVKALKEKIESEKGKDAFPVAGQKLIYAGKILNDDTALKEYKIDEKNFVVVMVTKPKAVSTPAPATTQQSAPASTTAVTSSTTTTVAQAPTPVPALAPTSTPASITPASATASSEPAPASAAKQEKPAEKPAETPVATSPTATDSTSGDSSRSNLFEDATSALVTGQSYENMVTEIMSMGYEREQVIAALRASFNNPDRAVEYLLMGIPGDRESQAVVDPPQAASTGAPQSSAVAAAAATTTATTTTTSSGGHPLEFLRNQPQFQQMRQIIQQNPSLLPALLQQIGRENPQLLQQISQHQEHFIQMLNEPVQEAGGQGGGGGGGSGGIAEAGSGHMNYIQVTPQEKEAIERLKALGFPEGLVIQAYFACEKNENLAANFLLQQNFDED",
   "taxon_ids": [
    "9606"
   ],
   "accessions": [
    "A0A2R9C4U2",
    "G3QJ86",
    "P54727"
   ],
   "assertions": [
    {
     "assertion_id": "05eb8f1a36ba790f88e4aa24996d86f186f71f1a54d3cf3bee948e9304d0ef93:GO:0070628",
     "go_id": "GO:0070628",
     "go_name": "proteasome binding",
     "go_definition": "\"Binding to a proteasome, a large multisubunit protein complex that catalyzes protein degradation.\" [GOC:mah]",
     "supporting_assertions": [
      {
       "accession": "A0A2R9C4U2",
       "symbol": "A0A2R9C4U2",
       "relation": "enables",
       "go_id": "GO:0070628",
       "references": [
        "GO_REF:0000118"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "PANTHER:PTN000065603"
       ],
       "object_name": "UV excision repair protein RAD23",
       "taxon": [
        "taxon:9597"
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
       "accession": "G3QJ86",
       "symbol": "RAD23B",
       "relation": "enables",
       "go_id": "GO:0070628",
       "references": [
        "GO_REF:0000033"
       ],
       "evidence_code": "IBA",
       "with_from": [
        "AGI_LocusCode:AT1G16190",
        "AGI_LocusCode:AT1G79650",
        "AGI_LocusCode:AT3G02540",
        "AGI_LocusCode:AT5G38470",
        "PANTHER:PTN000065572",
        "SGD:S000000763"
       ],
       "object_name": "UV excision repair protein RAD23",
       "taxon": [
        "taxon:9595"
       ],
       "annotation_date": "20180301",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000033.html"
       ]
      },
      {
       "accession": "P54727",
       "symbol": "RAD23B",
       "relation": "enables",
       "go_id": "GO:0070628",
       "references": [
        "GO_REF:0000033"
       ],
       "evidence_code": "IBA",
       "with_from": [
        "AGI_LocusCode:AT1G16190",
        "AGI_LocusCode:AT1G79650",
        "AGI_LocusCode:AT3G02540",
        "AGI_LocusCode:AT5G38470",
        "PANTHER:PTN000065572",
        "SGD:S000000763"
       ],
       "object_name": "Lysine-specific demethylase RAD23B",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20180301",
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
   "sequence_sha256": "b63a601ab2796248420d635d949b1626b4eec60203aaac43be6c4d5d667935f2",
   "sequence": "MSAAEEVDGLGVVRPHYGSVLDNERLTAEEMDERRRQNVAYEYLCHLEEAKRWMEACLGEDLPPTTELEEGLRNGVYLAKLGNFFSPKVVSLKKIYDREQTRYKATGLHFRHTDNVIQWLNAMDEIGLPKIFYPETTDIYDRKNMPRCIYCIHALSLYLFKLGLAPQIQDLYGKVDFTEEEINNMKIELEKYGIQMPAFSKIGGILANELSVDEAALHAAVIAINEAIDRRVAADTFTALKNPNAMLVNLEEGLAPTYQDVLYQAKQDKMTNAKNRTENSDRERDVYEELLTQAEIQGNVNKVNTSSALANISLALEQGCAVTLLKALQSLALGLRGLQTQNSDWYMKQLQSDLQQKRQSGQTDPLQKEEVQAGVDAANSAAQQYQRRLAAVAAINAAIQKGIAEKTVLELMNPEAQLPQVYPFAADLYQKELATLQQQSPEHSLTHPELTVAVEMLSSVALINRALESGDMTTVWKQLSSSVTGLTNIEEENCQRYLDELMKLKAQAHAENNAFITWNDIQACVDHVNLVVHEEHERILAIGLINEALDEGDAQKTLQALQIPAAKLEGVLAEVAQHYQDTLIRAKREKAQETQDESAVLWLDEIQGGIWQSNKDTQEAQRFALGISAINEAVDSGDVGRTLSALRSPDVGLYGVIPECGETYQSDLAEAKKKRLAAGDNNSKWVKHWVKGGYHYYHNLETQAGGWAEPPDFVQNSVQLSREEIQSSISGVTAAYNREQLWLANEGLITKLQACCRGYLVRQEFRSRMNFLKKQIPAITCIQSQWRGYKQKKAYQDRLAYLHSHKDEVVKIQSLARMHQARKRYRDRLQYFRDHINDIIKIQAFIRANKARDDYKTLINAEDPPMIVVRKFVHLLDQSDQDFQEELDLMKMREEVITLIRSNQQLENDLNLMDIKIGLLVKNKITLQDVVSHSKKLTKKNKEQLSDMMMINKQKGGLKALSKEKREKLEAYQHLFYLLQTNPTYLAKLIFQMPQNKSTKFMDSVIFTLYNYASNQREEYLLLRLFQTALQEEIKSKVDQIQEIVTGNPTVIKMVVSFNRGARGQNALRQILAPVVKEIMDDKSLNIKTDPVDIYKSWVNQMESQTGEASKLPYDVTPEQALSHEEVKTRLDNSIRNMRAVTDKFLSAIVSSVDKIPYGMRFIAKVLKDSLHEKFPDAGEDELLKIIGNLLYYRYMNPAIVAPDAFDIIDLSAGGQLTTDQRRNLGSIAKMLQHAASNKMFLGDNAHLSIINEYLSQSYQKFRRFFQVACDVPELQDKFNVDEYSDLVTLTKPVIYISIGEIINTHTLLLDHQDAIAPEHNDPIHELLDDLGEVPTIESLIGESCGNSNDPNKEALAKTEVSLTLTNKFDVPGDENAEMDARTILLNTKRLIVDVIRFQPGETLTEILETPATNEQEAEHQRAMQRRAIRDAKTPDKMKKSKPMKEDNNLSLQEKKEKIQTGLKKLTELGTVDPKNRYQELINDIAKDIRNQRRYRQRRKAELVKLQQTYSALNSKATFYGEQVDYYKSYIKTCLDNLASKGKVSKKPREMKGKKSKKISLKYTAARLHEKGVLLEIEDLQANQFKNVIFEIGPTEEVGDFEVKAKFMGVQMETFMLHYQDLLQLQYEGVAVMKLFDRAKVNVNLLIFLLNKKFYGK",
   "taxon_ids": [
    "10090"
   ],
   "accessions": [
    "Q9JKF1"
   ],
   "assertions": [
    {
     "assertion_id": "b63a601ab2796248420d635d949b1626b4eec60203aaac43be6c4d5d667935f2:GO:0044548",
     "go_id": "GO:0044548",
     "go_name": "S100 protein binding",
     "go_definition": "\"Binding to a S100 protein. S100 is a small calcium and zinc binding protein produced in astrocytes that is implicated in Alzheimer's disease, Down Syndrome and ALS.\" [GOC:jid]",
     "supporting_assertions": [
      {
       "accession": "Q9JKF1",
       "symbol": "Iqgap1",
       "relation": "enables",
       "go_id": "GO:0044548",
       "references": [
        "GO_REF:0000119"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "UniProtKB:P46940"
       ],
       "object_name": "Ras GTPase-activating-like protein IQGAP1",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20260723",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000119.html",
        "https://www.uniprot.org/uniprotkb/P46940/entry"
       ]
      },
      {
       "accession": "Q9JKF1",
       "symbol": "Iqgap1",
       "relation": "enables",
       "go_id": "GO:0044548",
       "references": [
        "GO_REF:0000107"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB:P46940",
        "ensembl:ENSP00000268182"
       ],
       "object_name": "Ras GTPase-activating-like protein IQGAP1",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20260724",
       "assigned_by": "Ensembl",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000107.html",
        "https://www.uniprot.org/uniprotkb/P46940/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "bc09f474c8d396de8c69d6095b0d5a30eacd08de85bbb1a2fbbda4aa6329dba8",
   "sequence": "MAACVPPGEAPRSASGTPTRRQVTIVRIYLDGVYGIGKSTTGRVMASAASGGSPTLYFPEPMAYWRTLFETDVISGIYDTQNRKQQGNLAVDDAALITAHYQSRFTTPYLILHDHTCTLFGGNSLQRGTQPDLTLVFDRHPVASTVCFPAARYLLGDMSMCALMAMVATLPREPQGGNIVVTTLNVEEHIRRLRTRARIGEQIDITLIATLRNVYFMLVNTCHFLRSGRVWRDGWGELPTSCGAYKHRATQMDAFQERVSPELGDTLFALFKTQELLDDRGVILEVHAWALDALMLKLRNLNVFSADLSGTPRQCAAVVESLLPLMSSTLSDFDSASALERAARTFNAEMGV",
   "taxon_ids": [
    "10333"
   ],
   "accessions": [
    "P24425"
   ],
   "assertions": [
    {
     "assertion_id": "bc09f474c8d396de8c69d6095b0d5a30eacd08de85bbb1a2fbbda4aa6329dba8:GO:0005524",
     "go_id": "GO:0005524",
     "go_name": "ATP binding",
     "go_definition": "\"Binding to ATP, adenosine 5'-triphosphate, a universally important coenzyme and enzyme regulator.\" [ISBN:0198506732]",
     "supporting_assertions": [
      {
       "accession": "P24425",
       "symbol": "TK",
       "relation": "enables",
       "go_id": "GO:0005524",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001889"
       ],
       "object_name": "Thymidine kinase",
       "taxon": [
        "taxon:10333"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001889/"
       ]
      },
      {
       "accession": "P24425",
       "symbol": "TK",
       "relation": "enables",
       "go_id": "GO:0005524",
       "references": [
        "GO_REF:0000043"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB-KW:KW-0067"
       ],
       "object_name": "Thymidine kinase",
       "taxon": [
        "taxon:10333"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000043.html"
       ]
      },
      {
       "accession": "P24425",
       "symbol": "TK",
       "relation": "enables",
       "go_id": "GO:0005524",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR000335150"
       ],
       "object_name": "Thymidine kinase",
       "taxon": [
        "taxon:10333"
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
   "sequence_sha256": "c31a46c028ae185f039a3d84753566a08576928895af4882c62345173f49fcb5",
   "sequence": "MIQNVGNHLRRGLASVFSNRTSRKSALRAGNDSAMADGEGYRNPTEVQMSQLVLPCHTNQRGELSVGQLLKWIDTTACLSAERHAGCPCVTASMDDIYFEHTISVGQVVNIKAKVNRAFNSSMEVGIQVASEDLCSEKQWNVCKALATFVARREITKVKLKQITPRTEEEKMEHSVAAERRRMRLVYADTIKDLLANCAIQGDLESRDCSRMVPAEKTRVESVELVLPPHANHQGNTFGGQIMAWMENVATIAASRLCRAHPTLKAIEMFHFRGPSQVGDRLVLKAIVNNAFKHSMEVGVCVEAYRQEAETHRRHINSAFMTFVVLDADDQPQLLPWIRPQPGDGERRYREASARKKIRLDRKYIVSCKQTEVPLSVPWDPSNQVYLSYNNVSSLKMLVAKDNWVLSSEISQVRLYTLEDDKFLSFHMEMVVHVDAAQAFLLLSDLRQRPEWDKHYRSVELVQQVDEDDAIYHVTSPALGGHTKPQDFVILASRRKPCDNGDPYVIALRSVTLPTHRETPEYRRGETLCSGFCLWREGDQLTKCCWVRVSLTELVSASGFYSWGLESRSKGRRSDGWNGKLAGGHLSTLKAIPVAKINSRFGYLQDT",
   "taxon_ids": [
    "9606"
   ],
   "accessions": [
    "Q8WXI4"
   ],
   "assertions": [
    {
     "assertion_id": "c31a46c028ae185f039a3d84753566a08576928895af4882c62345173f49fcb5:GO:0008289",
     "go_id": "GO:0008289",
     "go_name": "lipid binding",
     "go_definition": "\"Binding to a lipid.\" [GOC:ai]",
     "supporting_assertions": [
      {
       "accession": "Q8WXI4",
       "symbol": "ACOT11",
       "relation": "enables",
       "go_id": "GO:0008289",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR002913"
       ],
       "object_name": "Acyl-coenzyme A thioesterase 11",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR002913/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "6dd347e086bacaaab71164878a19ce11a66fd90f561371990b6d8b1a00ba0b39",
   "sequence": "MAANYSSTSSRKEHVKVTSEPQPGFLERLSETSGGMFVGLMTFLLSFYLIFTNEGRALKTATSLAEGLSLVVSPDSIHSVAPENEGRLVHIIGALRTSKLLSDPNYGVHLPAVKLRRHVEMYQWVETEESSEYTEDGQVKKETKYSYNTEWRSEIVNSRNFDREIGHKNPSAMAVESFTATAPFVQIGRFFLSAGLIDKIDNFKALSLAKLEDPHVDIIRRGDFFYHSENPKYPEVGDVRVSFSYAGLSSDDPDLGPAHVVTVIARQRGDQLIPYSTKSGDTLLLLHHGDFSAEEVFRREQKSNSMKTWGLRAAGWMAMFMGLNLMTRILYTLVDWFPVFRDLVNIGLKAFAFCVATSLTLLTVAAGWLFYRPLWAALIGCLALVPIIIARTRVPAKKLE",
   "taxon_ids": [
    "10090"
   ],
   "accessions": [
    "Q9DBS1"
   ],
   "assertions": [
    {
     "assertion_id": "6dd347e086bacaaab71164878a19ce11a66fd90f561371990b6d8b1a00ba0b39:GO:0005249",
     "go_id": "GO:0005249",
     "go_name": "voltage-gated potassium channel activity",
     "go_definition": "\"Enables the transmembrane transfer of a potassium ion by a voltage-gated channel. A voltage-gated channel is a channel whose open state is dependent on the voltage across the membrane in which it is embedded.\" [GOC:mtg_transport, ISBN:0815340729]",
     "supporting_assertions": [
      {
       "accession": "Q9DBS1",
       "symbol": "Tmem43",
       "relation": "enables",
       "go_id": "GO:0005249",
       "references": [
        "DOI:10.1101/2022.11.08.515259"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "UniProtKB:Q9BTV4"
       ],
       "object_name": "Transmembrane protein 43",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20250425",
       "assigned_by": "MGI",
       "annotation_extension": "part_of(GO:0006813)",
       "gene_product_form_id": "",
       "links": [
        "https://www.uniprot.org/uniprotkb/Q9BTV4/entry"
       ]
      },
      {
       "accession": "Q9DBS1",
       "symbol": "Tmem43",
       "relation": "enables",
       "go_id": "GO:0005249",
       "references": [
        "GO_REF:0000107"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB:Q9BTV4",
        "ensembl:ENSP00000303992"
       ],
       "object_name": "Transmembrane protein 43",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20260724",
       "assigned_by": "Ensembl",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000107.html",
        "https://www.uniprot.org/uniprotkb/Q9BTV4/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "c1874e4a2a96216109a3cabf02a9d5841209f7eb273f766971fc49bee6d42452",
   "sequence": "MKVEEILEKALELVIPDEEEVRKGREAEEELRRRLDELGVEYVFVGSYARNTWLKGSLEIDVFLLFPEEFSKEELRERGLEIGKAVLDSYEIRYAEHPYVHGVVKGVEVDVVPCYKLKEPKNIKSAVDRTPFHHKWLEGRIKGKENEVRLLKGFLKANGIYGAEYKVRGFSGYLCELLIVFYGSFLETVKNARRWTRRTVIDVAKGEVRKGEEFFVVDPVDEKRNVAANLSLDNLARFVHLCREFMEAPSLGFFKPKHPLEIEPERLRKIVEERGTAVFAVKFRKPDIVDDNLYPQLERASRKIFEFLERENFMPLRSAFKASEEFCYLLFECQIKEISRVFRRMGPQFEDERNVKKFLSRNRAFRPFIENGRWWAFEMRKFTTPEEGVRSYASTHWHTLGKNVGESIREYFEIISGEKLFKEPVTAELCEMMGVKD",
   "taxon_ids": [
    "224325"
   ],
   "accessions": [
    "O28126"
   ],
   "assertions": [
    {
     "assertion_id": "c1874e4a2a96216109a3cabf02a9d5841209f7eb273f766971fc49bee6d42452:GO:0000166",
     "go_id": "GO:0000166",
     "go_name": "nucleotide binding",
     "go_definition": "\"Binding to a nucleotide, any compound consisting of a nucleoside that is esterified with (ortho)phosphate or an oligophosphate at any hydroxyl group on the ribose or deoxyribose.\" [GOC:mah, ISBN:0198547684]",
     "supporting_assertions": [
      {
       "accession": "O28126",
       "symbol": "cca",
       "relation": "enables",
       "go_id": "GO:0000166",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR000101564"
       ],
       "object_name": "CCA-adding enzyme",
       "taxon": [
        "taxon:224325"
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
   "sequence_sha256": "3313530069f164536fd6b99bc8fe7812d1de792ac1167f930f191f4afed775ec",
   "sequence": "MYVNFFNKLNIPKCLRTNFNFKRRFSETKQNIQIKYIKENKKKYVFYTAVTSLSLSIALTYAYEKFVLYKWNRRYDYYYNPNIKFFEKKEKKEGKKDKKNTTKHIILVRHGQYERRYKDDENSKRLTKEGCKQADITGKKLKDILNNKKVSVIYHSDMIRAKETANIISKYFPDANLINDPNLNEGTPYLPDPLPRHSKFDAQKIKEDNKRINKAYETYFYKPSGDEDEYQLVICHGNVIRYFLCRALQIPLFAWLRFSSYNCGITWLVLDDEGSVVLREFGSVSHLPFESVTYF",
   "taxon_ids": [
    "36329"
   ],
   "accessions": [
    "Q8I1V2",
    "W7K9K2"
   ],
   "assertions": [
    {
     "assertion_id": "3313530069f164536fd6b99bc8fe7812d1de792ac1167f930f191f4afed775ec:GO:0004722",
     "go_id": "GO:0004722",
     "go_name": "protein serine/threonine phosphatase activity",
     "go_definition": "\"Catalysis of the reaction: protein serine phosphate + H2O = protein serine + phosphate, and protein threonine phosphate + H2O = protein threonine + phosphate.\" [GOC:bf]",
     "supporting_assertions": [
      {
       "accession": "Q8I1V2",
       "symbol": "PF3D7_0413500",
       "relation": "enables",
       "go_id": "GO:0004722",
       "references": [
        "GO_REF:0000033"
       ],
       "evidence_code": "IBA",
       "with_from": [
        "FB:FBgn0023517",
        "PANTHER:PTN008573713",
        "UniProtKB:Q96HS1"
       ],
       "object_name": "Serine/threonine-protein phosphatase PGAM5, mitochondrial",
       "taxon": [
        "taxon:36329"
       ],
       "annotation_date": "20240528",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000033.html",
        "https://www.uniprot.org/uniprotkb/Q96HS1/entry"
       ]
      },
      {
       "accession": "W7K9K2",
       "symbol": "PFNF54_01017",
       "relation": "enables",
       "go_id": "GO:0004722",
       "references": [
        "GO_REF:0000118"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "PANTHER:PTN001075958"
       ],
       "object_name": "Serine/threonine-protein phosphatase PGAM5, mitochondrial",
       "taxon": [
        "taxon:5843"
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
   "sequence_sha256": "a3bc566baab8685fe7d30aa1ed957a74debd409e26cbaa2da29ea3e3f2ab327e",
   "sequence": "MFADLDYDIEEDKLGIPTVPGKVTLQKDAQNLIGISIGGGAQYCPCLYIVQVFDNTPAALDGTVAAGDEITGVNGRSIKGKTKVEVAKMIQEVKGEVTIHYNKLQADPKQGMSLDIVLKKVKHRLVENMSSGTADALGLSRAILCNDGLVKRLEELERTAELYKGMTEHTKNLLRAFYELSQTHRAFGDVFSVIGVREPQPAASEAFVKFADAHRSIEKFGIRLLKTIKPMLTDLNTYLNKAIPDTRLTIKKYLDVKFEYLSYCLKVKEMDDEEYSCIALGEPLYRVSTGNYEYRLILRCRQEARARFSQMRKDVLEKMELLDQKHVQDIVFQLQRLVSTMSKYYNDCYAVLRDADVFPIEVDLAHTTLAYGLNQEEFTDGEEEEEEEDTAAGEPSRDTRGAAGPLDKGGSWCDS",
   "taxon_ids": [
    "9606"
   ],
   "accessions": [
    "A0A0D9R3D4",
    "A0A2K5WCA7",
    "A0A5F8AC16",
    "F6Z9S9",
    "H9G232",
    "Q9NRD5",
    "U3FMB7"
   ],
   "assertions": [
    {
     "assertion_id": "a3bc566baab8685fe7d30aa1ed957a74debd409e26cbaa2da29ea3e3f2ab327e:GO:0071933",
     "go_id": "GO:0071933",
     "go_name": "Arp2/3 complex binding",
     "go_definition": "\"Binding to an Arp2/3 complex, a protein complex that contains two actin-related proteins, Arp2 and Arp3, and five novel proteins (ARPC1-5).\" [GOC:mah]",
     "supporting_assertions": [
      {
       "accession": "Q9NRD5",
       "symbol": "PICK1",
       "relation": "enables",
       "go_id": "GO:0071933",
       "references": [
        "GO_REF:0000024"
       ],
       "evidence_code": "ISS",
       "with_from": [
        "UniProtKB:Q9EP80"
       ],
       "object_name": "PRKCA-binding protein",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20131129",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000024.html",
        "https://www.uniprot.org/uniprotkb/Q9EP80/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "6d5da28fb2e605b14e209af020bdc5711b85c03bed7c9caeb05cfffec20cf5b2",
   "sequence": "MLRVFVIFSLFIAAINASGEFTVLNRPKAISFKGNDALESHYVGDVLYASMGNAVSGDTNWNGLTINDPFNLAKGVILVHVQGIGHVTTAGNVKTYELTGSGTDASLNALAAELEAANEPVCDINFEQFDDGVQAWKSCFGDFEAPAAKPTKHLNPSLHTADKQFLQEVGFINSAADHLAEMAKPSNVLMLRVSVDGVAKAHGEKSVAVEEANKLLSAAISRLLAASQKSSDSVLFVQTTEKDVAASRAKRDTIAASTTNPYNLAVYYGSDYPVIFNIILWFMVVFGLSLLAICYAIAAMDPGRDSIIYRMTSTRIKKDN",
   "taxon_ids": [
    "7227"
   ],
   "accessions": [
    "Q9VHG4"
   ],
   "assertions": [
    {
     "assertion_id": "6d5da28fb2e605b14e209af020bdc5711b85c03bed7c9caeb05cfffec20cf5b2:GO:0038023",
     "go_id": "GO:0038023",
     "go_name": "signaling receptor activity",
     "go_definition": "\"Receiving a signal and transmitting it in the cell to initiate a change in cell activity. A signal is a physical entity or change in state that is used to transfer information in order to trigger a response.\" [GOC:bf, GOC:signaling]",
     "supporting_assertions": [
      {
       "accession": "Q9VHG4",
       "symbol": "ATP6AP2",
       "relation": "enables",
       "go_id": "GO:0038023",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR012493"
       ],
       "object_name": "ATPase H(+)-transporting accessory protein 2",
       "taxon": [
        "taxon:7227"
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
   "sequence_sha256": "9df783062a35d2ac9b279ebdee01d886b0282cb47a352c01a4c41c770548d8c7",
   "sequence": "MPNDPSDMNNELHNENTIRILISSDPHVGYGEKDPVRGNDSFVSFNEILEIARERDVDMILLGGDIFHDNKPSRKALYQALRSLRLNCLGDKPCELELLSDTSLTTGDTAVCNINYLDPNINVAIPVFSIHGNHDDPSGDGRYSALDILQVTGLVNYFGRVPENDNIVVSPILLQKGFTKLALYGISNVRDERLYHSFRENKVKFLRPDLYRDEWFNLLTVHQNHSAHTPTSYLPESFIQDFYDFVLWGHEHECLIDGSYNPTQKFTVVQPGSTIATSLSPGETAPKHCGILNITGKDFHLEKIRLRTVRPFIMKDIILSEVSSIPPMVENKKEVLTYLISKVEEAITEANAQWYEAQGTVPVVENEKPPLPLIRLRVDYTGGYQTENPQRFSNRFVGRVANATDVVQFYLKKKYTRSKRNDGLYTSAVEDIKINSLRVESLVNEYLKTNRLECLPEDSLGEAVVNFVEKDDRDAIKECVETQLNKQINLLVKKRVTEENLEQEISSIINDLPKISTTKRKDYEELPEEVSETSINIAEHTPVLKHTSSLLDHHSPLATSSSEHEMEATPSPALLKKTNKRRELPSSLTKKNTRTPQRSKEVKKVPARKLSQSTKKSDKNTQSTLLFYDPSSTTEAQYLDNEDDEILDD",
   "taxon_ids": [
    "284812"
   ],
   "accessions": [
    "Q09683"
   ],
   "assertions": [
    {
     "assertion_id": "9df783062a35d2ac9b279ebdee01d886b0282cb47a352c01a4c41c770548d8c7:GO:0004520",
     "go_id": "GO:0004520",
     "go_name": "DNA endonuclease activity",
     "go_definition": "\"Catalysis of the cleavage of ester linkages within deoxyribonucleic acid by creating internal breaks.\" [GOC:mah, ISBN:0198547684]",
     "supporting_assertions": [
      {
       "accession": "Q09683",
       "symbol": "rad32",
       "relation": "enables",
       "go_id": "GO:0004520",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR003701"
       ],
       "object_name": "Double-strand break repair protein rad32",
       "taxon": [
        "taxon:284812"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR003701/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "d0a06a66686e8769d886c37d3efa56005a237ab07fb6dbffc8bc5592bf4e148d",
   "sequence": "MKVLALRHSVAQVYADTQVYVHDDTKDSYENAFLISNLTTHNILYLNYSIKTLEILNKSGIAAIALQSLEELFTLIRCNFTYDYELDIIYLHDYSYYTNNEIRTDQHWITKTNIEEYLLPGWKLTYVGYNGSETRGHYNFSFKCQNAATDDDLIIEYIYSEALDFQNFMLKKIKERMTTSLPIARLSNRVFRDKLFPSLLKEHKNVVNVGPRNESMFTFLNYPTIKQFSNGAYLVKDTIKLKQERWLGKRISQFDIGQYKNMLNVLTAIYYYYNLYKSKPIIYMIGSAPSYWIYDVRHYSDFFFETWDPLDTPYSSIHHKELFFINDVKKLKDNSILYIDIRTDRGNADWKKWRKTVEEQTINNLDIAYEYLRTGKAKVCCVKMTAMDLELPISAKLLHHPTTEIRSEFYLLLDTWDLTNIRRFIPKGVLYSFINNIITENVFIQQPFKVKVLNDSYIVALYALSNDFNNRSEVIKLINNQKQSLITVRINNTFKDEPKVGFKNIYDWTFLPTDFDTKEAIITSYDGCLGLFGLSISLASKPTGNNHLFILSGTDKYYKLDQFANHTSISRRSHQIRFSESATSYSGYIFRDLSNNNFNLIGTNIENSVSGHVYNALIYYRYNYSFDLKRWIYLHSIDKVDIEGGKYYELAPIELIYACRSAKEFATLQDDLTVLRYSNEIENYINTVYSITYADDPNYFIGIQFRNIPYKYDVKIPHLTFGVLHISDNMVPDVIDILKIMKNELFKMDITTSYTYMLSDGIYVANVSGVLSTYFKIYNVFYKNQITFGQSRMFIPHITLSFNNMRTVRIETTKLQIKSIYLRKIKGDTVFDMVE",
   "taxon_ids": [
    "450149"
   ],
   "accessions": [
    "A2T3S5"
   ],
   "assertions": [
    {
     "assertion_id": "d0a06a66686e8769d886c37d3efa56005a237ab07fb6dbffc8bc5592bf4e148d:GO:0004482",
     "go_id": "GO:0004482",
     "go_name": "mRNA 5'-cap (guanine-N7-)-methyltransferase activity",
     "go_definition": "\"Catalysis of the reaction: S-adenosyl-L-methionine + G(5')pppR-RNA = S-adenosyl-L-homocysteine + m7G(5')pppR-RNA. m7G(5')pppR-RNA is mRNA containing an N7-methylguanine cap; R may be guanosine or adenosine.\" [EC:2.1.1.56]",
     "supporting_assertions": [
      {
       "accession": "A2T3S5",
       "symbol": "A2T3S5",
       "relation": "enables",
       "go_id": "GO:0004482",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR011181"
       ],
       "object_name": "Protein VP3",
       "taxon": [
        "taxon:450149"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR011181/"
       ]
      },
      {
       "accession": "A2T3S5",
       "symbol": "A2T3S5",
       "relation": "enables",
       "go_id": "GO:0004482",
       "references": [
        "GO_REF:0000003"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "EC:2.1.1.56"
       ],
       "object_name": "Protein VP3",
       "taxon": [
        "taxon:450149"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000003.html"
       ]
      },
      {
       "accession": "A2T3S5",
       "symbol": "A2T3S5",
       "relation": "enables",
       "go_id": "GO:0004482",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR000866602"
       ],
       "object_name": "Protein VP3",
       "taxon": [
        "taxon:450149"
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
   "sequence_sha256": "369958a465f10dfe6cc02cdd26ceed8fc8917c0c5a8f4e4fc4a2d322479aecb9",
   "sequence": "MPVEYKRRPTFWEIFKARCSTADLGPISLNWFEELSSEAPPYNSEPPEESEYKPHGYEPQLFKTPQRNPPYHQFASTPIMFKERSQTLPLDQSPFRELGKVVASSKHKTHSKKKTKVDPVVDVASPPLKSCLSESPLTLRCTQAVLQREKPVVSGSLFYTPKLKEGQTPKPISESLGVEVDPDMSWTSSLATPPTLSSTVLIARDEEARSSVTPADSPATLKSCFSNHNESPQKNDRSVPSVIDSENKNQQEAFSQGLGKMLGDSSGKRNSFKDCLRKPIPNILEDGETAVDTSEEDSFSLCFPKRRTRNLQKMRMGKTRKKIFSETRTDELSEEARRQTDDKNSFVFEMELRESDPLDPGVTSQKPFYSQNEEICNEAVQCSDSRWSQSNLSGLNETQTGKITLPHISSHSQNISEDFIDMKKEGTGSITSEKSLPHISSLPEPEKMFSEETVVDKEHEGQHFESLEDSIAGKQMVSRTSQAACLSPSIRKSIFKMREPLDETLGTVFSDSMTNSTFTEEHEASACGLGILTACSQREDSICPSSVDTGSWPTTLTDTSATVKNAGLISTLKNKKRKFIYSVSDDASLQGKKLQTHRQLELTNLSAQLEASAFEVPLTFTNVNSGIPDSSDKKRCLPNDPEEPSLTNSFGTATSKEISYIHALISQDLNDKEAIVIEEKPQPYTAREADFLLCLPERTCENDQKSPKVSNGKEKVLVSACLPSAVQLSSISFESQENPLGDHNGTSTLKLTPSSKLPLSKADMVSREKMCKMPEKLQCESCKVNIELSKNILEVNEICILSENSKTPGLLPPGENIIEVASSMKSQFNQNAKIVIQKDQKGSPFISEVAVNMNSEELFPDSGNNFAFQVTNKCNKPDLGSSVELQEEDLSHTQGPSLKNSPMAVDEDVDDAHAAQVLITKDSDSLAVVHDYTEKSRNNIEQHQKGTEDKDFKSNSSLNMKSDGNSDCSDKWSEFLDPVLNHNFGGSFRTASNKEIKLSEHNVKKSKMFFKDIEEQYPTRLACIDIVNTLPLANQKKLSEPHIFDLKSVTTVSTQSHNQSSVSHEDTDTAPQMLSSKQDFHSNNLTTSQKAEITELSTILEESGSQFEFTQFRKPSHIAQNTSEVPGNQMVVLSTASKEWKDTDLHLPVDPSVGQTDHSKQFEGSAGVKQSFPHLLEDTCNKNTSCFLPNINEMEFGGFCSALGTKLSVSNEALRKAMKLFSDIENSEEPSAKVGPRGFSSSAHHDSVASVFKIKKQNTEKSFDEKSSKCQVTLQNNIEMTTCIFVGRNPEKYIKNTKHEDSYTSSQRNNLENSDGSMSSTSGPVYIHKGDSDLPADQGSKCPESCTQYAREENTQIKENISDLTCLEIMKAEETCMKSSDKKQLPSDKMEQNIKEFNISFQTASGKNTRVSKESLNKSVNIFNRETDELTVISDSLNSKILHGINKDKMHTSCHKKAISIKKVFEDHFPIVTVSQLPAQQHPEYEIESTKEPTLLSFHTASGKKVKIMQESLDKVKNLFDETQYVRKTASFSQGSKPLKDSKKELTLAYEKIEVTASKCEEMQNFVSKETEMLPQQNYHMYRQTENLKTSNGTSSKVQENIENNVEKNPRICCICQSSYPVTEDSALAYYTEDSRKTCVRESSLSKGRKWLREQGDKLGTRNTIKIECVKEHTEDFAGNASYEHSLVIIRTEIDTNHVSENQVSTLLSDPNVCHSYLSQSSFCHCDDMHNDSGYFLKNKIDSDVPPDMKNAEGNTISPRVSATKERNLHPQTINEYCVQKLETNTSPHANKDVAIDPSLLDSRNCKVGSLVFITAHSQETERTKEIVTDNCYKIVEQNRQSKPDTCQTSCHKVLDDSKDFICPSSSGDVCINSRKDSFCPHNEQILQHNQSMSGLKKAATPPVGLETWDTSKSIREPPQAAHPSRTYGIFSTASGKAIQVSDASLEKARQVFSEMDGDAKQLSSMVSLEGNEKPHHSVKRENSVVHSTQGVLSLPKPLPGNVNSSVFSGFSTAGGKLVTVSESALHKVKGMLEEFDLIRTEHTLQHSPIPEDVSKILPQPCAEIRTPEYPVNSKLQKTYNDKSSLPSNYKESGSSGNTQSIEVSLQLSQMERNQDTQLVLGTKVSHSKANLLGKEQTLPQNIKVKTDEMKTFSDVPVKTNVGEYYSKESENYFETEAVESAKAFMEDDELTDSEQTHAKCSLFTCPQNETLFNSRTRKRGGVTVDAVGQPPIKRSLLNEFDRIIESKGKSLTPSKSTPDGTVKDRSLFTHHMSLEPVTCGPFCSSKERQGAQRPHLTSPAQELLSKGHPWRHSALEKSPSSPIVSILPAHDVSATRTERTRHSGKSTKVFVPPFKMKSQFHGDEHFNSKNVNLEGKNQKSTDGDREDGNDSHVRQFNKDLMSSLQSARDLQDMRIKNKERRHLRLQPGSLYLTKSSTLPRISLQAAVGDRAPSACSPKQLYIYGVSKECINVNSKNAEYFQFDIQDHFGKEDLCAGKGFQLADGGWLIPSNDGKAGKEEFYRALCDTPGVDPKLISSIWVANHYRWIVWKLAAMEFAFPKEFANRCLNPERVLLQLKYRYDVEIDNSRRSALKKILERDDTAAKTLVLCISDIISPSTKVSETSGGKTSGEDANKVDTIELTDGWYAVRAQLDPPLMALVKSGKLTVGQKIITQGAELVGSPDACAPLEAPDSLRLKISANSTRPARWHSRLGFFRDPRPFPLPLSSLFSDGGNVGCVDIIVQRVYPLQWVEKTVSGLYIFRSEREEEKEALRFAEAQQKKLEALFTKVHTEFKDHEEDTTQRCVLSRTLTRQQVHALQDGAELYAAVQYASDPDHLEACFSEEQLRALNNYRQMLNDKKQARIQSEFRKALESAEKEEGLSRDVTTVWKLRVTSYKKKEKSALLSIWRPSSDLSSLLTEGKRYRIYHLAVSKSKSKFERPSIQLTATKRTQYQQLPVSSETLLQVYQPRESLHFSRLSDPAFQPPCSEVDVVGVVVSVVKPIGLAPLVYLSDECLNLLVVKFGIDLNEDIKPRVLIAASNLQCQPESTSGVPTLFAGHFSIFSASPKEAYFQEKVNNLKHAIENIDTFYKEAEKKLIHVLEGDSPKWSTPNKDPTREPHAASTCCASDLLGSGGQFLRISPTGQQSYQSPLSHCTLKGKSMPLAHSAQMAAKSWSGENEIDDPKTCRKRRALDFLSRLPLPSPVSPICTFVSPAAQKAFQPPRSCGTKYATPIKKEPSSPRRRTPFQKTSGVSLPDCDSVADEELALLSTQALTPDSVGGNEQAFPGDSTRNPQPAQRPDQQVGPRSRKESLRDCRGDSSEKLAVES",
   "taxon_ids": [
    "10090"
   ],
   "accessions": [
    "P97929"
   ],
   "assertions": [
    {
     "assertion_id": "369958a465f10dfe6cc02cdd26ceed8fc8917c0c5a8f4e4fc4a2d322479aecb9:GO:0043015",
     "go_id": "GO:0043015",
     "go_name": "gamma-tubulin binding",
     "go_definition": "\"Binding to the microtubule constituent protein gamma-tubulin.\" [GOC:jl]",
     "supporting_assertions": [
      {
       "accession": "P97929",
       "symbol": "Brca2",
       "relation": "enables",
       "go_id": "GO:0043015",
       "references": [
        "GO_REF:0000119"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "UniProtKB:P51587"
       ],
       "object_name": "Breast cancer type 2 susceptibility protein homolog",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20260723",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000119.html",
        "https://www.uniprot.org/uniprotkb/P51587/entry"
       ]
      },
      {
       "accession": "P97929",
       "symbol": "Brca2",
       "relation": "enables",
       "go_id": "GO:0043015",
       "references": [
        "GO_REF:0000107"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB:P51587",
        "ensembl:ENSP00000369497"
       ],
       "object_name": "Breast cancer type 2 susceptibility protein homolog",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20260724",
       "assigned_by": "Ensembl",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000107.html",
        "https://www.uniprot.org/uniprotkb/P51587/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "a425fb2811dcee2c86b5b9a61b04d14e44daa96c2b12186ae8e497b6467dc5bd",
   "sequence": "MTIAIGRAPAERGWFDILDDWLKRDRFVFVGWSGILLFPCAYLALGGWLTGTTFVTSWYTHGLASSYLEGCNFLTVAVSTPANSMGHSLLLLWGPEAQGDFTRWCQLGGLWTFIALHGAFGLIGFMLRQFEIARLVGVRPYNAIAFSAPIAVFVSVFLIYPLGQSSWFFAPSFGVAAIFRFLLFFQGFHNWTLNPFHMMGVAGVLGGALLCAIHGATVENTLFQDGEGASTFRAFNPTQAEETYSMVTANRFWSQIFGIAFSNKRWLHFFMLFVPVTGLWMSAIGVVGLALNLRSYDFISQEIRAAEDPEFETFYTKNLLLNEGIRAWMAPQDQPHENFVFPEEVLPRGNAL",
   "taxon_ids": [
    "197221"
   ],
   "accessions": [
    "Q8CM25"
   ],
   "assertions": [
    {
     "assertion_id": "a425fb2811dcee2c86b5b9a61b04d14e44daa96c2b12186ae8e497b6467dc5bd:GO:0005506",
     "go_id": "GO:0005506",
     "go_name": "iron ion binding",
     "go_definition": "\"Binding to an iron (Fe) ion.\" [GOC:ai]",
     "supporting_assertions": [
      {
       "accession": "Q8CM25",
       "symbol": "psbD1",
       "relation": "enables",
       "go_id": "GO:0005506",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR000076254"
       ],
       "object_name": "Photosystem II D2 protein",
       "taxon": [
        "taxon:197221"
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
   "sequence_sha256": "7daf75aa27336830845f26d94819026d80b7c71a78d18ab20ae274bf10b1a7bf",
   "sequence": "MRALRAPVLVMGLVLLICTAAQSDASANKVEKTLRRIMEIMRQVENSADDESAQKTESAPEPKDTHHLKTASGETILEIFPRDLSRKEKFITILTGRRLSNPLRHQSCVFETSDLVVFHPSGPLYFGPKCKKDVYRLYHNTRDCTIPAHYKRCARLLTRLAGTRKCQEG",
   "taxon_ids": [
    "7955"
   ],
   "accessions": [
    "E7F5F0"
   ],
   "assertions": [
    {
     "assertion_id": "7daf75aa27336830845f26d94819026d80b7c71a78d18ab20ae274bf10b1a7bf:GO:0030971",
     "go_id": "GO:0030971",
     "go_name": "receptor tyrosine kinase binding",
     "go_definition": "\"Binding to a receptor that possesses protein tyrosine kinase activity.\" [GOC:mah]",
     "supporting_assertions": [
      {
       "accession": "E7F5F0",
       "symbol": "alkal2a",
       "relation": "enables",
       "go_id": "GO:0030971",
       "references": [
        "GO_REF:0000033"
       ],
       "evidence_code": "IBA",
       "with_from": [
        "PANTHER:PTN002010604",
        "UniProtKB:Q6UX46",
        "UniProtKB:Q6UXT8"
       ],
       "object_name": "ALK and LTK ligand 2a",
       "taxon": [
        "taxon:7955"
       ],
       "annotation_date": "20190822",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000033.html",
        "https://www.uniprot.org/uniprotkb/Q6UX46/entry",
        "https://www.uniprot.org/uniprotkb/Q6UXT8/entry"
       ]
      },
      {
       "accession": "E7F5F0",
       "symbol": "alkal2a",
       "relation": "enables",
       "go_id": "GO:0030971",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR029364"
       ],
       "object_name": "ALK and LTK ligand 2a",
       "taxon": [
        "taxon:7955"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR029364/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "c604b3a6fe9b42c9b2de2f55cf48c99bed5fd1244e0028f1687b02d1276283f6",
   "sequence": "MEDSEDDVVVVSCDTSMKEKVKAKLVEIRKFVPFIRRVRIDFQDTLSKVQGHRLDALVNLLDREDVSMSSLNKIEVIIDKLRTRFNPRIEIDTGEIIDITENTEYDSTASSSSGSSAQPPAQRAKSSDKGGLLNGSSLAGRLNASLNPASDATKDHLSTPKQLELERNSKVSCSTSKIKTSSISAKISVFPATSMDLNIPKTTSAKASDEGQRSPAEPRAALQAIVQDTKTPTIPEPTSPAALKHSSLRGSRGFLAVMQKALIEEKKQRASEQKTDKETNGVTQLETNFSRRSYTTSSQSTCRSSEISVRAENPDFKRRSTSLVQHAPLQEASPGQSKKDLPISLSVQGLPALVSASTASPANTLEEARKKLAALKYGLGTTVPSMPPLASNINDPRGRKGINLPETNNNKDNDLGIALQSPPPMRTPSPIPPPPRMKAGTWASFSNVPQESAFTGQHAVQRNSVPPGDSRAFGDALAHEPRSFYGTDSREPRDPRIWKSKTSQQQQHQQAPQAQIPPYSSDPRRSISTYSGFEESANRGQTNSKNNNNHNQNQNVTANGSAYIANWSGNPYGNGPNPRPPFPSNQNGNEYAGGFNGSHSFRGGFRGGHNKRGFGRHNDVPRTYGEHRKAKARAEAEAKAKAEAEAKAKAAAEAKAKAAAEVRQLETEVSREMEAQEKNKQQKEKPEESEAEKSTIAVTQVPELDTSYRNVNLGVLNKKLDFRIPKKTLPPATTITSTSPVNGNGENPSCPSNSPTSKSCDANQDKDTYKNKDRYLNKAKAKDKVDKGNEVSENNLDKSEKLEKSQDKKANDKENKSDKKEKKRLNREPEKKSKVENPLEIVDSNSVVSEESSENTDNVENEPPLSETNASPVPELATSTQDSQQDQSVSEELDILAKNRRMSGTRIKTPISSTGNPALKRRADDDVEDKLENPTKKNCAKWEAKPDKEKSEDDTIDKIKSMKVTKFADVEMKVTEESQSAEEEEITEQKESTEEEEGTEHKKSTEEKDKPPKISKIKIVLTPIAHTTQVVRPNDGFKNNQEKILDNMATDEHDDEEVPGPPAQFLRRIMQRRNSLAPTYMKPMVDKDKIASSSFTYEDLPEQKRGNQNARNLAIIFEKTSDNCSVSTQNIINGKRRTRGCETSFNETQLSRNIFGMGQINRSRPKATEGKATRAKISVPAKDSTRADDDPIPTPNLKRKRQAIHKETEDDVEMKPKKARLEAQEINGVSVTPDEQQVENNVEVTQKEVEAISSEPLLSSEVEPTRKPRTKPRKNELDKLNDDIAQMYYGEEVMRATSRRACTRRSRTSSHTRTSSQHSRTSSVSRTDSISTVSDISSIIVRNTARRGRGIRSSENGINRATFNASLNAKKPKLCRVRIKRCAALMEMIKDQEKEEQEKKEQKNKEPKKKKVGVQKKPLKSKPKRENSVILNTNPEWHSISKAVIKCVVCSKWVRRSPLSHYMMCHKEHYAARLPPDVLKELRAGRGNRPDYWVSQRGGYTLHFTCPFCQKPLLLCQKGMIEHLIGHMGESRFYCSNCNMPQNRLSRLLDHTASCGPGAKPLSSKTVCLPMSVHVCHICQFMQYSKENMDRHLTVQHGLTKEELESVEREELMLCDTTDVPYADSNKDGSAVGPEIQKNKPKSKASAALSNTTKKNKQQKKTNNSRFMKMVKKSVSLLKREREEQDDQNMTQANEGSEVPEIPPPPPEIEPLFVVNECLMTSEMDTDMEEVLEQPVQHMSLMVDEKPVTLLSGATEQLEPSVPDPEPVVPSAQDDGKDVNEDEDVDVEAVVDSLQSHTDQTATSMLAEVSLAELAGDVLDGIGSDASDYEMDDNSEQVDTTNKNGYGDDDDDALTDDWVDLETAKRNSKSAKSIFRVFNRFCSRLNKLPRSSRAVPSNGSENSDGSDNNDDDGDNPDPSELMPTMQPLEPEPEMGDSSTSTGAKSLSERVENVGFQKPSSDEDQNRVAASYYCVQPGCTFLFSNELEGLENHFALEHPLVRWSGKCGMCRQKITATETNLRISEELRHMRDVHMKDISTLPPPQSSAVESPAVIESCLNQREPVPESEPDPVPELPKLRVRRFTGDRLVVDSQAEKSQPVAIVVSDDDNPRNGMLRDLLAADPRPPNQQLDLQAAGLGEFLCAKPDSPSTEPVKQTPVIVGYSSGLGLKIGQVLSRTQISANSRLSPVVNDPLPEKSSAPAAVEENRNRFRCMATNCNFVAHKLMFMREHMKFHSYSFSSTGHLNCAYCSHVAVDVDDYLRHGVIIHDLAPRSELESSTGPPSVTQKIRDMLSQRENGRVPPPTPQVTLSDVVLGLLECTGYSEDKLYACPQKGCIVRLTDEQLVNHLRYHIRSTHQGSELVKCKFCTKAMHPPALRTHLQQYHARHSIFCGICLATSVNQRIMMYHMSTVHSKAYGRPNARLAFVSLPVKIDASKKNVESEFYVAVVEQPFGNLQMQDFQRKLFDEMDRRRSGTKTYFRSSEVHILPTQPTFQRPLYCTECPFSTTSRVNMQMHLYEHKDETIREASKLADLIVPATSSVLTVSASTLVAPPRPGKDSEKPSTSGQSGDAATEQLNPDVPGTHKPIKPPLRYVPPDQRYRCGFLRCSVLCFSESAVRKHMQANHKYSEVVRCPHCKNCQGQFGVDKYFDHLAMHKRHIFQCGACSRHNSRRVIERHIQERHNIQDVDMIVHRHNDSNKTTEARWLKAPKLARHSLMEYTCNLCLKYFPTTVQIMAHAASVHKRNYQYHCPYCEFGGNLATALIEHILREHPEREVQPVQIYQRIVCKNKQTLGFYCTTCHEVASSFQKIAMHCDKEHKSRNPVQCPHCIFGHLAERQVVLHIQEKHPHERGLAMVQFERVLNDIPNSISWEIGRPIEVEPEKEIPNNGESAFLPLSQRQVVTEVVDLLDSDDEADEYGEQDDAKIVEFACTHCDGTNTNLPDLRSQHWAREHPDQPFYFRVQPMLLCSECKRFRGNAKALREHLRATHSIRSIVAADIRRPMECAYCDYRYKNRHDLAKHISEIGHLPNDLKHVTDDEIDALMLLSASGSGGAVNEYYQCGLCSVVMPTKETIVQHGQVEHCKPDERFCFRQLVSPVIYHCSFCMFNSTDELTTLRHMVDHYSRFLVCHFCTRSQPGGFDEYIQHCYTYHRDDIKSFRDVHTFSDLKRYLSQVHYQFQNGLIITKSSLRYTRYKSDKCMLELDAELMAKAQRPPIPRLHIRLKSTGVQMQSPEGADVEKPVSLLRITKRRKTLNPGELLRSFREENEVQPQPPASSTSSGTAPSPAAGSVFNLFKRRNSLVVRPATSNLDQH",
   "taxon_ids": [
    "7227"
   ],
   "accessions": [
    "Q9W3W6"
   ],
   "assertions": [
    {
     "assertion_id": "c604b3a6fe9b42c9b2de2f55cf48c99bed5fd1244e0028f1687b02d1276283f6:GO:0003700",
     "go_id": "GO:0003700",
     "go_name": "DNA-binding transcription factor activity",
     "go_definition": "\"A transcription regulator activity that modulates transcription of gene sets via selective and non-covalent binding to a specific double-stranded genomic DNA sequence (sometimes referred to as a motif) within a cis-regulatory region. Regulatory regions include promoters (proximal and distal) and enhancers. Genes are transcriptional units, and include bacterial operons.\" [GOC:txnOH-2018]",
     "supporting_assertions": [
      {
       "accession": "Q9W3W6",
       "symbol": "sov",
       "relation": "enables",
       "go_id": "GO:0003700",
       "references": [
        "GO_REF:0000033"
       ],
       "evidence_code": "IBA",
       "with_from": [
        "MGI:MGI:104897",
        "PANTHER:PTN008617116",
        "RGD:621069",
        "UniProtKB:Q13127",
        "UniProtKB:Q9NQV6"
       ],
       "object_name": "C2H2-type domain-containing protein",
       "taxon": [
        "taxon:7227"
       ],
       "annotation_date": "20250620",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000033.html",
        "https://www.uniprot.org/uniprotkb/Q13127/entry",
        "https://www.uniprot.org/uniprotkb/Q9NQV6/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "356cddd146593ce5445e5a944cbd179758583451f4c69aff265d55bed936101e",
   "sequence": "MSEDSEEEDYSDRSISDDDDLDEDSFMKFVSDDIHPCTLLAADSIGDPFFPRTTQILLEYQLGRWVPRLRGPRDLYGVSSSGPLSPTRWPYHCEVIDEKVQHIEWTPFVPEPVYVPTGLEIEPVYPNSKEDTVVYLAEDDHLCKAYKEPCFVYSRVGGNRTSLKQPVDNCDNTLVFEARFESGNLQKVVKVADHEYELTVRPDLFTNKHTQWYYFQVTNTQAEIVYRFTIVNFTKPASLYNRGMKPLFYSEKEAKTHNIGWQRIGDQIKYYKNNLGQDGRHFFSLTWTFQFPHSQDTCYFAHCYPYTYSNLQEYLSGINSDPVRSKFCKIRVLCHTLARNMVYVLTITTPLKTSDSKRKAVILTARVHPGETNSSWIMKGFLDYILGDSSDARLLRDTFIFKVVPMLNPDGVIVGNYRCSLAGRDLNRNYTSLLKESFPSVWYTRNMINRLMEKREVILYCDLHGHSRKQNIFMYGCDGSSRSKTKGLYLQQRIFPLMLSKNCPNIFSFSACKFNVQKSKEGTGRVVMWKMGIRNSFTLEATFCGSTLGNKRGTHFGTKDLESMGYHFCDSLLDYCDPDRSKYYQCLKELEEMEKHLSSERVSDNTDTSLVEISLDVESSSRGSDSSESNDTQTYLLKVTSQARNKKKYLKTKRERNAILANCQNNMQEVYGKEHLLQRHDESNSDGNDPRIDAPDVYVAHCFRRPLPNQGVVKIPGQRFYPGKTWSSSQRMIKSLNKDHRTCILETCKNPIQEVQSRGINIHESCFKMAKCPMNKRPSHWIEKTRIPTESHHQLKSKAKRCSSFQSKRTGTNWTDDEKRIYRDKRIAQTQEILKYLLPIVESSQNRKSTQMNNLINPIANLQQHQLIPTACINRRRYSIPWTPTRNLPFKAQRNLMTDTSEWLQSVPLGSFESLLPLCNLQKKTKHFELWGKKAKDVQLATSQWEAVPLSSNMDASIIRGNSVLQPKEFTMRSSKQRIPYLTKTSKKPSESDGLLTFQLKIHRNS",
   "taxon_ids": [
    "10090"
   ],
   "accessions": [
    "Q8CDP0"
   ],
   "assertions": [
    {
     "assertion_id": "356cddd146593ce5445e5a944cbd179758583451f4c69aff265d55bed936101e:GO:0015631",
     "go_id": "GO:0015631",
     "go_name": "tubulin binding",
     "go_definition": "\"Binding to monomeric or multimeric forms of tubulin, including microtubules.\" [GOC:clt]",
     "supporting_assertions": [
      {
       "accession": "Q8CDP0",
       "symbol": "Agbl3",
       "relation": "enables",
       "go_id": "GO:0015631",
       "references": [
        "GO_REF:0000033"
       ],
       "evidence_code": "IBA",
       "with_from": [
        "MGI:MGI:1918244",
        "MGI:MGI:2159437",
        "MGI:MGI:2441745",
        "MGI:MGI:3646469",
        "PANTHER:PTN007589887"
       ],
       "object_name": "Cytosolic carboxypeptidase 3",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20240613",
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
   "sequence_sha256": "a97eafaf3a98fbcb89da68d5e4ff5ed954d04c39eb45be376a4733a6c2602f45",
   "sequence": "MSNTMQGAHYPQMAFDPTWQMHPAYATATAFNPYNVSFDFTRSPNTSGGSDESMSDGSKIDPKRSPKYLEKRMKNNEAAKKSRASRKHREQKNQTENELLKRKNAALEEELKQAKCELAQMQITIRDMSIEREAYRRENEMLKMVNNKFADSKFEPPQPLRDMTNCNLPYKYELLS",
   "taxon_ids": [],
   "accessions": [
    "P46505"
   ],
   "assertions": [
    {
     "assertion_id": "a97eafaf3a98fbcb89da68d5e4ff5ed954d04c39eb45be376a4733a6c2602f45:GO:0003700",
     "go_id": "GO:0003700",
     "go_name": "DNA-binding transcription factor activity",
     "go_definition": "\"A transcription regulator activity that modulates transcription of gene sets via selective and non-covalent binding to a specific double-stranded genomic DNA sequence (sometimes referred to as a motif) within a cis-regulatory region. Regulatory regions include promoters (proximal and distal) and enhancers. Genes are transcriptional units, and include bacterial operons.\" [GOC:txnOH-2018]",
     "supporting_assertions": [
      {
       "accession": "P46505",
       "symbol": "zip-8",
       "relation": "enables",
       "go_id": "GO:0003700",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR004827",
        "InterPro:IPR031106",
        "InterPro:IPR046347"
       ],
       "object_name": "bZIP transcription factor 8",
       "taxon": [
        "taxon:6239"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR004827/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR031106/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR046347/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "5571568bcd9d69d97def748aa3f69859223b4638c83761f7328fb76a3940917c",
   "sequence": "MAMRELVEGECGGANPLMKLATHFTQDKALRQEGLRPGPWPPGASAAETVSKPLGVGTEDELVSEFLQDQNATLVSRAPQTFKMDDLLAEMQEIEQSNFRQAPQRAPGVADLALSENWAQEFLAAGDAVDVAQDYNETDWSQEFIAEVTDPLSVSPARWAEEYLEQSEEKLWLGDQEGSSTADRWYDEYHPEEDLQHTASDFVSKVDDPKLANSEFLKFVRQIGEGQVSLESAAGSGGAQAEQWAAEFIQQQGTSEAWVDQFTRPGNKIAALQVEFERAKSAIESDVDFWDKLQAELEEMAKRDAEAHPWLSDYDDLTSASYDKGYQFEEENPLRDHPQPFEEGLHRLEEGDLPNAVLLFEAAVQQDPKHMEAWQYLGTTQAENEQELLAISALRRCLELKPDNRTALMALAVSFTNESLQRQACETLRDWLRYSPAYAHLVAPGEEGATGAGPSKRILGSLLSDSLFLEVKDLFLAAVRLDPTSIDPDVQCGLGVLFNLSGEYDKAVDCFTAALSVRPNDYLMWNKLGATLANGNQSEEAVAAYRRALELQPGYIRSRYNLGISCINLGAHREAVEHFLEALNMQRKSRGPRGEGGAMSENIWSTLRLALSMLGQSDAYGAADARDLSALLAMFGLPQ",
   "taxon_ids": [
    "10090"
   ],
   "accessions": [
    "O09012"
   ],
   "assertions": [
    {
     "assertion_id": "5571568bcd9d69d97def748aa3f69859223b4638c83761f7328fb76a3940917c:GO:0005515",
     "go_id": "GO:0005515",
     "go_name": "protein binding",
     "go_definition": "\"Binding to a protein.\" [GOC:curators]",
     "supporting_assertions": [
      {
       "accession": "O09012",
       "symbol": "Pex5",
       "relation": "enables",
       "go_id": "GO:0005515",
       "references": [
        "PMID:9668159"
       ],
       "evidence_code": "ISS",
       "with_from": [
        "UniProtKB:P50542-1"
       ],
       "object_name": "Peroxisomal targeting signal 1 receptor",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20070911",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "UniProtKB:O09012-1",
       "links": [
        "https://pubmed.ncbi.nlm.nih.gov/9668159/",
        "https://www.uniprot.org/uniprotkb/P50542-1/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "147d7dc64f2bd9410a6d14de72d572e3c4860bfbd7b452ee0b18940283c3be46",
   "sequence": "MPEGHTLHRLARLHQRRFAGAPVSVSSPQGRFADSASALNGRVLRRASAWGKHLFHHYVGGPVVHVHLGLYGTFTEWARPTDGWLPEPAGQVRMRMVGAEFGTDLRGPTVCESIDDGEVADVVARLGPDPLRSDANPSSAWSRITKSRRPIGALLMDQTVIAGVGNVYRNELLFRHRIDPQRPGRGIGEPEFDAAWNDLVSLMKVGLRRGKIIVVRPEHDHGLPSYLPDRPRTYVYRRAGEPCRVCGGVIRTALLEGRNVFWCPVCQT",
   "taxon_ids": [
    "83332"
   ],
   "accessions": [
    "A0ABX2VII1",
    "A0ACD3TSJ7",
    "A5U5G0",
    "P64159",
    "P9WNB8",
    "P9WNB9"
   ],
   "assertions": [
    {
     "assertion_id": "147d7dc64f2bd9410a6d14de72d572e3c4860bfbd7b452ee0b18940283c3be46:GO:0140078",
     "go_id": "GO:0140078",
     "go_name": "class I DNA-(apurinic or apyrimidinic site) endonuclease activity",
     "go_definition": "\"Catalysis of the cleavage of an AP site 3' of the baseless site by a beta-lyase mechanism, leaving an unsaturated aldehyde, termed a 3'-(4-hydroxy-5-phospho-2-pentenal) residue, and a 5'-phosphate.\" [PMID:1698278, RHEA:66592]",
     "supporting_assertions": [
      {
       "accession": "A0ABX2VII1",
       "symbol": "A7J32_10470",
       "relation": "enables",
       "go_id": "GO:0140078",
       "references": [
        "GO_REF:0000003"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "EC:4.2.99.18"
       ],
       "object_name": "Endonuclease 8 1",
       "taxon": [
        "taxon:1844474"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000003.html"
       ]
      },
      {
       "accession": "P9WNB9",
       "symbol": "nei1",
       "relation": "enables",
       "go_id": "GO:0140078",
       "references": [
        "GO_REF:0000003"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "EC:4.2.99.18"
       ],
       "object_name": "Endonuclease 8 1",
       "taxon": [
        "taxon:83332"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000003.html"
       ]
      },
      {
       "accession": "P9WNB9",
       "symbol": "nei1",
       "relation": "enables",
       "go_id": "GO:0140078",
       "references": [
        "GO_REF:0000116"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "RHEA:66592"
       ],
       "object_name": "Endonuclease 8 1",
       "taxon": [
        "taxon:83332"
       ],
       "annotation_date": "20260727",
       "assigned_by": "RHEA",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000116.html"
       ]
      },
      {
       "accession": "A0ACD3TSJ7",
       "symbol": "FHI80_09610",
       "relation": "enables",
       "go_id": "GO:0140078",
       "references": [
        "GO_REF:0000003"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "EC:4.2.99.18"
       ],
       "object_name": "Endonuclease 8 1",
       "taxon": [
        "taxon:1305738"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000003.html"
       ]
      },
      {
       "accession": "A5U5G0",
       "symbol": "MRA_2490",
       "relation": "enables",
       "go_id": "GO:0140078",
       "references": [
        "GO_REF:0000003"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "EC:4.2.99.18"
       ],
       "object_name": "Endonuclease 8 1",
       "taxon": [
        "taxon:419947"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000003.html"
       ]
      },
      {
       "accession": "P9WNB8",
       "symbol": "nei1",
       "relation": "enables",
       "go_id": "GO:0140078",
       "references": [
        "GO_REF:0000003"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "EC:4.2.99.18"
       ],
       "object_name": "Endonuclease 8 1",
       "taxon": [
        "taxon:83331"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000003.html"
       ]
      },
      {
       "accession": "P9WNB8",
       "symbol": "nei1",
       "relation": "enables",
       "go_id": "GO:0140078",
       "references": [
        "GO_REF:0000116"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "RHEA:66592"
       ],
       "object_name": "Endonuclease 8 1",
       "taxon": [
        "taxon:83331"
       ],
       "annotation_date": "20260727",
       "assigned_by": "RHEA",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000116.html"
       ]
      },
      {
       "accession": "P64159",
       "symbol": "nei1",
       "relation": "enables",
       "go_id": "GO:0140078",
       "references": [
        "GO_REF:0000003"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "EC:4.2.99.18"
       ],
       "object_name": "Endonuclease 8 1",
       "taxon": [
        "taxon:233413"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000003.html"
       ]
      },
      {
       "accession": "P64159",
       "symbol": "nei1",
       "relation": "enables",
       "go_id": "GO:0140078",
       "references": [
        "GO_REF:0000116"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "RHEA:66592"
       ],
       "object_name": "Endonuclease 8 1",
       "taxon": [
        "taxon:233413"
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
   "sequence_sha256": "8a07dda38060381f735a7faba509a27f650d654190a1f97d6ea197364db5495e",
   "sequence": "MRKLLLLFLFICGLVLHAQAQKNVELSSPNGEIKVSVRLTDKIYYDVACNDEVLLKDNYLQLKLKDKILGEHPKLSGQKQNKIKETLTPVVPLKFSTINNEYNQLLLNFKGNYSVEFRAFNDGIAYRFITRQKGEIEVLGEDLSLQFPTDYLLHLQQPGGFKTAYEEPYTHISSHSWKPSDKMSVLPVLIDTRKSYKILISESDLSDYPCMFLKGNGTNGISSIFPKVPLAFGEDGDRSLKIEKEAEYIAKTSGKRNFPWRYFVITKEDKQLVENTMTYKLADKNQLEDVSWIKPGQVSWEWWNGATPYGQDVTFKAGCNLETYKYFIDFASKFGIPYIIMDEGWAKSTRDPYTPNPDVDLHELIRYGKEKNVGIVLWLTWLTVEKNFDLFKTFNEWGIKGVKIDFMDRSDQWMVNYYERVAQEAARHHLFVDFHGSFKPAGLEYKYPNVLSYEGVRGMEQMGGCKPENSIYLPFMRNAVGPMDYTPGAMISMQPNIYRSERPNSASIGTRAYQMALFVIFESGLQMLADNPTLYYRNEDCTRFITQVPVTWDETVVLEAKVGEYVIVAKRKGEKWFIGGMTNDKENEREFEINLDFLKDGRTYRVTSFEDGINAGYQAMDYRMKSAAMNKNQKLSVKMARNGGFSAVLE",
   "taxon_ids": [],
   "accessions": [
    "Q8A2L1"
   ],
   "assertions": [
    {
     "assertion_id": "8a07dda38060381f735a7faba509a27f650d654190a1f97d6ea197364db5495e:GO:0030246",
     "go_id": "GO:0030246",
     "go_name": "carbohydrate binding",
     "go_definition": "\"Binding to a carbohydrate, which includes monosaccharides, oligosaccharides and polysaccharides as well as substances derived from monosaccharides by reduction of the carbonyl group (alditols), by oxidation of one or more hydroxy groups to afford the corresponding aldehydes, ketones, or carboxylic acids, or by replacement of one or more hydroxy group(s) by a hydrogen atom. Cyclitols are generally not regarded as carbohydrates.\" [GOC:mah]",
     "supporting_assertions": [
      {
       "accession": "Q8A2L1",
       "symbol": "BT_3294",
       "relation": "enables",
       "go_id": "GO:0030246",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR014718"
       ],
       "object_name": "Alpha-glucosidase",
       "taxon": [
        "taxon:226186"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR014718/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "888bc916e4d799e8d414bb20aa5d55fe1b1b6d8869151f16e6efc7e2d5c4953f",
   "sequence": "MDVTRLLLATLLVFLCFFTANSHLPPEEKLRDDRSLRSNSSVNLLDVPSVSIVALNKKSKQIGRKAAEKKRSSKKEASMKKVVRPRTPLSAPCVATRNSCKPPAPACCDPCASCQCRFFRSACSCRVLSLNC",
   "taxon_ids": [
    "9606"
   ],
   "accessions": [
    "P42127"
   ],
   "assertions": [
    {
     "assertion_id": "888bc916e4d799e8d414bb20aa5d55fe1b1b6d8869151f16e6efc7e2d5c4953f:GO:0070996",
     "go_id": "GO:0070996",
     "go_name": "type 1 melanocortin receptor binding",
     "go_definition": "\"Binding to a type 1 melanocortin receptor.\" [GOC:BHF, GOC:mah]",
     "supporting_assertions": [
      {
       "accession": "P42127",
       "symbol": "ASIP",
       "relation": "enables",
       "go_id": "GO:0070996",
       "references": [
        "GO_REF:0000024"
       ],
       "evidence_code": "ISS",
       "with_from": [
        "UniProtKB:Q03288"
       ],
       "object_name": "Agouti-signaling protein",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20250714",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000024.html",
        "https://www.uniprot.org/uniprotkb/Q03288/entry"
       ]
      },
      {
       "accession": "P42127",
       "symbol": "ASIP",
       "relation": "enables",
       "go_id": "GO:0070996",
       "references": [
        "GO_REF:0000107"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB:Q03288",
        "ensembl:ENSMUSP00000029123"
       ],
       "object_name": "Agouti-signaling protein",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20260724",
       "assigned_by": "Ensembl",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000107.html",
        "https://www.uniprot.org/uniprotkb/Q03288/entry"
       ]
      }
     ]
    }
   ]
  }
 ]
}
```
