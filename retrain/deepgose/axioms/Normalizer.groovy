// DeepGO-SE's axiom normalizer: GO-Plus (EL profile) in, the normal forms that DeepGO-SE's loader reads out.
// Usage: groovy Normalizer.groovy -i go-plus-el.owl -o go-plus.norm
//
// Based on deepgo2 groovy/Normalizer.groovy (github.com/bio-ontology-research-group/deepgo2, commit 6c77c98; BSD 3-Clause,
// ../LICENSE-deepgo2): jcel 0.24.1 normalization and Manchester rendering, as there. Three changes:
//   1. The two ObjectHasSelf axioms of the current GO-Plus (CL_0000540, UBERON_0001062; outside GO) are removed; any other
//      ObjectHasSelf axiom stops the run.
//   2. jcel gets only the logical axioms (without axiom annotations) and the non-datatype declarations: jcel fails on
//      Declaration(Datatype(xsd:date)) and on IRI-valued axiom annotations, neither of which carries logical content.
//   3. Auxiliary classes keep their own identity. The official script normalizes with a fresh object factory whose
//      auxiliary classes are numbered from jcel's first free identifier, the numbers the translator already gave the
//      ontology's own classes; such auxiliary classes are then rendered as unrelated real classes (visible in the authors'
//      2021 go-plus.norm, e.g. "GO_0022607 and UBERON_0002726 SubClassOf GO_0098775"), and the remaining statements with
//      an auxiliary class are skipped. Here the translator's own factory is used and every normalized statement is
//      written, auxiliary classes as AUX_<6-digit jcel id> (10 characters like GO_0000000, so the loader reads them like
//      any other class) and auxiliary properties as AUXR_<id> or the inverse of the original property. For the GO-Plus
//      of the GO release of 2026-07-26 this writes all 467,171 normalized statements (the official script skips 63,712).
@Grapes([
    @Grab(group="org.semanticweb.elk", module="elk-owlapi", version="0.4.3"),
    @Grab(group="net.sourceforge.owlapi", module="owlapi-api", version="4.2.5"),
    @Grab(group="net.sourceforge.owlapi", module="owlapi-apibinding", version="4.2.5"),
    @Grab(group="net.sourceforge.owlapi", module="owlapi-impl", version="4.2.5"),
    @Grab(group="net.sourceforge.owlapi", module="owlapi-parsers", version="4.2.5"),
    @Grab(group="de.tu-dresden.inf.lat.jcel", module="jcel-owlapi", version="0.24.1"),
    @Grab(group='commons-cli', module='commons-cli', version='1.4'),
    @GrabConfig(systemClassLoader=true)
])

import org.semanticweb.owlapi.apibinding.OWLManager
import org.semanticweb.owlapi.manchestersyntax.renderer.*
import org.semanticweb.owlapi.model.IRI
import org.semanticweb.owlapi.model.OWLDataFactory
import org.semanticweb.owlapi.model.OWLOntology
import org.semanticweb.owlapi.model.OWLOntologyManager
import org.semanticweb.owlapi.model.AxiomType
import org.semanticweb.owlapi.model.OWLObjectHasSelf

import de.tudresden.inf.lat.jcel.owlapi.main.*
import de.tudresden.inf.lat.jcel.reasoner.main.*

import de.tudresden.inf.lat.jcel.ontology.normalization.*
import de.tudresden.inf.lat.jcel.ontology.axiom.extension.*
import de.tudresden.inf.lat.jcel.coreontology.axiom.*
import de.tudresden.inf.lat.jcel.owlapi.translator.*
import de.tudresden.inf.lat.jcel.ontology.axiom.complex.*
import de.tudresden.inf.lat.jcel.coreontology.datatype.IntegerEntityManager
import org.semanticweb.owlapi.model.OWLAxiom
import org.semanticweb.owlapi.model.OWLClass
import org.semanticweb.owlapi.model.OWLObjectProperty
import org.semanticweb.owlapi.model.OWLObjectPropertyExpression

// Renders a normalized statement back to OWL; auxiliary entities are named on first use.
class UniqueAuxiliaryExporter implements NormalizedIntegerAxiomVisitor<OWLAxiom> {
    TranslationRepository repository
    IntegerEntityManager entities
    OWLDataFactory df
    Map<Integer, OWLClass> auxiliaryClasses = [:]
    Map<Integer, OWLObjectProperty> auxiliaryProperties = [:]
    OWLClass cls(Integer id) {
        def known = repository.getOptOWLClass(id)
        if (known.isPresent()) return known.get()
        assert entities.isAuxiliary(id) && id < 1000000: "Unknown class ${id}"
        auxiliaryClasses.computeIfAbsent(id, { df.getOWLClass(IRI.create('http://deepgose.local/AUX_' + String.format('%06d', id))) })
    }
    OWLObjectPropertyExpression prop(Integer id) {
        def known = repository.getOptOWLObjectProperty(id)
        if (known.isPresent()) return known.get()
        if (entities.getAuxiliaryInverseObjectProperties().contains(id)) {
            def original = entities.createOrGetInverseObjectPropertyOf(id)
            assert original != id
            return prop(original).getInverseProperty()
        }
        assert entities.isAuxiliary(id): "Unknown property ${id}"
        auxiliaryProperties.computeIfAbsent(id, { df.getOWLObjectProperty(IRI.create('http://deepgose.local/AUXR_' + id)) })
    }
    OWLAxiom visit(GCI0Axiom a) { df.getOWLSubClassOfAxiom(cls(a.subClass), cls(a.superClass)) }
    OWLAxiom visit(GCI1Axiom a) { df.getOWLSubClassOfAxiom(df.getOWLObjectIntersectionOf([cls(a.leftSubClass),cls(a.rightSubClass)] as Set),cls(a.superClass)) }
    OWLAxiom visit(GCI2Axiom a) { df.getOWLSubClassOfAxiom(cls(a.subClass),df.getOWLObjectSomeValuesFrom(prop(a.propertyInSuperClass),cls(a.classInSuperClass))) }
    OWLAxiom visit(GCI3Axiom a) { df.getOWLSubClassOfAxiom(df.getOWLObjectSomeValuesFrom(prop(a.propertyInSubClass),cls(a.classInSubClass)),cls(a.superClass)) }
    OWLAxiom visit(RI1Axiom a) { df.getOWLReflexiveObjectPropertyAxiom(prop(a.superProperty)) }
    OWLAxiom visit(RI2Axiom a) { df.getOWLSubObjectPropertyOfAxiom(prop(a.subProperty),prop(a.superProperty)) }
    OWLAxiom visit(RI3Axiom a) { df.getOWLSubPropertyChainOfAxiom([prop(a.leftSubProperty),prop(a.rightSubProperty)],prop(a.superProperty)) }
    OWLAxiom visit(RangeAxiom a) { df.getOWLObjectPropertyRangeAxiom(prop(a.property),cls(a.range)) }
    OWLAxiom visit(FunctObjectPropAxiom a) { df.getOWLFunctionalObjectPropertyAxiom(prop(a.property)) }
    OWLAxiom visit(NominalAxiom a) { df.getOWLEquivalentClassesAxiom([cls(a.classExpression),df.getOWLObjectOneOf(repository.getOWLNamedIndividual(a.individual))] as Set) }
}

def cli = new CliBuilder()
cli.with {
usage: 'Self'
  h longOpt:'help', 'this information'
  i longOpt:'input', 'input OWL file', args:1, required:true
  o longOpt:'output', 'output file containing normalized axioms',args:1, required:true
}
def opt = cli.parse(args)
if( !opt ) {
  return
}
if( opt.h ) {
    cli.usage()
    return
}

File ontologyFile = new File(opt.i)
PrintWriter fout = new PrintWriter(new BufferedWriter(new FileWriter(opt.o)))

OWLOntologyManager manager = OWLManager.createOWLOntologyManager()
OWLOntology ontology = manager.loadOntologyFromOntologyDocument(ontologyFile)
println "LOADED axioms=" + ontology.getAxiomCount()

// Change 1: exactly these two ObjectHasSelf axioms are removed.
def expectedSelf = [
    'EquivalentClasses(<http://purl.obolibrary.org/obo/CL_0000540> ObjectHasSelf(<http://purl.obolibrary.org/obo/RO_0020104>) )',
    'EquivalentClasses(<http://purl.obolibrary.org/obo/UBERON_0001062> ObjectHasSelf(<http://purl.obolibrary.org/obo/RO_0020105>) )'
] as Set
def selfAxioms = ontology.getAxioms().findAll { ax -> ax.getNestedClassExpressions().any { it instanceof OWLObjectHasSelf } }
def selfStrings = selfAxioms.collect { it.toString() }.toSet()
if (selfStrings != expectedSelf) {
    throw new IllegalStateException("ObjectHasSelf axioms differ from the expected two: " + selfStrings)
}
manager.removeAxioms(ontology, selfAxioms as Set)
// Change 2: logical axioms without their annotations, and the non-datatype declarations.
def logical = new HashSet(ontology.getLogicalAxioms().collect { it.getAxiomWithoutAnnotations() })
def declarations = ontology.getAxioms(AxiomType.DECLARATION).findAll { !it.getEntity().isOWLDatatype() }.collect { it.getAxiomWithoutAnnotations() }
def input = new HashSet(logical)
input.addAll(declarations)
ontology = manager.createOntology(input)
assert new HashSet(ontology.getLogicalAxioms()) == logical
println "PATCH removed ObjectHasSelf=" + selfAxioms.size() + "; jcel input: logical axioms=" + logical.size() +
        " declarations=" + declarations.size()

// The official normalization, with the translator's factory (change 3).
JcelReasoner ret = new JcelReasoner(ontology, false)
def translator = ret.getTranslator()
Set<ComplexIntegerAxiom> iOntology = ret.getIntegerOntology()
OntologyNormalizer normalizer = new OntologyNormalizer()
IntegerOntologyObjectFactory factory = translator.getOntologyObjectFactory()
Set<NormalizedIntegerAxiom> normalizedOntology = normalizer.normalize(iOntology, factory)
println "NORMALIZED statements=" + normalizedOntology.size()
def renderer = new ManchesterOWLSyntaxOWLObjectRendererImpl()

int written = 0
def exporter = new UniqueAuxiliaryExporter(repository: translator.getTranslationRepository(),
    entities: factory.getEntityManager(), df: manager.getOWLDataFactory())
normalizedOntology.each { ax ->
    fout.println(renderer.render(ax.accept(exporter)))  // no statement may fail here
    written++
}
println "AUXILIARY classes=" + exporter.auxiliaryClasses.size() + " properties=" + exporter.auxiliaryProperties.size()

fout.flush()
fout.close()
println "WRITTEN " + written
