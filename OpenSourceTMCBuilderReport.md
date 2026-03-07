# **Computational Architectures for Automated 3D Generation of Transition Metal Complexes: An Exhaustive Analysis of Permissive Open-Source Python Libraries**

## **Executive Summary**

The automated generation of three-dimensional (3D) structures for transition metal complexes represents one of the most significant challenges in modern computational chemistry and materials science. Unlike organic molecules, which follow predictable covalent bonding patterns and hybridization rules (e.g., ![][image1], ![][image2], ![][image3]), transition metal complexes exhibit a vast diversity of coordination numbers, oxidation states, and geometric isomers (e.g., octahedral, square planar, tetrahedral, trigonal bipyramidal). Furthermore, the bonding in these systems often involves dative (coordinate-covalent) interactions that defy standard valency checks inherent in traditional cheminformatics tools designed for organic chemistry.

As the demand for high-throughput virtual screening (HTVS) of catalysts, metal-organic frameworks (MOFs), and metallo-supramolecular cages grows, so does the need for robust software capable of constructing these complex architectures algorithmically. For industrial applications and integrated research workflows, the software license is a critical determinant of viability. "Permissive" licenses—specifically MIT, BSD, and Apache—allow for the unrestricted integration of these libraries into commercial pipelines, proprietary databases, and larger software suites without the "copyleft" obligations associated with the GNU General Public License (GPL).

This report provides a comprehensive, expert-level analysis of the open-source Python ecosystem for transition metal complex generation, strictly limiting the scope to libraries released under permissive licenses. The analysis identifies **RDKit** (BSD) as the foundational engine for ligand handling, but highlights its limitations in metal coordination. It establishes **stk** (Supramolecular Toolkit, MIT) as the premier topological assembler for complex supramolecular architectures, **cgbind** (MIT) as the specialized standard for template-based metallocage screening, and **Molassembler** (BSD) as the emerging leader in rigorous graph-theoretical stereochemical construction. **mBuild** (MIT) is evaluated for its utility in hierarchical soft-matter interfaces. The report details the algorithmic underpinnings, architectural logic, and integration strategies for these tools, offering a definitive guide for computational chemists and software architects seeking to build automated inorganic discovery platforms.

## ---

**1\. The Computational Landscape of Transition Metal Stereochemistry**

To evaluate the efficacy of any software tool in this domain, one must first delineate the specific computational hurdles imposed by transition metal chemistry. These challenges render standard organic cheminformatics algorithms insufficient and necessitate the specialized architectures found in the libraries analyzed in this report.

### **1.1. The Dative Bond and Valence Violations**

The fundamental data structure in cheminformatics is the molecular graph, where atoms are nodes and bonds are edges. In organic chemistry, the valence of an atom is strictly defined by its connectivity and bond orders. A neutral carbon atom with four single bonds is stable; a neutral nitrogen atom with four single bonds is hypervalent and chemically invalid in most contexts.

* **The Conflict:** In a transition metal complex, a neutral ligand (e.g., pyridine or amine) donates a lone pair to a metal center. In a standard connectivity table, this appears as a single bond. If the nitrogen is neutral, RDKit and similar tools interpret this as a valence violation (4 bonds to a neutral N), triggering "Sanitization Errors" that halt processing.1  
* **Software Responses:** Libraries must either bypass these checks (disabling sanitization), use explicit charge separation models (e.g., ![][image4]), or implement specialized bond types. Recent updates to RDKit (2020.09+) have introduced BondType.DATIVE to formally represent these interactions without affecting the valence count of the donor atom, a critical development utilized by downstream libraries like **stk**.2

### **1.2. Coordination Geometry and Isomerism**

The geometry of a transition metal center is not solely determined by steric repulsion (VSEPR theory) as in organic chemistry, but by electronic factors including d-orbital occupancy and Ligand Field Stabilization Energy (LFSE).

* **Geometric Diversity:** A four-coordinate metal center can be tetrahedral (e.g., ![][image5]) or square planar (e.g., ![][image6]). Six-coordinate complexes are typically octahedral but can exhibit significant distortion (Jahn-Teller effect).  
* **Stereochemical Complexity:** An octahedral complex with the formula ![][image7] can exist in multiple stereoisomeric forms (e.g., *all-cis*, *all-trans*, *cis-trans-cis*). Distinguishing and selectively generating these isomers is a non-trivial graph isomorphism problem. "Distance Geometry" (DG) algorithms used for organic conformers often fail to find specific isomers without explicit constraints or templates.4

### **1.3. The "Bite Angle" Constraint in Polydentate Ligands**

A defining feature of stable metal complexes is the chelate effect, often achieved by polydentate ligands (e.g., bipyridine, porphyrins, ethylenediamine). The geometric compatibility between the metal's preferred bond angles and the ligand's structural constraints is governed by the "bite angle."

* **Algorithmic Challenge:** A generative algorithm must ensure that a bidentate ligand is placed such that its donor atoms map to adjacent (cis) sites on the metal. It must also verify that the ligand backbone is flexible enough to span the required distance without inducing unphysical strain. Failure to account for bite angles leads to fragmented or highly distorted structures during optimization.6

## ---

**2\. Foundational Infrastructure: RDKit (BSD License)**

**RDKit** stands as the ubiquitous foundation of open-source cheminformatics. Released under the **BSD-3-Clause license**, it provides the essential data structures for molecular representation. While primarily optimized for organic chemistry, its permissive license and extensive Python API make it the base layer upon which almost all specialized metal-complex builders are constructed.1

### **2.1. Native Capabilities and Limitations**

RDKit provides robust algorithms for parsing SMILES strings, canonicalizing molecular graphs, and generating 3D conformers via Distance Geometry (DG) and ETKDG (Experimental-Torsion Knowledge Distance Geometry).1

* **Inorganic SMILES Parsing:** RDKit can parse inorganic SMILES (e.g., \[Fe+2\]), but it treats the metal as just another atom. It does not inherently "know" that Iron(II) prefers octahedral coordination.  
* **Conformer Generation:** When asked to embed a metal complex from a connectivity table, RDKit's ETKDG algorithm—trained on organic crystal structures—often produces geometrically erratic results. It may force a square-planar platinum complex into a tetrahedral geometry because it lacks the specific torsion rules for square-planar metal centers.4

### **2.2. The "Dummy Atom" Protocol**

To overcome RDKit's native limitations in metal coordination, advanced workflows (such as those in **stk** and **cgbind**) employ a "Dummy Atom" protocol using RDKit as a geometry engine.

* **Mechanism:**  
  1. **Ligand Conformers:** The organic ligand is built and optimized in RDKit as a standalone molecule.  
  2. **Vector Definition:** The coordination vector is defined relative to the donor atom (e.g., the lone pair direction on a pyridine nitrogen).  
  3. **Dummy Placement:** A dummy atom (atomic number 0 or a placeholder like Du) is placed at a specific distance along this vector, representing the ideal position of the metal center.12  
  4. **Alignment:** The system aligns the ligand's dummy atom with a target vertex on a predefined metal template (e.g., the vertices of a perfect octahedron) using the Kabsch algorithm.  
  5. **Assembly:** Once aligned, the dummy atoms are deleted, and the actual bond to the metal is created in the connectivity table.13

### **2.3. Handling Dative Bonds in RDKit**

Recent versions of RDKit (2020.09 and later) have explicitly addressed the issue of coordinate bonds.

* **The \-\> Notation:** SMILES strings can now represent dative bonds using the \-\> symbol (e.g., N-\>\[Fe\]). This allows RDKit to construct the graph without assigning formal covalent valence to the bond, effectively solving the "hypervalent nitrogen" sanitization error that plagued earlier workflows.2  
* **Impact on Downstream Tools:** This update has forced updates in dependent libraries like **stk**, which now utilize this feature to ensure valid graph representations for metal complexes, improving the stability of descriptors and property calculations derived from these models.2

## ---

**3\. The Topological Assembly Engine: stk (MIT License)**

**stk** (Supramolecular Toolkit) represents the state-of-the-art in permissive, algorithmic molecular assembly. Developed by the Jelfs Materials Group (Imperial College London) and Turcani et al., it abstracts the chemical assembly process into a graph-topological problem. Licensed under the **MIT License**, it is designed for extreme modularity and extensibility, making it the primary recommendation for complex, multi-component transition metal systems.2

### **3.1. Architectural Philosophy: The Topology Graph**

The core innovation of **stk** is the separation of atomic details from structural topology.

* **Topology Graph Class:** A molecule is defined by a TopologyGraph, which consists of Vertex objects (nodes) and Edge objects (connections). For a metal complex, the topology graph defines the ideal geometry (e.g., Octahedral, SquarePlanar) independent of the specific atoms involved.16  
* **Vertex Responsibilities:** Each vertex in the graph has a specific coordinate in 3D space. Crucially, the vertex object dictates how a building block placed upon it should be aligned. An Octahedral topology has one central metal vertex and six peripheral ligand vertices.15  
* **Edge Responsibilities:** Edges define the connectivity. In **stk**, the creation of a bond is the realization of an edge between two vertices.

### **3.2. Building Blocks and Functional Groups**

To populate the topology graph, users define BuildingBlock objects.

* **Functional Group Detection:** **stk** uses SMARTS pattern matching (via RDKit) to automatically identify "functional groups"—the atoms capable of forming bonds. For an amine ligand, the \-NH2 group is detected. For a metal atom, the metal itself is treated as a single-atom functional group.16  
* **Example Construction:**  
  Python  
  \# Conceptual example of stk workflow based on snippets  
  metal \= stk.BuildingBlock('\[Fe+2\]', functional\_groups=)  
  ligand \= stk.BuildingBlock('NCCN', functional\_groups=\[stk.PrimaryAminoFactory()\])  
  complex \= stk.ConstructedMolecule(  
      topology\_graph=stk.metal\_complex.Octahedral(  
          metals=metal,  
          ligands=ligand  
      )  
  )

  This code snippet demonstrates the high level of abstraction: the user does not manually specify coordinates; the Octahedral topology graph handles the placement logic.7

### **3.3. Specialized Metal Complex Topologies**

**stk** includes a dedicated stk.metal\_complex module containing predefined topologies for common coordination geometries.15

* **Supported Geometries:** The library natively supports SquarePlanar, Octahedral, Tetrahedral, Porphyrin, and Paddlewheel geometries.  
* **Isomer Specification:** The topology graph allows for the explicit assignment of building blocks to specific vertices. This enables the deterministic generation of isomers (e.g., placing Ligand A at vertices 0 and 1 for *cis*, or 0 and 5 for *trans*).7  
* **Polydentate Support:** **stk** supports bidentate ligands by mapping a single building block to multiple vertices on the topology graph. The topology enforces the geometric constraint that the donor atoms must align with the target vertices, implicitly handling the "bite angle" requirement (though severe mismatches will result in high-energy structures).6

### **3.4. Optimization and stko (stk optimizer)**

A structure generated by **stk** is a "placement" model—atoms are at the correct vectors, but bond lengths are idealized and unrelaxed. To convert this into a physically chemically valid structure, the **stko** (stk optimizer) library is used. Also **MIT Licensed**, **stko** provides wrappers for external optimization engines.7

* **The "Collapser":** One of the unique tools in **stko** is the Collapser or MCHammer (Monte Carlo Hammer). This algorithm relaxes the long bonds created during the graph assembly (e.g., when a ligand is placed slightly too far from the metal) without exploding the structure, a common issue with standard gradient descent minimizers when starting far from equilibrium.7  
* **Integration:** **stko** seamlessly integrates with **RDKit** (for UFF/MMFF optimization), **GULP** (General Utility Lattice Program, for specialized force fields), and **xTB** (Extended Tight Binding, for semi-empirical quantum mechanics). This allows the user to build a pipeline: stk.Build ![][image8] stko.Collapser ![][image8] stko.XTB.13

### **3.5. Evolutionary Algorithms (EA)**

A distinguishing feature of **stk** is its built-in Evolutionary Algorithm module. This allows for the "Generative Design" of metal complexes.

* **Mechanism:** Users can define a "fitness function" (e.g., cavity size, pore diameter, or energy). The EA then iteratively mutates the population of molecules—swapping ligands, changing metals, or altering the topology—to evolve structures with optimized properties. This moves **stk** beyond a simple builder into the realm of automated discovery.15

## ---

**4\. The Template-Based Specialist: cgbind (MIT License)**

While **stk** is a generalist builder, **cgbind** is a specialized tool developed by the Duarte Group (University of Oxford) specifically for the construction and screening of Metallo-Organic Cages (MOCs). Licensed under the **MIT License**, it is optimized for high-throughput screening tasks where speed and template adherence are paramount.20

### **4.1. The Template Superposition Algorithm**

**cgbind** does not build cages from scratch using topology graphs in the same way **stk** does. Instead, it relies on a library of "Templates" derived from experimental crystal structures (e.g., ![][image9] helicates, ![][image10] tetrahedra, ![][image11] cuboctahedra).20

* **Workflow:**  
  1. **Conformer Generation:** The library generates a set of low-energy conformers for a user-provided organic linker using RDKit.20  
  2. **Motif Extraction:** It identifies the donor atoms (typically N-donors in pyridyl ligands) in the linker.  
  3. **Kabsch Alignment:** The core of **cgbind** is the **Kabsch algorithm**. It calculates the optimal rotation and translation matrix to superimpose the donor atoms of the linker conformer onto the linker vectors of the cage template.21  
  4. **Cost Function:** The quality of the fit is quantified by a cost function based on the Root Mean Square Deviation (RMSD). If the linker is too short, too long, or bent at the wrong angle to bridge the metal centers defined in the template, the cost is high.20

### **4.2. Handling Bite Angles and Geometric Constraints**

**cgbind** explicitly addresses the bite angle problem as a filtering criterion.

* **Pre-Screening:** Before attempting to build a full cage, **cgbind** checks the "donor-donor distance" and the angle of the donor vectors relative to the linker backbone. If these geometric descriptors do not match the requirements of the template (e.g., a specific Pd-Pd distance in a paddlewheel), the ligand is rejected. This makes it an efficient virtual screening tool, capable of processing thousands of ligands in seconds.21  
* **Linker Pre-organization:** The library favors linker conformers that are "pre-organized" for binding (i.e., those where the low-energy solution conformation closely matches the bio-active/cage-active conformation), reducing the entropic penalty of assembly.20

### **4.3. Host-Guest Binding Characterization**

Beyond structure generation, **cgbind** includes unique functionality for characterizing the interior of the generated cages.

* **Cavity Analysis:** It calculates the cavity volume and pore size.  
* **Binding Affinity Prediction:** It uses a simplified force-field approach to estimate the binding affinity of guest molecules within the cage, enabling the automated discovery of hosts for specific substrates (e.g., drug delivery vectors or catalytic substrates).20

### **4.4. Comparison with stk**

* **Speed:** **cgbind** is generally faster for specific cage topologies because it uses rigid-body fitting rather than full graph reconstruction.  
* **Flexibility:** It is less flexible than **stk**. If a user wants to build a novel topology not in the library, they must manually define the template vectors, whereas **stk** allows for programmatic graph definition.  
* **Application:** **cgbind** is the tool of choice for "screening" existing libraries against known cage types. **stk** is the tool of choice for "exploring" new topological spaces.23

## ---

**5\. The Graph-Theoretical Architect: Molassembler (BSD License)**

**Molassembler** represents the most rigorous chemical logic among the analyzed libraries. Developed by the Reiher Group (ETH Zurich) as part of the SCINE (Software for Chemical Interaction Networks) project, it is released under the **BSD-3-Clause license**. It focuses on the fundamental graph theory of molecules, treating stereochemistry as a primary graph property rather than a secondary geometric feature.25

### **5.1. Rigorous Stereocenter Permutation**

Standard tools often struggle to enumerate all stereoisomers of a metal complex (e.g., distinguishing *fac* vs. *mer* isomers of ![][image12]). **Molassembler** solves this using advanced permutation group theory.

* **Stereopermutators:** The library assigns "Stereopermutators" to atoms. For a six-coordinate metal center, the library calculates all unique permutations of the ligands based on their graph identity. It can rigorously determine that an octahedral complex with distinct ligands has exactly ![][image13] stereoisomers, and it can generate the 3D structure for each specific one.27  
* **Polydentate Handling:** The graph representation explicitly understands denticity. It knows that a bidentate ligand forms a cycle with the metal, restricting the possible permutations (e.g., a bidentate ligand can span cis sites but usually not trans sites in an octahedron).28

### **5.2. Distance Geometry (DG) Implementation**

**Molassembler** uses a highly specialized Distance Geometry implementation for 3D generation.

* **4D Smoothing:** Unlike standard 3D embedding, **Molassembler** often embeds in 4 dimensions initially to resolve steric clashes and chiral constraints, then projects back to 3D. This results in high-quality initial structures that rarely require the heavy "collapsing" optimization needed by **stk**.26  
* **Inorganic Focus:** The weighting schemes in its embedding algorithm are tuned to handle the variable bond lengths and angles found in transition metal chemistry, avoiding the "organic bias" of RDKit's ETKDG.29

### **5.3. Integration and C++ Core**

**Molassembler** is a Python library wrapping a high-performance C++ core. This provides significant speed advantages for complex graph operations. It serves as a bridge between high-level Python scripting and low-level quantum chemical calculations (e.g., automated input generation for ORCA or Turbomole).5

## ---

**6\. Hierarchical Systems Interface: mBuild (MIT License)**

**mBuild**, part of the MoSDeF (Molecular Simulation Design Framework) ecosystem, offers a different paradigm: hierarchical, "LEGO-like" assembly. Licensed under the **MIT License**, it is primarily designed for soft matter (polymers, surfaces) but possesses unique capabilities for metal-organic systems.30

### **6.1. The Port-Based Assembly System**

**mBuild** does not use topology graphs (stk) or templates (cgbind). Instead, it uses **Ports**.

* **Concept:** A molecule (or fragment) in **mBuild** has "Ports"—vectors defining attachment points.  
* **Metal Implementation:** To build a metal complex, a user defines a "Metal Hub" compound—a single metal atom with, for example, 6 Ports pointing along the Cartesian axes (![][image14]) for an octahedral geometry.  
* **Ligand Attachment:** Ligands are defined with a single Port on the donor atom. The assembly command force\_overlap(port1, port2) snaps the ligand to the metal by aligning the two ports and translating the ligand.33  
* **Flexibility:** This system is incredibly flexible. It allows for mixed-ligand systems (heteroleptic complexes) to be built simply by snapping different ligands to different ports. It naturally handles hierarchical complexity; a "ligand" can itself be a complex supramolecular assembly built from smaller parts.35

### **6.2. Interfacing with Molecular Dynamics (MD)**

**mBuild**'s strength lies in its integration with **Foyer** (atom-typing) and **GMSO** (General Molecular Simulation Object).

* **Automated Parameterization:** Once a metal complex is built in **mBuild**, it can be passed to Foyer to apply a force field. While standard force fields often lack metal parameters, the MoSDeF ecosystem allows for the easy definition of custom XML-based force fields (e.g., UFF or specific metal-organic parameters) to type the generated structure for MD simulations.36  
* **Soft Matter Integration:** **mBuild** is the tool of choice when the metal complex is part of a larger system—for example, a metal catalyst embedded in a polymer matrix or a MOF node connected to a periodic lattice.38

## ---

**7\. Comparative Analysis and Selection Strategy**

To aid in selection, the libraries are benchmarked against critical technical criteria.

### **7.1. Feature Comparison Matrix**

| Feature | RDKit | stk | cgbind | Molassembler | mBuild |
| :---- | :---- | :---- | :---- | :---- | :---- |
| **License** | BSD-3-Clause | MIT | MIT | BSD-3-Clause | MIT |
| **Core Paradigm** | Cheminformatics / SMILES | Topology Graphs | Template Matching | Graph Theory / Stereochemistry | Hierarchical Ports |
| **Metal Handling** | Poor (Requires Dummy Atoms) | Excellent (Dedicated Classes) | Excellent (Dedicated Templates) | Excellent (Rigorous Permutations) | Good (Manual Port Definition) |
| **Isomer Control** | Low | High (Graph Mapping) | Medium (Template Dependent) | Very High (Enumeration) | High (Manual Placement) |
| **Optimization** | DG / ETKDG / MMFF | via **stko** (GULP, xTB, etc.) | via RDKit / xTB | DG (4D Projection) | via Energy Minimizers |
| **Polydentate Support** | Low | High (Graph Edges) | High (Distance Filters) | High (Cycle Detection) | Medium (Multiple Ports) |
| **Best For...** | Ligand Prep / Scripting | Supramolecular Assembly / Polymers | Metallocage Screening | Rigorous Isomer Generation | Interfacing with MD / Soft Matter |

### **7.2. The Contrast with GPL Tools (Contextual Analysis)**

It is instructive to compare these permissive tools with restricted (GPL) alternatives like **molSimplify** and **epic-mace** to understand the "cost" of the license constraint.

* **molSimplify (GPL):** Offers an "all-in-one" GUI and automated DFT prep. It has a massive built-in database of ligands and metals. Users of permissive tools (stk/cgbind) must often manually define the "building blocks" that molSimplify provides out-of-the-box. The permissive tools are "frameworks" for building, whereas molSimplify is a "finished product".40  
* **epic-mace (GPL):** Specialized for stereoisomer enumeration of octahedral complexes. **Molassembler** is the permissive functional equivalent, offering similar rigor in stereochemistry without the viral license.4

## ---

**8\. Strategic Implementation: Building a Permissive Pipeline**

For a professional research environment requiring a permissive license stack, no single tool functions as a complete "turnkey" solution. Instead, a pipeline approach is required. The following architecture represents the optimal integration of these libraries.

### **8.1. Step 1: Ligand Preparation (RDKit)**

The pipeline begins with **RDKit**.

* **Input:** Ligand SMILES (e.g., c1cc(ccn1)c2ccccn2 for bipyridine).  
* **Processing:**  
  * Canonicalize SMILES.  
  * Generate 3D conformers using AllChem.EmbedMultipleConfs.  
  * **Crucial Step:** Identify the donor atoms (e.g., SMARTS query \`\`).  
  * **Crucial Step:** If using dative bonds, update the graph to use BondType.DATIVE to prevent valence errors in downstream steps.1

### **8.2. Step 2: Complex Assembly (stk or Molassembler)**

* **Scenario A: Supramolecular/Polymeric Assembly (Use stk):**  
  * Initialize stk.BuildingBlock with the RDKit molecule.  
  * Define the metal building block (e.g., \[Pd+2\]).  
  * Select the Topology (e.g., stk.metal\_complex.SquarePlanar).  
  * Map ligands to vertices. **stk** handles the translation and rotation (via Kabsch in the background) to form the bonds.15  
* **Scenario B: Exhaustive Isomer Enumeration (Use Molassembler):**  
  * Convert the RDKit molecule to a **Molassembler** molecule.  
  * Define the complex connectivity.  
  * Ask **Molassembler** to enumerate all stereopermutations.  
  * Generate 3D coordinates for the specific isomer of interest (e.g., *mer-isomer*).26

### **8.3. Step 3: Geometry Optimization (stko)**

The raw structure from Step 2 will have correct connectivity but likely unphysical bond lengths (e.g., stretched bonds in a strained cage).

* **Relaxation:** Pass the structure to **stko**.  
* **Protocol:**  
  1. Run stko.Collapser to resolve gross steric overlaps and long bonds using a Monte Carlo approach.  
  2. Run stko.MMFF (using RDKit's force field) for a quick minimization.  
  3. (Optional) If an external optimizer is available (e.g., GULP), use the stko.GulpUFFOptimizer for metal-specific parameters.13

### **8.4. Step 4: Property Calculation (cgbind)**

If the target is a cage or host-guest system:

* Pass the optimized structure to **cgbind**.  
* Use cgbind.Cage.get\_cavity\_volume() to analyze the pore.  
* Use its automated docking routine to screen for guest binding.20

## ---

**9\. Future Outlook and Emerging Trends**

The domain of open-source inorganic discovery is rapidly evolving.

* **Machine Learning Integration:** **stk**'s evolutionary algorithm demonstrates the shift towards *generative* chemistry. Future iterations will likely integrate Graph Neural Networks (GNNs) directly into the fitness functions to predict stability without expensive quantum calculations.15  
* **Database Integration:** Tools like **cgbind** are increasingly connected to large crystallographic databases (CSD), allowing for the automated extraction of templates rather than manual definition.  
* **Force Field Evolution:** The limitations of UFF/MMFF for metals are being addressed by machine-learning potentials (MLPs). We expect to see **stko** and **mBuild** offering interfaces to MLP engines (like NequIP or Allegro) to provide DFT-level accuracy at force-field speeds within these permissive pipelines.

## **10\. Conclusion**

The requirement for MIT, BSD, or Apache licensing in transition metal complex generation does not force a compromise on capability. While "all-in-one" GPL tools exist, the modular ecosystem of **stk**, **cgbind**, **Molassembler**, and **RDKit** offers a superior level of flexibility and architectural rigor for industrial and advanced research applications.

**stk** provides the topological framework for assembling the unbuildable; **cgbind** offers the speed for screening the massive; **Molassembler** ensures the stereochemical validity of the complex; and **RDKit** remains the indispensable linguistic layer for chemical information. By strategically integrating these libraries, researchers can construct automated, high-throughput, and commercially viable discovery engines for the next generation of catalysts and functional materials.

### ---

**Table 1: License and Capability Summary of Analyzed Libraries**

| Library | License | Primary Domain | Key Class/Module | Optimization Backend |
| :---- | :---- | :---- | :---- | :---- |
| **stk** | MIT | Supramolecular Assembly | stk.metal\_complex | **stko** (RDKit, GULP, xTB) |
| **cgbind** | MIT | Metallocage Screening | cgbind.Cage | RDKit, xTB |
| **Molassembler** | BSD-3 | Stereoisomer Enumeration | scine\_molassembler | Distance Geometry (Internal) |
| **mBuild** | MIT | Hierarchical/Soft Matter | mbuild.Port | Foyer (Force Fields) |
| **RDKit** | BSD-3 | Ligand Cheminformatics | AllChem | ETKDG, MMFF94, UFF |

### **Table 2: Comparison of Alignment Algorithms**

| Method | Used By | Mechanism | Pros | Cons |
| :---- | :---- | :---- | :---- | :---- |
| **Kabsch Algorithm** | **cgbind**, **stk** | Minimizes RMSD between donor atoms and template vectors via rotation matrix. | Fast, robust for rigid ligands. | Fails for flexible ligands if conformer is wrong. |
| **Graph-Based DG** | **Molassembler** | Encodes stereochemistry as graph invariants; projects from 4D to 3D. | rigorous handling of chirality and isomers. | Computationally more expensive than Kabsch. |
| **Port Matching** | **mBuild** | Aligns coordinate vectors of two "Ports" and translates. | Intuitive, flexible "snap-together" logic. | Requires manual definition of ports on all fragments. |

#### **Works cited**

1. Getting Started with the RDKit in Python, accessed February 18, 2026, [https://www.rdkit.org/docs/GettingStartedInPython.html](https://www.rdkit.org/docs/GettingStartedInPython.html)  
2. lukasturcani/stk: A Python library which allows construction and manipulation of complex molecules, as well as automatic molecular design and the creation of molecular databases. \- GitHub, accessed February 18, 2026, [https://github.com/lukasturcani/stk](https://github.com/lukasturcani/stk)  
3. rdkit/ReleaseNotes.md at master \- GitHub, accessed February 18, 2026, [https://github.com/rdkit/rdkit/blob/master/ReleaseNotes.md](https://github.com/rdkit/rdkit/blob/master/ReleaseNotes.md)  
4. EPiCs-group/epic-mace: Python package for the automated ... \- GitHub, accessed February 18, 2026, [https://github.com/EPiCs-group/epic-mace](https://github.com/EPiCs-group/epic-mace)  
5. Identifying Dynamic Metal–Ligand Coordination Modes with Ensemble Learning | Journal of the American Chemical Society, accessed February 18, 2026, [https://pubs.acs.org/doi/10.1021/jacs.5c17169](https://pubs.acs.org/doi/10.1021/jacs.5c17169)  
6. Unlocking the computational design of metal-organic cages \- ResearchGate, accessed February 18, 2026, [https://www.researchgate.net/publication/358870540\_Unlocking\_the\_computational\_design\_of\_metal-organic\_cages](https://www.researchgate.net/publication/358870540_Unlocking_the_computational_design_of_metal-organic_cages)  
7. Unlocking the computational design of metal–organic cages \- PMC, accessed February 18, 2026, [https://pmc.ncbi.nlm.nih.gov/articles/PMC8932387/](https://pmc.ncbi.nlm.nih.gov/articles/PMC8932387/)  
8. The official sources for the RDKit library \- GitHub, accessed February 18, 2026, [https://github.com/rdkit/rdkit](https://github.com/rdkit/rdkit)  
9. RDKit, accessed February 18, 2026, [https://www.rdkit.org/](https://www.rdkit.org/)  
10. RDKit Cookbook — The RDKit 2025.09.5 documentation, accessed February 18, 2026, [https://www.rdkit.org/docs/Cookbook.html](https://www.rdkit.org/docs/Cookbook.html)  
11. stk: A python toolkit for supramolecular assembly \- PMC \- NIH, accessed February 18, 2026, [https://pmc.ncbi.nlm.nih.gov/articles/PMC6585955/](https://pmc.ncbi.nlm.nih.gov/articles/PMC6585955/)  
12. Inverse Design of Metal-Organic Polyhedra through Molecular Fragmentation and Evolutionary Optimisation \- Computational Modelling Group, accessed February 18, 2026, [https://como.ceb.cam.ac.uk/media/preprints/c4e-preprint-340.pdf](https://como.ceb.cam.ac.uk/media/preprints/c4e-preprint-340.pdf)  
13. stko documentation \- Read the Docs, accessed February 18, 2026, [https://stko-docs.readthedocs.io/en/stable/\_autosummary/stko.html](https://stko-docs.readthedocs.io/en/stable/_autosummary/stko.html)  
14. Introduction — stk documentation, accessed February 18, 2026, [https://stk.readthedocs.io/](https://stk.readthedocs.io/)  
15. stk: An extendable Python framework for automated molecular and supramolecular structure assembly and discovery | The Journal of Chemical Physics | AIP Publishing, accessed February 18, 2026, [https://pubs.aip.org/aip/jcp/article/154/21/214102/595430/stk-An-extendable-Python-framework-for-automated](https://pubs.aip.org/aip/jcp/article/154/21/214102/595430/stk-An-extendable-Python-framework-for-automated)  
16. Stk: An Extendable Python Framework for Automated Molecular and Supramolecular Structure Assembly and Discovery \- ChemRxiv, accessed February 18, 2026, [https://chemrxiv.org/doi/pdf/10.26434/chemrxiv.14179022.v2](https://chemrxiv.org/doi/pdf/10.26434/chemrxiv.14179022.v2)  
17. Basic Examples — stk documentation, accessed February 18, 2026, [https://stk.readthedocs.io/en/stable/basic\_examples.html](https://stk.readthedocs.io/en/stable/basic_examples.html)  
18. JelfsMaterialsGroup/stko: A collection of molecular optimisers and property calculators for use with stk. \- GitHub, accessed February 18, 2026, [https://github.com/JelfsMaterialsGroup/stko](https://github.com/JelfsMaterialsGroup/stko)  
19. lukasturcani/basic\_ea: A basic example of an stk evolutionary algorithm for molecular design. \- GitHub, accessed February 18, 2026, [https://github.com/lukasturcani/basic\_ea](https://github.com/lukasturcani/basic_ea)  
20. Computational Modeling of Supramolecular Metallo-organic Cages–Challenges and Opportunities \- PMC, accessed February 18, 2026, [https://pmc.ncbi.nlm.nih.gov/articles/PMC9127791/](https://pmc.ncbi.nlm.nih.gov/articles/PMC9127791/)  
21. cgbind: A Python Module and Web App for Automated Metallocage Construction and Host- Guest Characterization \- Semantic Scholar, accessed February 18, 2026, [https://pdfs.semanticscholar.org/fcf4/91f1f91e0fa83dad5d76d16ce45cd323dab4.pdf](https://pdfs.semanticscholar.org/fcf4/91f1f91e0fa83dad5d76d16ce45cd323dab4.pdf)  
22. stko.RmsdMappedCalculator — stko documentation, accessed February 18, 2026, [https://stko-docs.readthedocs.io/en/stable/\_autosummary/stko.RmsdMappedCalculator.html](https://stko-docs.readthedocs.io/en/stable/_autosummary/stko.RmsdMappedCalculator.html)  
23. Adjacent backbone interactions control self-sorting of chiral heteroleptic Pd3A2B4 isosceles triangles and Pd4A4C4 pseudo-tetrahedra | Request PDF \- ResearchGate, accessed February 18, 2026, [https://www.researchgate.net/publication/396602561\_Adjacent\_backbone\_interactions\_control\_self-sorting\_of\_chiral\_heteroleptic\_Pd3A2B4\_isosceles\_triangles\_and\_Pd4A4C4\_pseudo-tetrahedra](https://www.researchgate.net/publication/396602561_Adjacent_backbone_interactions_control_self-sorting_of_chiral_heteroleptic_Pd3A2B4_isosceles_triangles_and_Pd4A4C4_pseudo-tetrahedra)  
24. Modeling Kinetics and Thermodynamics of Guest Encapsulation into the \[M4L6\]12– Supramolecular Organometallic Cage \- ACS Publications \- American Chemical Society, accessed February 18, 2026, [https://pubs.acs.org/doi/10.1021/acs.jcim.1c00348](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00348)  
25. qcscine/molassembler: Chemoinformatics toolkit with support for inorganic molecules \- GitHub, accessed February 18, 2026, [https://github.com/qcscine/molassembler](https://github.com/qcscine/molassembler)  
26. Molassembler: Molecular Graph Construction, Modification, and Conformer Generation for Inorganic and Organic Molecules | Journal of Chemical Information and Modeling \- ACS Publications, accessed February 18, 2026, [https://pubs.acs.org/doi/10.1021/acs.jcim.0c00503](https://pubs.acs.org/doi/10.1021/acs.jcim.0c00503)  
27. Computational Discovery of Transition-metal Complexes: From High-throughput Screening to Machine Learning | Chemical Reviews \- ACS Publications, accessed February 18, 2026, [https://pubs.acs.org/doi/10.1021/acs.chemrev.1c00347](https://pubs.acs.org/doi/10.1021/acs.chemrev.1c00347)  
28. Graph neural networks for predicting metal–ligand coordination of transition metal complexes \- ChemRxiv, accessed February 18, 2026, [https://chemrxiv.org/doi/pdf/10.26434/chemrxiv-2024-nzk5q](https://chemrxiv.org/doi/pdf/10.26434/chemrxiv-2024-nzk5q)  
29. Graph neural networks for predicting metal–ligand coordination of transition metal complexes | PNAS, accessed February 18, 2026, [https://www.pnas.org/doi/10.1073/pnas.2415658122](https://www.pnas.org/doi/10.1073/pnas.2415658122)  
30. Towards Molecular Simulations that are Transparent, Reproducible, Usable By Others, and Extensible (TRUE) \- PMC, accessed February 18, 2026, [https://pmc.ncbi.nlm.nih.gov/articles/PMC7576934/](https://pmc.ncbi.nlm.nih.gov/articles/PMC7576934/)  
31. Quick Start — mbuild 1.3.0 documentation \- MoSDeF, accessed February 18, 2026, [https://mbuild.mosdef.org/en/stable/getting\_started/quick\_start/quick\_start.html](https://mbuild.mosdef.org/en/stable/getting_started/quick_start/quick_start.html)  
32. Example System — mbuild 1.3.0 documentation \- MoSDeF, accessed February 18, 2026, [https://mbuild.mosdef.org/en/stable/getting\_started/example\_system.html](https://mbuild.mosdef.org/en/stable/getting_started/example_system.html)  
33. Source code for mbuild.compound \- MoSDeF, accessed February 18, 2026, [https://mbuild.mosdef.org/en/0.18.0/\_modules/mbuild/compound.html](https://mbuild.mosdef.org/en/0.18.0/_modules/mbuild/compound.html)  
34. Data Structures — mbuild 0.17.0 documentation, accessed February 18, 2026, [https://mbuild.mosdef.org/en/0.17.1/topic\_guides/data\_structures.html](https://mbuild.mosdef.org/en/0.17.1/topic_guides/data_structures.html)  
35. 364785 Mbuild: A Hierarchical, Component Based Molecule Builder, accessed February 18, 2026, [https://www.researchgate.net/publication/267308275\_364785\_Mbuild\_A\_Hierarchical\_Component\_Based\_Molecule\_Builder](https://www.researchgate.net/publication/267308275_364785_Mbuild_A_Hierarchical_Component_Based_Molecule_Builder)  
36. About \- MoSDeF, accessed February 18, 2026, [https://mosdef.org/pages/about.html](https://mosdef.org/pages/about.html)  
37. 1 High-Throughput Screening of Tribological Properties of Monolayer Films using Molecular Dynamics and Machine Learning Co D. Qu, accessed February 18, 2026, [https://par.nsf.gov/servlets/purl/10332707](https://par.nsf.gov/servlets/purl/10332707)  
38. Elucidating the Mechanisms of Ion Permeation through Sub-Nanometer Graphene Pores: Uncovering Free Energy Barriers via High-Throughput Molecular Simulations \- PMC, accessed February 18, 2026, [https://pmc.ncbi.nlm.nih.gov/articles/PMC12752709/](https://pmc.ncbi.nlm.nih.gov/articles/PMC12752709/)  
39. hoobas Documentation, accessed February 18, 2026, [https://hoobas.readthedocs.io/\_/downloads/en/stable/pdf/](https://hoobas.readthedocs.io/_/downloads/en/stable/pdf/)  
40. Graph neural networks for predicting metal–ligand coordination of transition metal complexes \- PMC, accessed February 18, 2026, [https://pmc.ncbi.nlm.nih.gov/articles/PMC12541316/](https://pmc.ncbi.nlm.nih.gov/articles/PMC12541316/)  
41. User's Manual molSimplify version 1.0 \- Kulik Group, accessed February 18, 2026, [https://hjkgrp.mit.edu/tutorials/2016-12-02-molsimplify-tutorial-2-slab-builder/molSimplify\_v1.pdf](https://hjkgrp.mit.edu/tutorials/2016-12-02-molsimplify-tutorial-2-slab-builder/molSimplify_v1.pdf)

[image1]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABMAAAAZCAYAAADTyxWqAAAB7ElEQVR4AezTXSjeURzA8cdmq60ttdXW3ttLu1h7qbW1tVqrrXa3ra1dbBfi1oUoKSQ33FCKXEi58pK4kAtKFHGHXJAiUiIRIi685e37ffIvevCkXCh0Pv/fOefv/J5zfs95LoRO8O882fGLeQZr5pGfUalMFCIPX2C7zuMX8vEHr5GEYqThFkImMOoOj2wsoQEXkYAH+ATfxxBTkI5x1OAzchAbJLvEwMl3xFaMIBaX4bv7xCks4iaa0I4u+L//iVeDZH7iFSaeoBzJqINHMUkL/V5YhkniAJaxDY94jRgXJFtn4I56iO4wi1iGD1jBBOJwD/2Yhc311s/SzDpwVy958w3f8Qal8DgfiUFz1zcY9GEONsvizuoZbJjMmlQySIXb9QgV9DswCpsf+JiOX4p183iW5Sdz0yjBpsm26AyhEbfxFu9hYZuJNhdaL3fxigmP/5vobjOIrt822SYDv9o14g98xQw86gLR9pDHI5jc+v6l/wK5aMMGwvfMLbsLL6mXsoAXtQjqQjf0lIf3rJtYBC+qd3KQvpshhMLJwp0jHt6357y3tl6T8C4YRzSPGTG5Z8IaJjL+h7uIhzWNIUa0aMlWWTEMvy2P5s9snrGlIexv0ZJ5DTpZUrWrmjiGA1u0ZAcuOmzy9CbbAQAA//+9UsdyAAAABklEQVQDAAqTXTOWT8+sAAAAAElFTkSuQmCC>

[image2]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABsAAAAZCAYAAADAHFVeAAACkElEQVR4AdyVSUhVURjH72seaIbmgWgRTTTQXERRrSJaBRFtiqIioqIWRe0KGmkTRdEiohYt2gRRm6KoTUSjiCIqouKAijigOKH+fhfe9W18z6cXF8r3u/9zzz2e75zv+855o4Jh/Ivb2TjWPgPmwiwYA5HF6Wwss26Gi3AXbsMWsB8JgjidzWfGPZAHp6EazsM8CC1OZx3MWAPNYLsUXQqTILQ4nVUx4zP4ABNhFeRCE4QWp7NwQh7maDc6Gl5DHYQWt7MEs66HBaCjCjTyETXodCXL0avwEG7BNnCM3w7Tvg/HYSachTtwA5aBtoLHRjB3neg+mAahOVHY4LEQLkELvIUJcBQWw3ZwtQ3odXgKJfAeFoEON6Cn4Aq40JfoXvDsIUFU+nZspccQfEMLYTzY34Xuh1+gM0vZxXzh/Tf8A3czB30AB8AoHEFdvEeAZp+zBG/JCnpM+wIY88toJbyBHFgHjfAO2sHwTkeNwmTUHOWjScpoG04kiHbmP36l5wfsAEP1BPUG6EY9qG2ou/+OOh4JdGKeWnkxT0j/Zs4SfF4JO+EQrIXnsATWgN+RwLa7MHy+i8k3n4b9rx3p0Jl5ecEgrxarzINoCD/SVw49oBlCQ+3Evnsz7KLhriwGbw5e+zedOVkBQz7DFLCqTLiJt1joCtydxWMoPbDu3kIwGvcYYFUi6U1nJvAmw1zhQdTL1MS+ol0L2lQeOjCE9bStNJ1aTEbBBdOd3nTmQEPj2fCnwYNraacmfDXTzIafoINr6Dn4AwM2nQ1ksKE1R/8Z7LlDsrdMzqy2k0x7DPwFPoFugkFZJmfms4iZH8EZ+ATRLU47K8vkzKLxsFsESYqz8pAyOJOzlKFDb45cZ70AAAD//4S8nugAAAAGSURBVAMApDt7My9OLmYAAAAASUVORK5CYII=>

[image3]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABsAAAAZCAYAAADAHFVeAAAClklEQVR4AdyVR2hVQRSGbxQ7dhHEtlCwVxQroqILEREEQcSNIqKIqLhRdKcLC+pCLOBCRBciIoiCG0XRjYjdEAhJCKRBCiEhvef7JnkkEPIKXLJION/9z53MnTNz5sy8EdEQ/sUdbBRznwLTYQaMhSwIFmewkYy4Hs7BGbgDe2AMBIsz2ERG3ApNcBdqYR9MhWBxBmtlxJ/wF0znaLQK2iBYnMEaGfEj/IedYPreozUQLM5gYUAe7l0zWg6LYQIEizPYNEbcBZPhLZTCcZgDwfoHc0aLaL0IVtI1dBPYx/8dxL8FR8GBT6HX4QosAAc9gK4Ey90968AXJIocKDg87HwebYBX4Bk5jM6DzTAbzP9l9BEUwjuYCwa0QH7gL4QTMBOcdAkaLBHMWWykZQ18gTxwg21vx98NDmSwWfhO5hNq9f1B14GTfYk+hNdwG96Ak0f6VpbF2zhYBvfhLDyHC1AGL+AfrAbPj4O04JtebwyzoNbTVtGLBeKZ6+I9WGJlfviZlm+wBUzVA3QDdEIOWGGu/iu+/ZHIIEtwLHuD4A5uBnNVS+ni6d+ProLHMB9WgP9HIn1nb/p8FyvP/TTtv21IhsHclyd08j6zyrLxTaEHshg/kQZTaKodmOZoPI9t4KqeonWQ1AzmYLn08vR7v63Fd8PdeIuF18jVWTymcjsNrn4vajZuolYlktwM5t11lW7O0ItzB34RPINK0CbxMIAprMY/BAa1mMyCE6YpuRnMjqbGQ3yD7h5cS7v/hi+n3XPzHTXAJfQ0/IK0zWDpdDa17pE3uucunW8G9EkVzGrzfjvCl/4uHUP9gUQyt1TB3M98hr0HJ+ED+BuFZG6pglk0HnaLIEFB5mF6vkgVrKdXTM/hG6wbAAD//9DuxrAAAAAGSURBVAMAcCF/M/Ol4woAAAAASUVORK5CYII=>

[image4]: <data:image/png;base64,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>

[image5]: <data:image/png;base64,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>

[image6]: <data:image/png;base64,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>

[image7]: <data:image/png;base64,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>

[image8]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABMAAAAYCAYAAAAYl8YPAAAA5ElEQVR4AeyRzwoBURjFEUopG1GysbdgxRNYyMJa3sGbeANbK2s2ypKyVzZ2XkD+lITfmXSlWeh+KZuZzm/O3Ome0zd3ErEfXlGZ/2H+/czSzFwAOfaWZbIs8S6U4EOWsiMNFShDEpwsZVfSa2hCHpwsZQrPuaWgBTkIpDLRZtXzoMNefWIfH0AD4irCTdqTUmEd10+JqezOYgpjDybsVW6BD2EFD5Xh3qqRqMIMlnCCYDK5Lzp4lWwI3iCQZbIiyQxs4QJOlrID6RHs4AFOljJNo6Kza3k9WMpe0bBFZeEz+fbmCQAA//9K63PWAAAABklEQVQDAI6iJzHqA5E3AAAAAElFTkSuQmCC>

[image9]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAC8AAAAZCAYAAAChBHccAAAEFElEQVR4AeyXWahOURSAj3meyTxkJjwYQ2ZCUaRIyBDK8IB4E3kSRXgSDxQyJhIioiRDkswyz2We5+n7Duc6fn/uf/R3dcttff/aZ5+991l7rXXW2bdgkI///hv/r4KX6vlKGNIb+v2gNbokpBPn1uBGNL4v7eZQDHKTmgzoDtFzIt2NvoxFA+KDq3MxEhbAdpgPdSCdlKdzHGyEZWC7PboE5Cb1GTAEVsBumAE+dxA6Y0k1/iwzZ8MWuANVIJ3ni9DfDFrAB9gKGr8a/Qxyk0MMmANn4AnMglEwHTKWVOOdWI2fz3AaKkIpSBW93obOm6Dxl9CvIYk0YHBlOAUP4CskknTGm8eFWcVFzd9ytOPj9Lq5+Zh+U+Q++jIklSZMMLIn0W8gscSNcrLGVqDxAq7DOzASRdGRuDm9dp4O08b0ukI7iRRgcFPQ+BPorBiv4WVY7AY8hJegsW6KZmBEptDYBvY1Qmu4Y2lmLKZLQ0Y/gqvwERJLquctlRrvgqbFK1aMe34E19fgFui5L2hTxneEZsZSj5GWywvop/BXks54X8bIeD1v+Yy8PJinrILi0A70nC8rzURiqdT4c8x6DnExujrGqMb7f2vHjddAw2m+m4OR5zXenJ/HbEvae7Tl0w+Y6aL36ArFhy6m5Uu4AV0LUsV81zDvWZpTS6uRH8+kTxCXXlxoQ210KHHjrSp+kCJj9IgLmzYDGH0QjAgqsIRaLfS6m7SvLD8TYQd0AMfuQRslVI64ni/8XXostalGNqbf6FswaIZilMbSctOo7xIZX4hLF9V4J+kd664PsM534b4fLnPcUtmRayuRnqMZipGzdLqO6+2l1+OC0MyRurQ0/iLayKFCMZrem8nVTojEzZuiZoNpGvUHGl+aK88lE9CeUzTMkHIZ3OPHl3Ml2nRqhe4KnkXciN72S2taacgk7q0HU6st2jLqfJqB3uwcBEEPsNLocaOnY0TPruFeSzgAik7w2sIRZYT9IRpvrR3NlSXRyuH5xAl0BX6oPOccD4LAhw1FTwajcgzdCQaCDkDliBvS0Ln06C0jaVSncu2mjqAtj1avabSlJ1oH6CjTlctAm6rSsHS7AZo/ReNNk2F0eVCSMbR3gXKUn00QLea5xzFxFnLf8wkqFFPCSK7lyu+BaedmPW4Mpy8+N117EWMU08WPoJs2wmod7Tth6oZp48Bsoaf6sNhtMGr90aYLKrEYLY8evvgWE9fxndIZOj2rxusd89ZU8HhrupliUdSSWv+WCTpgH9pNWBB8hhtwY1k13hw+zIOWg5XJ/weW0E562mTKL2LKWSz8v2EpdyzP4Rc9dD8d2RA95bdgHYuZ77KZ9h8ko1tWNavefkbrEKOho7LqedbOW8mm5/PWcp6Wr43/BgAA///u9mcLAAAABklEQVQDAAx2xzOXo+CrAAAAAElFTkSuQmCC>

[image10]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAC8AAAAZCAYAAAChBHccAAAEE0lEQVR4AeyWWchNURTHj3meKVPIPMtY5iHFA0VCQoZ4wQPiTeRJyItHHigvSKJM8UBRhiSZ5ynzPM/T73dydNzvfN17btcn9d3W/6599vjfa6299i4b/Me/UvL/ynmZlq8HkeFg5C/0QFcFSeLYxjRE/UdQ7ggqgWzShA5DQLROpAdTl7NIIN65ER9TwEqwE6wAzUCS1KZyJtgC1gLLfdBVQDZpSYdxYB3YCxYC1x2DzlkyyZ9j5BKwDdwFDUCS5StQ3wF0Bp/BdiD5jeiXIJscpsNScBY8B4vBVLAA5CyZ5B3YkL9v4AyoC6qBTNHqPam8DSR/Gf0OpJFWdK4PToPH4AdIJUnkjePyzOKkxm8tyvF+Wt3YfEa9IfIAfRWklXYM0LOn0O9BaomTcrBk61B4DW6Cj0BPVERH4ua02gUqDBvD6xrlNFKGzu2B5E+iC0Je4jWY7BZ4At4AybopioEemUthB7CuDVri9qWYsxguren9FFwHX0BqybS8qVLyTmhYvGXGuOUn830D3AFa7jvakPGMUMxZWtDTdHkR/QLkJUnkPYwReS1v+oysPJZVNoDKoDfQch5WiqnEVCn584x6BeKidzWMXo3XFynHyUtQdxrvxmBkeckb88sZbUr7hDZ9eoEZLlqPqiKynppRIFOMd4k1pcHUnJla9fws2r4CRc/PpuD6y9C/7504ebOKDREZLeLEDpbEIQbqEVRgCjVbaHU3aV0EyY3no7gLx/k88PfoY6qNSPIZSlv+9b4JQ2MaqiYOL8LutA0AGjOIyJejwkkl7yAJmHddwDw/kHYvLmPcVNmXbyfUchT/kOZ8efC9fJyHzz/EdslfolbPoULRm7Yt4ms3UCTrZWhGk7A3/i4avFtC8tX58F2ia3ynSEyXUh3c58/DaQgYTl35HgR8i7iRmpSd3IkpBnrPO8DFP1ChAVChaM3+lIYCM40W13saRsygfhPoAg4CY9/NdKPshegNrBc8Lxo7JG+unUYHU6KZw/eJE1AVeFH5zjkRBIGLTUDPAZI6ju4HRgMN4AYkd4VvD7J9Isur9eo82nqBo8D0KJn5lMUwtJ7QUIarHjb+9fAR2lYDL0XDUS+F5A2TiTT4UBLTKe8ByjH+tgInQwW+e+wTxyoaDBGtZAp0wU7UuSGtZFZxIz43JlEfH5tUXkMfxc3pbQ3h+XMOU7dR4cZC8nYsBFzMW1cPGT5uQguZVvOZ33m8LD13zqH3NIiH3LUKSt6FzEiGmlYyjMwWYXzmw54xvjoNWe8U7xgzm1FhKi8oedYKxUzwiJJnYx9a16PyEs/Afka6AROHT24zXHijR6mS9oKJceldsZkZD4CHoBjJWq0HHa83fcCZXr0kw4F/g3w4cUn8lZIvCSsnrfETAAD//9VRWhYAAAAGSURBVAMAvtzPMzzuNAQAAAAASUVORK5CYII=>

[image11]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAD0AAAAZCAYAAACCXybJAAAE2UlEQVR4AeyXZagWXRDH17e7u/tVUVRssQsTDAQRsUUQ64ug2KD4wQBRQURBUSxsQVTEwsDEQgzsxkDF7t/vwd27Ptz75blHEK4P89+ZM7vP7Jk5M3POvhUVwd8bp4vKomev9Lc4Xhc0eoEK8A9AQfQrN+qA+PkyyIWl/zHQEMQ25Y5LowtC2U7/hNW2YAxYAcaDb0B+pL4TN+aBCaALqAEKS+Ux0B3MBctAB9AOVAVBKNvpfVgdCpaA8+B78B7IpndQGPkS8GdgFugKJoHC0hwMjASXwBFgYMVU5CCU7bRGTdkHCAfAF+BjkE0/o/gR3AOPwV5wF4Sikhj6FKwH2oeFo4KcfsgrdLoY3NVOP/cJOlf4KdxVvgg/B56AUFQWQ5+DzcB3wMJR2hmtvsvlO3AHHAX3we/gbSAZBOv+DwY6+xv8ELgGQpHZVRxjrrDl9sqd/oqX2a1duQvIpq9Ox8Ex+o3R7wCmn86HdvpvbFs6lsx15OAUOxMbdpXfZ3AaXAbW6S9wnxOVkSUd/QvBej8ONziwIOSW5da5E2tmGiws6Uja4g8MXOmTcLunztjYTG9ruzP6GeAz8B84C+zyIVPQfmHwd2E722nfa717n9u5Udpp6/lrzNi5b8GvgtsgXunhyO7HN+EGwBVxlS0DVAlp04xI79lOthdPrAEbQXOQH7n3/8ONuKdY1wwT+hepGTAY1n5/5LVgA2gAsslt1XlXTN9wgvH4SwTTyr0RMbIbW9uudDUUV8A2IBlpa0+nzQh1oiaXicD92maEmKGWXH2XB59ByO65+R02LBmDvJ9nDDwsIUvJ7LqBxnvt4cqt4CPAIpB2zqC0QefcYXnkRByZvnZlYcrapdWbum5NPRiMBpKHFSdnGei0maFebOLSB/i/2DbDDGnb5mfgzJYWGe3LFwNpkO3aPuNd52KmeFJrimILkCwpd4+PGPhe+08TZMnDUykEF9L+hJhHTsx917NtN9Sec42Mhw+G0Rkupq8rp1HrzWfroddZV9w0t/mhKpBmcmcAcAKmqBO1ZlFFBtzVrcVAGBzEyNq1RNT1i6LI0tLRg8jSZC5DgJlmVjkHO75BMnCW6lbum7GwPNJp66gjKreJE/BKwBMRLNrNZSzwkOCz1pPNzPrfg74KaARcQVhCTi4ZpAQn4jl6JbrlQNKWHzZmiOXlmd+tsCc3dbY33A8ZD0uLkR+BNMXzn49yFXAurrKl6fnBIKDOIx05xdBaaw0X1ooNh2GkY9MQXFUj5oeIz6ThCtj0eCwhXyQSBYIpWhvu3jsMbpbAIpvSUoS0zfxkPzz8JuDRhDw31GfkQWkgXJuWSDlkyaZooM0ue4K6SKczwiu4pJ22/qvzDlPX87SHIFMXVc70If/0s9aVXo3sdmtztAzd4w2yQbGUdNg657EouNN/YtWtwxp1u7BUTDe7el/uWR6D4WaMmYOYExXjX/YWbXpCtHtr09R39zHNj/GMdW7P0ulkgROBB0KQ+7qdexTGFgJryonYEKcz9hPUb2S/we24qHIie4b9Zwr/jr+7tWnjQpUhV9w9fByj7cCPKFgUfKXdy2041t463uCWZs0eRl4AZqdgHTLMmXyPzStt06OzBg2KvcMGbJDdJj3weC+40xmjr/sldHq/7v5m5lcknX4OAAD//zGtGX4AAAAGSURBVAMA7NHrM6TagdgAAAAASUVORK5CYII=>

[image12]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAEIAAAAZCAYAAACFHfjcAAAFOklEQVR4AeyXV6hVVxCGl+m994SEJKRX0kwCyUMKKaBgQ7EriqhgAx9UFH0SRcE3ERQroiKK5UEFewe7YsPeK/Zev29z9uW4PffsfY8K9+Fc5j8ze5VZs2bNzFr3oVD+izxQdkTkhhDKjqjEEa/Q/if4N4fv4U+CQlSDxnfAHyAe/wXyoyCNHmPAz8D1YFUm5//ArHjdQvw3+t8CmSgZEW8yqzEYCKaDAeBtUIheorEFmACGgDbAzWkkYlH6kd7h4FeQtIGmVHqCER5Ad/g0oJ3N4E1y0BZt6sn3uyCVkkZsZEYfMBkcBK8BF4XdQZ76Z7R8BW6AiaA1GAEugGKkY9sx4AOg45M20JxKZxkxGAwDp8AooCNitOJ7OWgLGoFUKmSE4XSNmRuAp17IEa/S9yk4ARy7FX4OpNHjDDCMPaWbyOoxxRCrTM8w4xOgDdvh6oNFpC07kdzf+/BUcmBykCd2i8Z14BHwMsgf9xTf1o4rcDdh5OxDdg6sKH1Nr1G2An4aWCPUgVhlep4ZRuR++A6QTx7e3zScBKYOrDjlb9CR5vcLCOeBHr0EtyCaCohBo19H+BC4+Edwx2kMYlHScGvIHkatBFfBvUSE+izO6vMgUBfRG/z2AB5YV/hCkEpJR7zIDEPOjR1DvgiMkNgRcvNwDu16XYfoiCN8p5FR9DSDjIaj8HtxxMPMN4VdvxbybLAoByPAw+vG9xTgYcKKU9IRnpDG7mLacWDhs6CZInyG2vwYbp6CuWe7Y90UXZWShn1D72bg6VnsriObJkYZYpXIQ/iOGYdBf2DxjTGI7y+B0ZAlXRka7npQ6QgjYje9RoQp4iaMBAucV+tI+nSW+W4k6AiaKiXneucbbXMZ5S1zBm6R1RHJw6ArlXzbeAUfYuQCsCUPU5Fdpy68JchE+UZYH7wljAIdYETIDUE3o+f7otV0MT894QN86zRYpfQePZ1Ab6A+T8mo+JZv602+DRbPzrT7xhgLbwrMddgdZFtNWrwxrFWIFWSkeaUarUaNHal6842wSHr625wJLgOjwuLzF/J64JUKCzrMq8to8NawrRBig2fSWSMBT42mYOrJPYjfETS+F3wJqA9cB1ZB2mybkas93l4VnQgemreJDldHJr0qZW6w+GiQjthLg0bDgqFnGlgbhtoAfAu4kPe27wcXpPkuUrdR4+Y84eQAT9M0sRjbpz4db5hrjzbYn9Rvn89n02uTE3MwArzqfV16cN4Wk+jLpFdj9ew/TPBZ6t3rszc2Tkd4g/iCM2UsQj5tfRSZ43Fk6BxUVJD5qzH9aNGR1h7EiEyHmkimFyy45i8IbtwT9Gms3o9pWwqMOlhE6tWx6naDz9KqU8T/yD6pvS0WI3cEFnFTJU1v0BEWLP9nsBZ4Ff6EAjcMC2v58R3vlefC9fjWYc5bg+y74D/4cyAmN+QT1yruC88TtVjG/Z8j+By31sxANszbw70JYMEUdbyO9tFlaIfcn/M6IHtbrIa7dhe48FrXLm8N107WjmJ6I0fo8YYoc5OiOfIsIPnwGYfghmDBE3ZMPjxBC6v9wo1rbDxGA/1fxD4xnx8Njfvlrul9bxRY2EYzZhmwXWebDnwGHeb4QmjAAK/MMXAdCIvIlEnTGzkiGl0Nfkwh/0Gqgy0ab0SYcnGk0FwSZdJrKJWk/QFMcuNeq4awNcj0nMc63lYWTcSSKJPe6uQIU8Pn8Xi2q0N0gkU67Z3C8KKUSW91coS1xYeQRdAoWMX2dIhVH7FkyqS3Ojmi5J3ej4llR+S8WHZEzhG3AQAA//+Hc9tkAAAABklEQVQDAHEnGEJCKkQXAAAAAElFTkSuQmCC>

[image13]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAAZCAYAAAA8CX6UAAABtElEQVR4AeySyyumYRiHXzPT1ByamebQNLPCitgQYScrKSt2NhayUJT/wY6ytLBQ8gfISsmKECE5xUJyzEKRHHK+rk+v3odEsrDw9bu++3nu9+n33vf9vO+iF/q9bqM0uvwN5QnyWH+CWF9ZFEHyTDb7tGRrGv0jWQet0AttkAOxvrOohHbogCYohMDoksQsNIMm48QssIL3RLXJXwsMQCc0QA9cJitin5Jv/cCqG87BN/4nxjpj8QUGYQ0sILpr5Ju/8fAChmEMNHJWLFOy/WNWu5AyId4z+kzSt68QV2EIfkEx/ADlcH1+4CbmbkXeSjoP5+EEJmEJSsB5ESKH/ySjDE4vgtLQ9vLZyEdiJuyA7RFulKzI67c1+47LNo5wdB1KoQzM7ROv4FZJI00sf46n8SGj7U2R06iauA17EChp5Bes0XRwIoo22I+ClVYRD+FBI6/9Lwe8EefCMpBGM2ROwTaPiIGs6CeZCmgEb6SAqCnhVgusJsB5OWiWoTTyumtJ/4FlqIFcSMpK+kl0wRbck0YO0t5j6jnl508I5JfeR8YvmhBKozDzzN2b0eODe30zugYAAP//5WffyAAAAAZJREFUAwDuDVUz8AYGMQAAAABJRU5ErkJggg==>

[image14]: <data:image/png;base64,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>
