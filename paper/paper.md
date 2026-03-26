<img width="619" height="392" alt="image" src="https://github.com/user-attachments/assets/290af12b-1226-4d71-ac77-b6503cc40d90" />1) Overview

Anubis: An Interactive Platform for Collagen Data of Ancient Egyptian Species

Paper Authors
1. El-Azab, Hagar
2. M.Hassan, Nawal.

Paper Author Roles and Affiliations
1. First author, and Corresponding author, Department of Biotechnology, Faculty of Agriculture, Cairo University, Giza, Egypt.
2. Undergraduate Studies at the Department of Biotechnology, Faculty of Agriculture,
Ain Shams University, Cairo, Egypt.

Abstract
Zooarchaeology by Mass Spectrometry (ZooMS) is a method of identifying fragmentary archaeological remains through collagen fingerprinting. However, researchers lack specialized databases for unique peptides and often face complex data retrieval from UniProt. We present Anubis, a Shiny application designed for the analysis of ancient Egyptian bone, leather, and skin samples. Anubis integrates Swiss-Prot and TrEMBL sequences, employing a weighted scoring model to ensure reliability while maximizing coverage. The application utilizes automated in silico trypsin digestion to identify unique peptide markers for key species, including cattle and sheep. This tool enhances taxonomic precision and facilitates the workflow for archaeological proteomics and Egyptological research.

Keywords
Anubis, ZooMS, Paleoproteomic profiling, tandem mass spectrometry (LC-MS/MS), clothomics, ancient Egyptian.


Introduction 
The advent of collagen fingerprinting, formally designated as Zooarchaeology by Mass Spectrometry (ZooMS), has evolved into a highly prevalent methodology for enhancing the precision of taxonomic data derived from fragmentary bone and skin tissues [1].  Bone is composed mainly of type I collagen, a three-stranded molecule with a triple helical structure The collagen triple-helix structure [2]. Tetrapods generally possess two identical a1(I) chains and one a2(I) chain [3]. A repeated glycine-X-Y amino acid motif is present in both strands, where amino acids X and Y are typically proline and hydroxyproline, respectively [4]. This structural element displays notable conservation across multiple species, particularly within the a1(I) chain. These developments have facilitated the creation of ZooMS [5], a collagen peptide mass fingerprinting technique capable of extracting taxonomic data from fragmentary bones, frequently down to the genus level.
The major basis of ancient Egyptian beliefs and rituals about life after death was the hope for an ongoing existence in the next life [6]. Preserving the corpse was essential, since it served as the place of residence for the ka and the ba, both components of human existence. This ensured the transformation of the earthly body into an eternal one [7]. Mummification was absolutely essential for the incredible transformation of the corpse into a mummy, and the deceased into a glorified spirit [8]. Anubis is identified as the Golden Jackal (Canis aureus) the inventor of mummification, first performing the rite on the god Osiris. In addition, priests known as stm priests wore jackal masks to become Anubis during the wrapping and the "Opening of the Mouth" ceremony [9].
In consideration of collagen fingerprinting in ZooMS, a Shiny application was developed. This application, named Anubis, is intended for researchers interested in studying ancient Egyptian samples from bone, leather, or skin. These samples are based on unique peptides associated with species utilized by the ancient Egyptians.

Application overview 



Data collection module
The application consists of three main modules, the first module is the data collection module in which collagen associated sequences are downloaded from the uniprot [10, 11] database including both swissprot [12], and TREMBL databases [13]. For swissprot, the data is manually curated and reviewed, while for TREMBL a score is assigned to each protein based on the protein existence values (Protein uncertain; Protein predicted; Protein inferred from homology; Experimental evidence at transcript level and Experimental evidence at protein level) as defined below

S =L/5
In this model, L denotes the evidence assigned values, with protein uncertainty set at 1. The values assigned to experimental evidence at the transcript level and experimental evidence at the protein level are 5. Further filtration is performed based on the level of function and expression to ensure sequence relatedness. This filtration is performed based on expression/function in bone, leather, or skin.
Data preprocessing module
Subsequently, the sequences that have been downloaded from the data collection step are transferred to the data preprocessing module. In this module, the sequences undergo in silico digestion using trypsin, along with various parameters for missed cleavages. Only peptides that are unique to a particular species are subsequently subjected to further analysis. Statistics for all species are calculated to ensure the availability of peptides for the species utilized by the ancient Egyptians. 
Data visualization module
Statistics for species, proteins and peptides are augmented into a shiny web application [14]. The application, developed using the Shiny framework, offers researchers a visual representation of the protein and peptide statistics for each species. This facilitates the analysis of collagen fingerprinting data for ZooMS research. The multi-tab application provides visualization for proteins, peptides per species with their corresponding scores derived from the data collection module. The complete process for all the modules is shown in Figure 1 below. 
<img width="619" height="392" alt="image" src="https://github.com/user-attachments/assets/05bddd8c-076a-47eb-a854-12aad32c6ad9" />

Figure 1: Workflow of Anubis database. The data is passed through the main modules: data collection, data preprocessing and data visualization. 

Statement of Need
Protein databases are utilized extensively, with notable examples including Uniprot, Swissprot, and TREMBL. To the best of our knowledge, there are no databases for unique peptides associated with collagen used in ZooMS fingerprinting due to its stability and preservation. The process of downloading related sequences associated with the specified tissues is a complex task, primarily due to the indirect nature of the information provided by the Uniprot database concerning the existence and expression of proteins in specific tissues. The process of trypsin digestion can be executed through the employment of varying parameters, contingent upon the expertise of the researcher. Anubis database provides a resource for researchers in the ZooMS field. 
Anubis provides peptides that are unique to goat, cattle, bovine, and sheep species. These peptides were utilized by the ancient Egyptians. Furthermore, the method of data collection varies among researchers due to the preference of most researchers to utilize the Swissprot database. This practice ensures the reliability of the data, yet it can result in a limited number of peptides for each species, thereby potentially impacting the precision of species classification. Aubis utilized both Swissprot and TREMBL with assigning a score for TREMBL to make the best use of the available resources in the Uniprot database. 

Implementation and architecture
The preprocessing is performed using Python programming language (RRID: SCR_008394, version: 3.12.12) and Jupyter Notebook (RRID: SCR_018315, version: 7.5.4) [15]. The web application is developed using the shiny platform (RRID:SCR_001626, version: 1.13.0). 
The main page contains an overview of the different modules, beginning with an explanation of the mummification and its relation to collagen fingerprinting. The left panel displays Aubis performing the mummification on one of the ancient kings. As illustrated in the right panel, the statistical data from the databases is presented in a visual format for interpretation. This analysis indicates that 17 species are linked to 67 distinct proteins and 24,804 peptides, of which 13,840 are unique (representing 55.79% of the total) as explained in Figure 2.
 <img width="674" height="375" alt="image" src="https://github.com/user-attachments/assets/ab0f0e0f-1403-48e1-a245-ebe381a5779d" />

Figure 2: Homepage of the Anubis database. The home page contains a series of tabs, each displaying visualizations and statistics of the database content.
(2)Availability  
The source code can be accessed at https://github.com/Hagarelazab/Anubis. The application can be executed in two distinct environments: either on a local host via RStudio or on a standalone server. A server-based version of the web application is currently operational at https://anubisdb.shinyapps.io/anubis_/. A comprehensive user guide is accessible at the web application's repository, providing links to data and processing codes. Experts are invited to utilize the software, report any errors they may encounter, and suggest, write, or edit functions. These contributions can be made via the GitHub repository or by contacting the author directly.
Operating system
Windows, Linux and Mac
Programming language
R version 4.5.3.
Additional system requirements
Minimum ram 4GB
Dependencies
RStudio 2026.01.1+403
Software location:
Archive (e.g. institutional repository, general repository) (required – please see instructions on journal website for depositing archive copy of software in a suitable repository) 
Name: Anubis
Persistent identifier: 10.5281/zenodo.18998816.
Licence: GNU GENERAL PUBLIC LICENSE
Publisher: Hagar El-Azab
Version published: 1.0
Date published: 13/03/2026
Code repository (e.g. SourceForge, GitHub etc.)
Name: Anubis
Identifier: https://github.com/Hagarelazab/Anubis 
Licence: GNU GENERAL PUBLIC LICENSE
Date published: 09/03/2026
Language
English

(3) Reuse potential 
The objective of this project is to provide other researchers with a user-friendly web application customized for ZooMS collagen fingerprinting task.  This application ensures accessibility and facilitates repeatability, thereby guaranteeing reusability and citability of the software.

Competing interests 
The authors declare that they have no competing interests.
References 
1. 	Oldfield, E-M, Dunstan, M S, Chowdhury, M P, Slimak, L, and Buckley, M 2025 AutoZooMS: Integrating robotics into high-throughput ZooMS for the species identification of palaeontological remains at Grotte Mandrin, France. Archaeological and Anthropological Sciences, 17(1): 12. [Online; accessed 10-March-2026]
2. 	Brodsky, B and Ramshaw, J A 1997 The collagen triple-helix structure. Matrix Biology: Journal of the International Society for Matrix Biology, 15(8-9): 545–554. [Online; accessed 10-March-2026]
3. 	Karsdal, M 2019 Biochemistry of collagens, laminins and elastin: Structure, function and biomarkers. Academic Press. [Online; accessed 10-March-2026]
4. 	Chu, M L, Wet, W de, Bernard, M, Ding, J F, Morabito, M, Myers, J, Williams, C, and Ramirez, F 1984 Human pro alpha 1(I) collagen gene structure reveals evolutionary conservation of a pattern of introns and exons. Nature, 310(5975): 337–340. [Online; accessed 10-March-2026]
5. 	Buckley, M, Collins, M, Thomas-Oates, J, and Wilson, J C 2009 Species identification by analysis of bone collagen using matrix-assisted laser desorption/ionisation time-of-flight mass spectrometry: Species identification of bone collagen using MALDI-TOF-MS. Rapid Communications in Mass Spectrometry: RCM, 23(23): 3843–3854. [Online; accessed 10-March-2026]
6. 	Smith, Mark. "Traversing eternity." Texts for the Afterlife from Ptolemaic and Roman Egypt (2009). [Online; accessed 10-March-2026]
7. 	Taylor, J H 2001 Death and the afterlife in ancient Egypt. British Museum Press. [Online; accessed 10-March-2026]
8. 	Zesch, S, Panzer, S, Paladin, A, Sutherland, M L, Lindauer, S, Friedrich, R, Pommerening, T, Zink, A, and Rosendahl, W 2024 The multifaceted nature of Egyptian mummification: Paleoradiological insights into child mummies. PloS One, 19(12): e0316018. [Online; accessed 10-March-2026]
9. 	Osborn, D J and Osbornová, J 2024 The mammals of ancient Egypt. Casemate. [Online; accessed 10-March-2026]
10. 	UniProt Consortium 2015 UniProt: a hub for protein information. Nucleic Acids Research, 43(Database issue): D204–12. [Online; accessed 10-March-2026]
11. 	Butler, D 2002 NIH pledges cash for global protein database. Nature, 419(6903): 101. [Online; accessed 10-March-2026]
12. 	Boeckmann, B, Bairoch, A, Apweiler, R, Blatter, M-C, Estreicher, A, Gasteiger, E, Martin, M J, Michoud, K, O’Donovan, C, Phan, I, Pilbout, S, and Schneider, M 2003 The SWISS-PROT protein knowledgebase and its supplement TrEMBL in 2003. Nucleic Acids Research, 31(1): 365–370. [Online; accessed 10-March-2026]
13. 	O’Donovan, C, Martin, M J, Gattiker, A, Gasteiger, E, Bairoch, A, and Apweiler, R 2002 High-quality protein knowledge resource: SWISS-PROT and TrEMBL. Briefings in Bioinformatics, 3(3): 275–284. [Online; accessed 10-March-2026]
14. 	Chang, W, Cheng, J, Allaire, J J, Sievert, C, Schloerke, B, Aden-Buie, G, Xie, Y, Allen, J, McPherson, J, Dipert, A, and Borges, B 2012 shiny: Web Application Framework for R. DOI: https://doi.org/10.32614/cran.package.shiny
15. 	Sanner, M F 1999 Python: a programming language for software integration and development. Journal of Molecular Graphics & Modelling, 17(1): 57–61. [Online; accessed 10-March-2026]

