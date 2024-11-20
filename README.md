
# OMAnnotation

Authors: Sade Bates, Yannis Nevers

OMAnnotation is a workflow for combining protein-coding gene structural annotation from different sources, using evolutionary information from other species. It repurposes the OMA orthology software (https://omabrowser.org/standalone/) to form a consensus annotation by combining evidence from multiple annotation methods with existing homology information from the OMA Browser. This strategy is made easy to use through a Python toolkit available in this repository.

## Installation

To use the OMAnnotation toolkit, clone this GitHub repository. Then move to the repository folder and install the associated environment with:
```conda env create -f OMAnnotation_conda.yaml```
Then activate the environment:
```conda activate OMAnnotation_env```
Then, you can try running the toolkit using:
```python OMAnnotation/OMAnnotation.py -h```

## Using OMAnnotation

OMAnnotation is run in four steps.
1. Download precomputed protein All-against-All orthology data from the OMA Database from https://omabrowser.org/oma/export/ and extract the obtained archive.
2. Run the prepare_data step of the OMAnnotation toolkit to prepare the annotation data to be run through OMA Standalone with the ```prepare_data``` module (see ```python OMAnnotation/OMAnnotation.py prepare_data --help``` for usage details).
3. Run the OMA Standalone software
4. Extract the consensus annotation from the OMA Standalone results with the OMAnnotation toolkit using the ```extract_consensus``` module (see ```python OMAnnotation/OMAnnotation.py extract_consensus --help``` for usage details).

These steps are detailed below.

### 1. Exporting data from the OMA Browser

OMAnnotation uses the OMA Standalone orthology inference software, repurposed to combine gene models from different annotation sources into a consensus annotation, using evolutionary information as a tiebreaker. It uses the OMA algorithm to integrate orthology data from other species as external evidence that a given gene exists. Since computing these relations can be computationally expensive, we recommend using precomputed data from the OMA Browser as input.

You can do this through the ```https://omabrowser.org/oma/export/``` page. Select a set of species using the interactive tree or the search bar and click ```submit``` to download the all-against-all archive file. 

Your species selection depends on your species of interest; we recommend selecting 5-20 species that are recognised as having a high-quality annotation and that cover some diversity in the taxonomic range of your species of interest. 

_Note_: the number of species you select is a balance between the compute time you have available and the extra evolutionary information they bring. Selecting ten species is usually optimal.

#### Example step 1: exporting data from the OMA Browser
If you are interested in doing a test run of OMAnnotation, this example demonstrates how to reproduce a consensus annotation of *Drosophila melanogaster*.

First, we will download the OMA orthology data for the following five species using the OMA export AllAll utility at https://omabrowser.org/oma/export/. 
- *Drosophila simulans*
- *Lucilia cuprina*
- *Megaselia scalaris*
- *Aedes aegyptii*
- *Bombyx mori*

Select them in the tree and then click ``submit```. Once the archive is ready, download it, move it to your working directory, and unzip it.
This will now be the working directory for OMAnnotation.

### 2. Preparing the data

The two mandatory inputs needed for OMAnnotation are: 
1. The genome file of the species to be annotated.
2. The GFF files from each source annotation in a single folder. 

With these two inputs and the exported data from OMA Standalone, you can use the OMAnnotation ```prepare_data``` module. 

#### Usage

Required arguments ```-a (--gff_annot_folder)```, ```-g (--genome_file)```, ```-d (--db_folder)```, ```-f (--fasta_folder)```, ```-s (--splice_folder)```
```usage: OMAnnotation.py prepare_data [-h] -a GFF_ANNOT_FOLDER -f FASTA_FOLDER -d DB_FOLDER -s SPLICE_FOLDER -g GENOME_FILE [-t FEATURE_TYPE]```

#### Arguments

 Flag                                                         | Default         | Description                                                                                                                                                 |
|:-------------------------------------------------------------|:----------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------|
| ``-a`` ``--gff_annot_folde``                  |                 |Path to the folder containing the source GFF annotation files.                                                                                                           |
| ``-f`` ``--fasta_folder``  |            |Path to the output folder that will contain FASTA sequences corresponding to each annotation.           |
| ``-d`` ```--db_folder```                |                 | Path to the OMA DB Folder that will be used by OMA Standalone        |
| ``-s`` ``--splice_folder``  |  |Path to the output folder that will contain splice files corresponding to each annotation in which splice variants are explicitely included.                                                      |
| ``-g`` ``--genome_file``                |             | Path to the genome file to be annotated, in the FASTA format.                                                                                                 |
| ``-t`` ``--feature_type``         | None            | Path to a tsv file indicating which feature should be used in GFF file to define gene, transcript or CDS. 4 columns by line in order: filename, gene, transcript, cds. (Optionnal)|

#### Example step 2: preparing the data

Collect the genome file to annotate (here: Droso.fasta) and the GFF files of three different annotations placed in the same directory (here: GFF).

To prepare your annotation data, run the following command:
```python OMAnnotation/OMAnnotation.py prepare_data -g Droso.fasta -a GFF -f FASTA -s Splice -d OMA.2.6.0/DB/```
This script will extract protein sequences according to the GFF files, as well as splicing isoform information if it exists, and feed this information into the OMA database to prepare for the OMA Standalone run.

You will also need to provide OMA with a species tree that shows the relationship between the annotation species and all other species in OMA.
Running the `prepare_data` module writes a tree for the species files you exported from OMA browser into the parameters.drw file, on the ```SpeciesTree``` line. Identify the node where your species should be and add a polytomy containing all the prefixes of your annotation file names, e.g. the annotation 'RNA_seq.gff' is entered as 'RNA_seq'.

In our case, the tree is:
```(BOMMO,(AEDAE,(MEGSC,(DROSI,LUCCU))));```
The polytomy from our example annotation prefixes is ```(rna,homology,abinitio)```. This is added as a sister species of *Drosophila simulans*, as it is the closest species to *Drosophila melanogaster, so the final modified tree becomes:
```(BOMMO,(AEDAE,(MEGSC,((DROSI,(rna,homology,abinitio)),LUCCU))));```

### 3. Running OMA Standalone

After preparing the data as explained above, we can run OMA Standalone.

Refer to OMA Standalone Cheatsheet for a quick overview of this step (https://lab.dessimoz.org/blog/media/2020/04/omastandalone_cheat_sheet.pdf), and the OMA Standalone user guide for more detailed information (https://omabrowser.org/standalone/).

_Note_: the run time for the OMA Standalone step depends on the number of precomputed species you are using and the number of source annotations you are making a consensus from.

#### Example step 3: running OMA Standalone
If you can access a High-Performance Computer (HPC) with a SLURM scheduler, you can copy the OMA job scripts we have included in the main folder of this repository into your OMA folder. If not, you will need to use the information from the OMA Standalone user guide to write your own. Run scripts 1-3 (```oma_part1.sh```, ```oma_part2.sh```, ```oma_part3.sh```) sequentially (it is important to wait for each script to complete before moving on to the next). Check your job logs for any errors after running each script, and see https://lab.dessimoz.org/blog/media/2020/04/omastandalone_cheat_sheet.pdf for the list of output files your OMA run should generate.

### 4. Extract the consensus annotation

Once the OMA Standalone run is complete, the final step is to extract the consensus gene annotation using the resulting OrthoXML. In principle, OMAnnotation will select one gene model from each of what OMA infers as "ancestral Hierarchical Orthologous Groups" for the annotation. These are genes supported by more than one annotation method or only one annotation method and homologs from the other user-selected related species.

#### Usage
Required arguments ```-a (--gff_annot_folder)```, ```-x (--orthoxml)```, ```-st (--species_tree) ```, ```-f (--fasta_folder)```, ```-o (--output_prefix)```

```usage: OMAnnotation.py extract_consensus [-h] -a GFF_ANNOT_FOLDER -f FASTA_FOLDER -x ORTHOXML -st SPECIES_TREE -o OUTPUT_PREFIX [-t FEATURE_TYPE]```

#### Arguments

 Flag                                                         | Default         | Description                                                                                                                                                 |
|:-------------------------------------------------------------|:----------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------|
| ``-a`` ``--gff_annot_folde``                  |                 |Path to the folder containing the source GFF annotation files.                                                                                                           |
| ``-f`` ``--fasta_folder``  |            |Path to the folder that contains FASTA sequences corresponding to each annotation. (Generated in the prepare_data step)         |
| ``-x`` ```--orthoxml```                |                 | Path to the orthoXML file from the OMA Standalone run        |
| ``-st`` ``--species_tree``  |  | Path to the species tree used in the OMA Standalone run, in newick format. A copy can be found in the output folder of OMA as ManualSpeciesTree.nwk                                                      |
| ``-o` ``	``--output_prefix``                |             | Output file path and prefix. Two output files will be created: a GFF annotation file and a FASTA file.'format.                                                                                                 |
| ``-t`` ``--feature_type``         | None            | Path to a tsv file indicating which feature should be used in GFF file to define gene, transcript or CDS. 4 columns by line in order: filename, gene, transcript, cds. (Optionnal)|

Once this is done, the resulting consensus annotation and protein FASTA file should be available at the provided output path.

#### Example step 4: extract the consensus annotation

Run this command, substituting the paths corresponding to your directory structure if necessary:

```python OMAnnotation/OMAnnotation.py extract_consensus -g Droso.fasta -a GFF -f FASTA -st OMA.2.6.0/ManualSpeciesTree.nwk -x OMA.2.6.0/Orthoxml.oxml  -o ConsensusAnnotation```

The resulting files should be available as ConsensusAnnotation.gff and ConsensusAnnotation.fa
