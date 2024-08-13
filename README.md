
# OMAnnotation

Authors: Sade Bates, Yannis Nevers

OMAnnotation is a workflow for combining protein-coding gene structural annotation from different sources, using evolutionary information from other species. It uses the OMA orthology annotation software, repurposed to combineevidence from multiple annotation methods with existing homology information from the OMA Browser. This strategy is made easy to use through a Python toolkit available in this repository.
## Installation

To use the OMAnnotation toolkit, clone this GitHub repository. Then move to the repository folder and install the associated environment with:
```conda env create -f OMAnnotation_conda.yaml```
Then activate the environment:
```conda activate OMAnnotation_env```
Then, you can try running the toolkit using:
```python OMAnnotation/OMAnnotation.py -h```

## Using OMAnnotation

OMAnnotation is run in four steps.
1. Downloading precomputed protein All-against-All from the OMA Database from https://omabrowser.org/oma/export/ and extracting the obtained archive,
2. Running the prepare_data step of the OMAnnotation toolkit to prepare the annotation data to be run through OMAStandalone with the ```prepare_data``` module (```python OMAnnotation/OMAnnotation.py prepare_data --help```)
3. Run the OMA Standalone software
4. Extracting the consensus annotation from the OMA Standalone results with the OMAnnotation toolkit with the ```extract_consensus``` module (```python OMAnnotation/OMAnnotation.py extract_consensus --help```)

These different steps are detailed below.

### Exporting data from the OMA Browser

OMAnnotation makes use of the OMAStandalone orthology inference software, repurposed to combine gene models from different annotation type, using evolutionary information as tiebreaker. It can integrate data from other species, as external evidence that a given gene exists. Since computing these relations can be computationally expensive, we recommand using precomputed data from the OMABrowser as input.
You can do it through the ```https://omabrowser.org/oma/export/``` page. Select your species of interest, and after a few minutes it will download an archive. 
The species selection depends of your species of interest, we recommend selecting between 5-20 species that are recognized as having high quality annotation and cover some diversity in the taxonomic range of your species of interest.

#### Example
If you are interested in doing a test run of OMAnnotation, we will give indication to be able to reproduce a consensus annotation of *Drosophila melanogaster*.
First, we download OMA data for four species: 
- *Drosophila simulans*
- *Lucilia cuprina*
- *Megaselia scalaris*
- *Aedes aegyptii*
- *Bombyx mori*

Select them in the tree and then submit. When it is ready, download the archive.
Once it is done, move the archive to your working directory, and unzip it.
This will be  the working directory for OMAnnotation.

### Preparing the data.

The two mandatory inputs needed for OMAnnotation are: 
1. the genome file of the species to annotation 
2. 2. the GFF files of the composite annotation in a shared folder. 

With these two inputs and the exported data from OMA Standalone, you can use the OMAnnotation ```prepare_data``` module. 

### Usage

Required arguments ```-a (--gff_annot_folder)```, ```-g (--genome_file)```, ```-d (--db_folder)```, ```-f (--fasta_folder)```, ```-s (--splice_folder)```
```usage: OMAnnotation.py prepare_data [-h] -a GFF_ANNOT_FOLDER -f FASTA_FOLDER -d DB_FOLDER -s SPLICE_FOLDER -g GENOME_FILE [-t FEATURE_TYPE]```

### Arguments

 Flag                                                         | Default         | Description                                                                                                                                                 |
|:-------------------------------------------------------------|:----------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------|
| ``-a`` ``--gff_annot_folde``                  |                 |Path to the folder containing the source GFF annotation files.                                                                                                           |
| ``-f`` ``--fasta_folder``  |            |Path to the output folder that will contain FASTA sequences corresponding to each annotation.           |
| ``-d`` ```--db_folder```                |                 | Path to the OMA DB Folder that will be used by OMA Standalone        |
| ``-s`` ``--splice_folder``  |  |Path to the output folder that will contain splice files corresponding to each annotation in which splice variants are explicitely included.                                                      |
| ``-g`` ``--genome_file``                |             | Path to the genome file to be annotated, in the FASTA format.                                                                                                 |
| ``-t`` ``--feature_type``         | None            | Path to a tsv file indicating which feature should be used in GFF file to define gene, transcript or CDS. 4 columns by line in order: filename, gene, transcript, cds. (Optionnal)|

### Example

Collect the genome file to annotate (here: Droso.fasta) and the GFF files of three different annotation, placed in the same directory (here: GFF).
To prepare the data, run the following command:
```python OMAnnotation/OMAnnotation.py prepare_data -g Droso.fasta -a GFF -f FASTA -s Splice -d OMA.2.6.0/DB/```

This script will extract protein sequences corresponding to the GFF file, as well as splicing isoform information if it exists, and feed this information into the OMA database to prepare for the OMA Standalone run.

For this, you will also need OMA with a species tree that indicate the relations between the species we are annotating to all other species in OMA.
Open the parameters.drw file, and find the ```SpeciesTree``` line. 
Here, the tree of the downloaded species is found in the newick file. Identify the node where your species should be and add a polytomy containing all the prefix of your annotation files as prefix in this place.
In our case, the tree is:
```(BOMMO,(AEDAE,(MEGSC,(DROSI,LUCCU))));```
We add our annotation prefix as polytomy ```(rna,homology,abinitio)``` as a sister species of *Drosophila simulans*, as it is the closest species from *Drosophila melanogaster. The modfied tree is :
```(BOMMO,(AEDAE,(MEGSC,((DROSI,(rna,homology,abinitio)),LUCCU))));```

## Running OMA Standalone

After preparing the data as explained above, we can run OMA Standalone.

Refer to the OMA Standaone user guide for more information : https://omabrowser.org/standalone/

Note that this take may take time depending of the number of species you are using and the number of annotation you are making a consensus for.

## Extract consensus

Once the OMA Standalone run is complete, the last step is to extract the consensus gene annotation using the resulting OrthoXML. In principle, OMAnnotation will select one gene model by what OMA infer as "ancestral Hierarchical Orthologous Groups" for the annotation: genes that are supported by either two or more annotation method, or by only one annotation method and homologs from other species. 


### Usage
Required arguments ```-a (--gff_annot_folder)```, ```-x (--orthoxml)```, ```-st (--species_tree) ```, ```-f (--fasta_folder)```, ```-o (--output_prefix)```

```usage: OMAnnotation.py extract_consensus [-h] -a GFF_ANNOT_FOLDER -f FASTA_FOLDER -x ORTHOXML -st SPECIES_TREE -o OUTPUT_PREFIX [-t FEATURE_TYPE]```

### Arguments

 Flag                                                         | Default         | Description                                                                                                                                                 |
|:-------------------------------------------------------------|:----------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------|
| ``-a`` ``--gff_annot_folde``                  |                 |Path to the folder containing the source GFF annotation files.                                                                                                           |
| ``-f`` ``--fasta_folder``  |            |Path to the folder that contains FASTA sequences corresponding to each annotation. (Generated in the prepare_data step)         |
| ``-x`` ```--orthoxml```                |                 | Path to the orthoXML file from the OMA Standalone run        |
| ``-st`` ``--species_tree``  |  | Path to the species tree used in the OMA Standalone run, in newick format. A copy can be found in the output folder of OMA as ManualSpeciesTree.nwk                                                      |
| ``-o` ``	``--output_prefix``                |             | Output file path and prefix. Two output files will be created: a GFF annotation file and a FASTA file.'format.                                                                                                 |
| ``-t`` ``--feature_type``         | None            | Path to a tsv file indicating which feature should be used in GFF file to define gene, transcript or CDS. 4 columns by line in order: filename, gene, transcript, cds. (Optionnal)|

Once this is done, the resulting consensus annotation and protein FASTA file should be available at the provided output path.

### Example

Run this command, replacing if needed by the path corresponding to your directory structure

```python OMAnnotation/OMAnnotation.py extract_consensus-g Droso.fasta -a GFF -f FASTA -st OMA.2.6.0/ManualSpeciesTree.nwk -x OMA.2.6.0/Orthoxml.oxml  -o ConsensusAnnotation```

Your resulting files should be available as ConsensusAnnotation.gff and ConsensusAnnotation.fa