import os
from pathlib import Path
from collections import OrderedDict
from random import sample
from typing import Optional, List, Iterable

import fire
import humanize
import numpy as np
from rich.console import Console
from rich.table import Table
from rich.progress import Progress
import tables as tb
from textdistance import gotoh, smith_waterman
from pyexcel_ods3 import save_data as save_ods


from ppi_origami.sources import (
    string,
    uniprot,
    uniref,
    biogrid,
    hippie,
    common,
    rapppid,
    dscript,
    pipr,
    deepppi,
    oma,
    intrepppid
)
from ppi_origami import analysis
from ppi_origami.sources import utils


class Process(object):
    @staticmethod
    def uniref(raw_folder: Path, processed_folder: Path, threshold: int) -> None:
        """
        Process `UniRef <https://www.uniprot.org/uniref>`_ dataset. Namely, parse the UniRef XML file into a few new formats useful to PPI Origami:

        1. A CSV file mapping UniRef identifiers to UPKB accessions (``uniref{threshold}_members_upkb.csv.gz``)
        2. A CSV file mapping UniRef identifiers to UniParc identifiers (``uniref{threshold}_members_uniparc.csv.gz``)
        3. A CSV file mapping UniRef identifiers to their amino acid sequences (``uniref{threshold}_sequences.csv.gz``)
        4. A LevelDB databse with UPKB accessions keys and UniRef identifier values (``uniref{threshold}_members_upkb.leveldb``)
        5. A LevelDB databse with UniParc identifier keys and UniRef identifier values (``uniref{threshold}_members_uniparc.leveldb``)
        6. A LevelDB databse with UniRef identifier keys and values that correspond to the amino acid sequence of the representative protein of the cluster (``uniref{threshold}_sequences.leveldb``)

        The files generated are particularly useful for identifying whether to proteins, as identified by their UPKB accessions, belong to the same UniRef cluster. PPI Origami uses this to identify similar proteins.

        PPI Origami will parse the UniRef XML file by streaming it, avoiding having to load the whole file into memory.

        The CSV files are simple to use, but the LevelDB databases are very fast without requiring you to load the whole database to memory.

        You can call this function from the CLI using::
        
            ppi_origami process uniref RAW_FOLDER PROCESSED_FOLDER THRESHOLD

        :param pathlib.Path raw_folder: The folder datasets have been downloaded to.
        :param pathlib.Path processed_folder: The folder to output processed data to.
        :param int threshold: The UniRef identity threshold. Must be one of 50, 90, 100.
        :return: None
        """
        console = Console()

        if threshold not in [50, 90, 100]:
            console.log("ERROR: Threshold must be one of 50, 90, or 100.")
            return None

        console.log(f'Running the "uniref{threshold}" processing task.')
        uniref.process(Path(raw_folder), Path(processed_folder), threshold, console)

    @staticmethod
    def string_upkb(
            raw_folder: Path, processed_folder: Path, version: str = "12.0", taxon: Optional[int] = None
    ) -> None:
        """
        Add a UniprotKB accession column to STRING rows.

        This can be a long process, as a STRING ID to UniprotKB accession map must be built.

        Any protein pair with a STRING ID for which PPI Origami can't find the corresponding UniprotKB accession is
        omitted from the dataset.

        You can call this function from the CLI using::
        
            ppi_origami process string_upkb RAW_FOLDER PROCESSED_FOLDER VERSION TAXON

        :param version: The version of STRING DB to process.
        :param raw_folder: The folder datasets have been downloaded to.
        :param processed_folder: The folder to output processed data.
        :param taxon: The NCBI taxon ID of the organism whose links you wish to download. Omit for all organisms.
        :return: None
        """
        console = Console()
        console.log('Running the "string_upkb" processing task.')
        string.process_string_upkb(
            Path(raw_folder), Path(processed_folder), taxon, version, console
        )

    @staticmethod
    def string_uniref(
            processed_folder: Path, threshold: int, version: str = "12.0", taxon: Optional[int] = None
    ) -> None:
        """
        Add a UniRef ID column to STRING rows.

        You can call this function from the CLI using::
        
            ppi_origami process string_uniref PROCESSED_FOLDER THRESHOLD --version 12.0 --taxon 9606

        :param version: The version of STRING DB to process.
        :param processed_folder: The folder to output processed data.
        :param threshold: The UniRef identity threshold. Must be one of 50, 90, 100.
        :param taxon: The NCBI taxon ID of the organism whose links you wish to download. Omit for all organisms.
        :return: None
        """
        console = Console()
        console.log('Running the "string_uniref" processing task.')
        string.process_string_uniref(
            Path(processed_folder), threshold, taxon, version, console
        )

    @staticmethod
    def uniprot_id_mapping(raw_folder: Path, processed_folder: Path) -> None:
        """
        Process UniProt ID mappings. This creates a new CSV file mapping each UPKB accession code to a different identifier. These files will have filenames of the format ``uniprot_idmapping_{id_type}.csv.gz``.

        You can call this function from the CLI using::
        
            ppi_origami process uniprot_id_mapping RAW_FOLDER PROCESSED_FOLDER

        :param pathlib.Path raw_folder: The folder to download the dataset to.
        :param pathlib.Path processed_folder: The folder to output processed data.
        :return: None
        """
        uniprot.process("id_mapping", Path(raw_folder), Path(processed_folder))

    @staticmethod
    def common_format(
            processed_folder: Path, source: str, taxon: Optional[int], identifier: str, version: str
    ) -> None:
        """
        Converts raw source dataset into a "common" format from which PPI Origami can further transform.

        Currently, only the ``string_links_detailed`` and ``dscript`` dataset are supported.

        The identifier argument specifies which identifier to use when referring to proteins in the file. Must be one of:

        1. ``upkb`` -- The UniProtKB accession
        2. ``uniref50`` -- The UniRef50 identifier
        3. ``uniref90`` -- The UniRef90 identifier
        4. ``uniref100`` -- The UniRef100 identifier

        You'll have had to have added a column with the identifier you selected using the appropriate PPI Origami command/function.

        For instance, if you specify ``upkb`` as the identifier, and ``string_links_detailed`` as the source, then, you'll have had to run the PPI Origami command ``string_upkb``.

        You can call this function from the CLI using::
        
            ppi_origami process common_format PROCESSED_FOLDER SOURCE TAXON IDENTIFIER VERSION

        :param pathlib.Path processed_folder: The path to the processed folder.
        :param str source: The name of the source. Must be one of ``string_links_detailed`` or ``dscript``.
        :param Optional[int] taxon: The NCBI Taxon number of the databasae to be converted. Set to ``None`` if there is no specific organism.
        :param str identifier: The identifier to use to refer to proteins. Must be one of ``upkb``, ``uniref50``, ``uniref90``, or ``uniref100``.
        :param int version: The version of the source to use.
        :return: None
        """
        common.to_common(Path(processed_folder), source, taxon, identifier, version)

    @staticmethod
    def hippie_upkb(raw_folder: Path, processed_folder: Path) -> None:
        """
        Add a UPKB acession column to the HIPPIE dataset.

        You can call this function from the CLI using::
        
            ppi_origami process hippie_upkb RAW_FOLDER PROCESSED_FOLDER

        :param pathlib.Path raw_folder: The raw folder where datasets are downloaded to.
        :param pathlib.Path processed_folder: The processed folder, where you will deposit the new dataset.
        :return: None
        """
        console = Console()
        hippie.join_upkb(Path(raw_folder), Path(processed_folder), console)

    @staticmethod
    def common_to_rapppid(
            processed_folder: Path,
            common_path: Path,
            c_types: List[int],
            train_proportion: float = 0.8,
            val_proportion: float = 0.1,
            test_proportion: float = 0.1,
            neg_proportion: float = 1,
            uniref_threshold: int = 90,
            score_key: Optional[str] = None,
            score_threshold: Optional[str] = None,
            preloaded_protein_splits_path: Optional[Path] = None,
            seed: int = 8675309,
            trim_unseen_proteins: bool = False,
            negatives_path: Optional[Path] = None,
            taxon: Optional[int] = None,
            weighted_random: bool = False,
            scramble_proteins: bool = False,
            exclude_preloaded_from_neg: bool = True
    ) -> None:
        """
        Convert a dataset in the PPI Origami "common" format to the RAPPPID HDF5 format.

        You can call this function from the CLI using::
        
            ppi_origami process common_to_rapppid PROCESSED_FOLDER COMMON_PATH C_TYPES --train_proportion 0.8 --val_proportion 0.1 \\
                --test_proportion 0.1 --neg_proportion 1 --uniref_threshold 90 --score_key string_combined_score --score_threshold 950 \\
                --preloaded_protein_splits_path dataset.h5 --seed 8675309 --trim_unseen_proteins False --negatives_path None --taxon 9606 \\
                --weighted_random False --scramble_proteins False --exclude_preloaded_from_neg True

        :param pathlib.Path processed_folder: The processed folder, where you will deposit the new dataset.
        :param pathlib.Path common_path: The path to the common file.
        :param List[int] c_types: The different Park & Marcotte C-type levels to generate. Takes a list. *e.g.*: [1,2]
        :param float train_proportion: The proportion of interactions to assign to the training fold. Defaults to 0.8.
        :param float val_proportion: The proportion of interactions to assign to the validation fold. Defaults to 0.1.
        :param float test_proportion: The proportion of interactions to assign to the testing fold. Defaults to 0.1.
        :param float neg_proportion: The proportion of interactions that will be negative interactions. Defaults to 1.
        :param int uniref_threshold: The UniRef threshold rate to use to ensure proteins between splits are not too similar. Defaults to 90.
        :param Optional[str] score_key: The scoring key to theshold by, if any. Defaults to ``None``.
        :param Optional[str] score_threshold: The value to threshold the score key by. Values below this value will be filtered out. Defaults to ``None``.
        :param Optional[pathlib.Path] preloaded_protein_splits_path: Load protein splits from another RAPPPID dataset. Defaults to ``None``.
        :param int seed: An integer that will serve as the random seed for datasets. Defaults to 8675309.
        :param bool trim_unseen_proteins: If true, when a protein loaded from the preloaded_protein_splits is not found in the common dataset, it is not included in the dataset. Defaults to ``False``.
        :param Optional[pathlib.Path] negatives_path: Optional, path to file with negative interactions.  Defaults to ``None``.
        :param Optional[int] taxon: Optional, restrict a dataset to a certain organism.  Defaults to ``None``.
        :param bool weighted_random: If true, negative samples will be sampled in such a way as to maintain the same protein degree as the positive samples.  Defaults to ``False``.
        :param bool scramble_proteins: Scramble the association between protein ids and their sequences.  Defaults to ``False``.
        :param bool exclude_preloaded_from_neg: Set this to true if you don't want preloaded proteins to leak into negative.  Defaults to ``True``.
        :return: None
        """
        console = Console()

        trim_unseen_proteins = bool(trim_unseen_proteins)
        weighted_random = bool(weighted_random)
        scramble_proteins = bool(scramble_proteins)
        exclude_preloaded_from_neg = bool(exclude_preloaded_from_neg)

        rapppid.common_to_rapppid(
            console,
            Path(processed_folder),
            Path(common_path),
            c_types,
            train_proportion,
            val_proportion,
            test_proportion,
            neg_proportion,
            uniref_threshold,
            score_key,
            score_threshold,
            preloaded_protein_splits_path,
            seed,
            trim_unseen_proteins,
            negatives_path,
            taxon,
            weighted_random,
            scramble_proteins,
            exclude_preloaded_from_neg
        )

    @staticmethod
    def multimerge_rapppid(processed_folder: Path, dataset_paths: str) -> None:
        """
        Combine multiple RAPPPID datasets into one.

        It's important to note that if the datasets you are merging have different protein splits, then you will have
        data leakage. Only perform this operation if the two datasets have the same protein splits.

        You can call this function from the CLI using::
        
            ppi_origami process multimerge_rapppid PROCESSED_FOLDER DATASET_PATHS

        :param pathlib.Path processed_folder: The processed folder, where you will deposit the new dataset.
        :param str dataset_paths: Comma-seperated paths to the datasets
        :return: None
        """
        console = Console()
        processed_folder = Path(processed_folder)
        dataset_paths = [
            Path(dataset_path) for dataset_path in dataset_paths.split(',')
        ]

        print(dataset_paths)
        rapppid.multimerge_rapppid(console, processed_folder, dataset_paths)

    @staticmethod
    def train_sentencepiece_model(
            processed_folder: Path, dataset_file: Path, seed: int, vocab_size: int
    ) -> None:
        """
        Train a Uniword model using SentencePiece given a RAPPPID dataset.

        You can call this function from the CLI using::
        
            ppi_origami process train_sentencepiece_model PROCESSED_FOLDER DATASET_FILE SEED VOCAB_SIZE

        :param pathlib.Path processed_folder: The processed folder, where you will deposit the new dataset.
        :param pathlib.Path dataset_file: The RAPPPID file from which to train a Sentencepiece model.
        :param int seed: The random seed to use.
        :param int vocab_size: The number of tokens to learn.
        :return: None
        """
        rapppid.train_sentencepiece_model(
            processed_folder, dataset_file, seed, vocab_size
        )

    @staticmethod
    def dscript_upkb(raw_folder: Path, processed_folder: Path, taxon: Optional[int] = None) -> None:
        """
        Add a column of UPKB accessions to a D-SCRIPT-formatted datasets.

        You can call this function from the CLI using::
        
            ppi_origami process dscript_upkb RAW_FOLDER PROCESSED_FOLDER --taxon 9606

        :param pathlib.Path raw_folder: The folder datasets have been downloaded to.
        :param pathlib.Path processed_folder: The folder to output processed data.
        :param Optional[int] taxon: The NCBI taxon ID of the organism whose links you wish to download. Omit for all organisms. Defaults to ``None``.
        :return: None
        """
        console = Console()
        console.log('Running the "dscript_upkb" processing task.')
        dscript.process_dscript_upkb(
            Path(raw_folder), Path(processed_folder), taxon, console
        )

    @staticmethod
    def dscript_uniref(processed_folder: Path, threshold: int, taxon: Optional[int] = None) -> None:
        """
        Add a column of UniRef IDs to a D-SCRIPT-formatted datasets.

        You can call this function from the CLI using::
        
            ppi_origami process dscript_uniref PROCESSED_FOLDER THRESHOLD --taxon 9606

        :param pathlib.Path processed_folder: The folder to output processed data.
        :param int threshold: The UniRef identity threshold. Must be one of 50, 90, 100.
        :param Optional[int] taxon: The NCBI taxon ID of the organism whose links you wish to download. Omit for all organisms. Defaults to ``None``.
        :return: None
        """

        console = Console()
        console.log('Running the "dscript_uniref" processing task.')
        dscript.process_dscript_uniref(Path(processed_folder), taxon, threshold, console)

    @staticmethod
    def rapppid_to_dscript(
            rapppid_path: Path, dscript_folder: Path, c_types: List[int], trunc_len: int = 1500
    ) -> None:
        """
        Convert a RAPPPID HDF5 dataset to the D-SCRIPT dataset format.

        You can call this function from the CLI using::
        
            ppi_origami process rapppid_to_dscript RAPPPID_PATH DSCRIPT_FOLDER C_TYPES --trunc_len 1500

        :param pathlib.Path rapppid_path: Path to the RAPPPID dataset.
        :param pathlib.Path dscript_folder: Path to the folder into which we write the D-SCRIPT dataset.
        :param List[int] c_types: The Park & Marcotte C-type datasets to generate.
        :param int trunc_len: Length at which to truncate amino acid sequences.
        :return: None
        """
        for c_type in c_types:
            os.makedirs(Path(dscript_folder) / f"C{c_type}/", exist_ok=True)

            seq_path = Path(dscript_folder) / f"C{c_type}/seqs.fasta"
            train_path = Path(dscript_folder) / f"C{c_type}/train.tsv"
            test_path = Path(dscript_folder) / f"C{c_type}/test.tsv"
            val_path = Path(dscript_folder) / f"C{c_type}/val.tsv"
            metadata_path = Path(dscript_folder) / "metadata.json"

            rapppid.rapppid_to_dscript(
                rapppid_path,
                c_type,
                seq_path,
                train_path,
                test_path,
                val_path,
                metadata_path,
                trunc_len
            )

    @staticmethod
    def rapppid_to_sprint(rapppid_path: Path, sprint_folder: Path, c_types: List[int]) -> None:
        """
        Convert a RAPPPID HDF5 dataset to the SPRINT dataset format.

        You can call this function from the CLI using::
        
            ppi_origami process rapppid_to_sprint RAPPPID_PATH SPRINT_FOLDER C_TYPES

        :param pathlib.Path rapppid_path: Path to the RAPPPID dataset.
        :param pathlib.Path sprint_folder: Path to the folder into which we write the SPRINT dataset.
        :param List[int] c_types: The Park & Marcotte C-type datasets to generate.
        :return: None
        """
        for c_type in c_types:
            os.makedirs(Path(sprint_folder) / f"C{c_type}/", exist_ok=True)

            seq_path = Path(sprint_folder) / f"C{c_type}/seqs.fasta"
            train_pos_path = Path(sprint_folder) / f"C{c_type}/train.pos.txt"
            train_neg_path = Path(sprint_folder) / f"C{c_type}/train.neg.txt"
            train_seq_path = Path(sprint_folder) / f"C{c_type}/train.seq.txt"
            val_pos_path = Path(sprint_folder) / f"C{c_type}/val.pos.txt"
            val_neg_path = Path(sprint_folder) / f"C{c_type}/val.neg.txt"
            test_pos_path = Path(sprint_folder) / f"C{c_type}/test.pos.txt"
            test_neg_path = Path(sprint_folder) / f"C{c_type}/test.neg.txt"
            metadata_path = Path(sprint_folder) / "metadata.json"

            rapppid.rapppid_to_sprint(
                rapppid_path,
                c_type,
                seq_path,
                train_pos_path,
                train_neg_path,
                train_seq_path,
                val_pos_path,
                val_neg_path,
                test_pos_path,
                test_neg_path,
                metadata_path,
            )

    @staticmethod
    def rapppid_to_pipr(rapppid_path: Path, pipr_folder: Path, c_types: List[int]):
        """
        Convert a RAPPPID HDF5 dataset to the PIPR dataset format.

        You can call this function from the CLI using::
        
            ppi_origami process rapppid_to_pipr RAPPPID_PATH PIPR_FOLDER C_TYPES

        :param pathlib.Path rapppid_path: Path to the RAPPPID dataset.
        :param pathlib.Path pipr_folder: Path to the folder into which we write the PIPR dataset.
        :param List[int] c_types: The Park & Marcotte C-type datasets to generate.
        :return: None
        """

        for c_type in c_types:
            os.makedirs(Path(pipr_folder) / f"C{c_type}/", exist_ok=True)

            seq_path = Path(pipr_folder) / f"C{c_type}/seqs.tsv"
            train_path = Path(pipr_folder) / f"C{c_type}/train.tsv"
            test_path = Path(pipr_folder) / f"C{c_type}/test.tsv"
            val_path = Path(pipr_folder) / f"C{c_type}/val.tsv"
            metadata_path = Path(pipr_folder) / "metadata.json"

            pipr.rapppid_to_pipr(
                rapppid_path,
                c_type,
                seq_path,
                train_path,
                test_path,
                val_path,
                metadata_path,
            )

    @staticmethod
    def rapppid_to_deepppi(
            rapppid_path: Path, deepppi_folder: Path, c_types: List[int]
    ):
        """
        Convert a RAPPPID HDF5 dataset to the DeepPPI dataset format.

        You can call this function from the CLI using::
        
            ppi_origami process rapppid_to_deepppi RAPPPID_PATH DEEPPPI_FOLDER C_TYPES

        :param pathlib.Path rapppid_path: Path to the RAPPPID dataset.
        :param pathlib.Path deepppi_folder: Path to the folder into which we write the DeepPPI dataset.
        :param List[int] c_types: The Park & Marcotte C-type datasets to generate.
        :return: None
        """
        for c_type in c_types:
            os.makedirs(Path(deepppi_folder) / f"C{c_type}/", exist_ok=True)

            seq_path = Path(deepppi_folder) / f"C{c_type}/seqs.tsv"
            train_path = Path(deepppi_folder) / f"C{c_type}/train.tsv"
            test_path = Path(deepppi_folder) / f"C{c_type}/test.tsv"
            val_path = Path(deepppi_folder) / f"C{c_type}/val.tsv"
            metadata_path = Path(deepppi_folder) / "metadata.json"

            deepppi.rapppid_to_deepppi(
                rapppid_path,
                c_type,
                train_path,
                test_path,
                val_path,
                metadata_path,
            )

    @staticmethod
    def oma_upkb_groups(
            raw_path: Path, processed_path: Path, limit_taxons: Optional[List[int]] = None
    ):
        """
        Create a LevelDB database mapping UniProt accession codes to OMA Group IDs. To do this, the raw OMA XML file is parsed. This does, however mean, that you'll have had to run the command: ::

            ppi_origami download oma RAW_FOLDER

        You can call this function from the CLI using::
        
            ppi_origami process oma_upkb_groups RAW_PATH PROCESSED_PATH --limit_taxons [9606,10090]

        :param pathlib.Path raw_path: The folder datasets have been downloaded to.
        :param pathlib.Path processed_path: The folder to output processed data.
        :param Optional[List[int]] limit_taxons: A list of NCBI Taxon IDs. Will only parse IDs that belong to these taxa. Defaults to `None` (*i.e.*: No restriction on taxa).
        :return: None
        """


        oma.upkb_groups(Path(raw_path), Path(processed_path), limit_taxons)

    @staticmethod
    def rapppid_to_intrepppid(
        processed_path: Path,
        rapppid_path: Path,
        intrepppid_path: Path,
        c_types: List[int],
        allowlist_taxon: Optional[List[int]] = None,
        denylist_taxon: Optional[List[int]] = None,
        scramble_interactions: bool = False,
        scramble_orthologs: bool = False,
        uniref_threshold: int = 90
    ):
        """
        Convert a RAPPPID dataset to the INTREPPPID format. This primarily involves adding orthology data to the dataset.

        You can call this function from the CLI using::
        
            ppi_origami process rapppid_to_intrepppid PROCESSED_PATH RAPPPID_PATH INTREPPPID_PATH C_TYPES --allowlist_taxon [9606,10090] --denylist_taxon None --scramble_interactions False --scramble_orthologs False --uniref_threshold 90

        :param pathlib.Path processed_path: The folder where processed data is kept.
        :param pathlib.Path rapppid_path: The path to the RAPPPID dataset to convert to INTREPPPID.
        :param pathlib.Path intrepppid_path: The path to save the INTREPPPID dataset.
        :param List[int] c_types: The different Park & Marcotte C-type levels to generate. Takes a list. *e.g.*: [1,2]
        :param Optional[List[int]] allowlist_taxon: The NCBI Taxon IDs of organism for which orthologues are allowed to be from. Orthologues from orgnaism not in this list will be omitted. Cannot to be used with denylist. Defaults to ``None``.
        :param Optional[List[int]] denylist_taxon: The NCBI Taxon IDs of organism for which orthologues are not allowed to be from. Orthologues from orgnaism in this list will be omitted. Cannot to be used with denylist. Defaults to ``None``.
        :param bool scramble_interactions: If True, protein IDs will be scrambled, ablating the biological meaning of the interaction network. Defaults to ``False``.
        :param bool scramble_orthologs: If True, orthologue IDs will be scrambled, ablating the biological meaning of the orthologue data. Defaults to ``False``.
        :param int uniref_threshold: What uniref threshold to use when determining the similarity of proteins. Must be one of 50, 90, or 100. Defaults to 90.
        :return: None
        """

        processed_path = Path(processed_path)
        rapppid_path = Path(rapppid_path)
        intrepppid_path = Path(intrepppid_path)
        intrepppid.rapppid_to_intrepppid(processed_path, rapppid_path, intrepppid_path, c_types, allowlist_taxon,
                                         denylist_taxon, scramble_interactions, scramble_orthologs, uniref_threshold)

    def __rapppid_inflate_eval_negatives(self, processed_path: Path, rapppid_path: Path, negative_ratio: float, c_types: List[int]):
        processed_path = Path(processed_path)
        rapppid.inflate_eval_negatives(processed_path, rapppid_path, negative_ratio, c_types)


class Download(object):
    @staticmethod
    def uniref(raw_folder: Path, threshold: int):
        """
        Download the `UniRef <https://www.uniprot.org/uniref>`_ dataset from `UniProt <https://www.uniprot.org/>`_. You must specify a similarity threshold among the three available options: 50%, 90%, and 100%.

        You can call this function from the CLI using::
        
            ppi_origami download uniref RAW_FOLDER THRESHOLD

        :param pathlib.Path raw_folder: The folder to download the dataset to.
        :param int threshold: The UniRef identity threshold. Must be one of 50, 90, 100.
        :return: None
        """
        uniref.download(Path(raw_folder), threshold)

    @staticmethod
    def uniprot_sec_ac(raw_folder: Path):
        """
        Download UniProt secondary accessions to the raw folder.

        You can call this function from the CLI using::
        
            ppi_origami download uniprot_sec_ac RAW_FOLDER

        :param pathlib.Path raw_folder: The folder to download the dataset to.
        :return: None
        """
        uniprot.download("sec_ac", None, Path(raw_folder))

    @staticmethod
    def uniprot_delac(raw_folder: Path):
        """
        Download `UniProt <https://www.uniprot.org/>`_ deleted accessions to the raw folder. More info can be found in the `UniProtKB Manual <https://www.uniprot.org/help/deleted_accessions>`_.

        You can call this function from the CLI using::
        
            ppi_origami download uniprot_delac RAW_FOLDER

        :param pathlib.Path raw_folder: The folder to download the dataset to.
        :return: None
        """
        uniprot.download("delac", None, Path(raw_folder))

    @staticmethod
    def uniprot_id_mapping(raw_folder: Path):
        """
        Download `UniProt <https://www.uniprot.org/>`_ ID mappings to the raw folder. More info can be found on `UniProt.org <https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README>`_.

        You can call this function from the CLI using::
        
            ppi_origami download uniprot_id_mapping RAW_FOLDER

        :param pathlib.Path raw_folder: The folder to download the dataset to.
        :return: None
        """
        uniprot.download("id_mapping", None, Path(raw_folder))

    @staticmethod
    def uniprot_seqs_db(processed_folder: Path, taxon: Optional[int] = None):
        """
        Download `UniProt <https://www.uniprot.org/>`_ sequences and saves them to a `LevelDB <https://en.wikipedia.org/wiki/LevelDB>`_ file in the processed folder. Will download sequences for the specified taxon if taxon is not ``None``.

        You can call this function from the CLI using::
        
            ppi_origami download uniprot_seqs_db PROCESSED_FOLDER --taxon 9606

        :param pathlib.Path processed_folder: The folder to save the database to.
        :param int taxon: The NCBI taxon ID of the organism whose links you wish to download. Omit to download sequences for all organisms. Defaults to ``None``.
        :return: None
        """
        uniprot.make_seq_db(processed_folder=Path(processed_folder), taxon=taxon)

    @staticmethod
    def uniprot_seqs_fasta(raw_folder: Path):
        """
        Download `UniProt <https://www.uniprot.org/>`_ sequences in the FASTA format to the raw folder.

        You can call this function from the CLI using::
        
            ppi_origami download uniprot_seqs_fasta RAW_FOLDER

        :param pathlib.Path raw_folder: The folder to download the dataset to.
        :return: None
        """
        uniprot.download("uniprot_seq_fasta", None, Path(raw_folder))

    @staticmethod
    def string_links(raw_folder: Path, version: str = "12.0", taxon: Optional[int] = None):
        """
        Download the `STRING <https://string-db.org/>`_ links dataset. Downloads the ``protein.links.{version}.txt.gz`` file.

        You can call this function from the CLI using::
        
            ppi_origami download string_links RAW_FOLDER --version 12.0 --taxon 9606

        :param pathlib.Path raw_folder: The folder to download the dataset to.
        :param str version: The version of the STRING database to download.
        :param Optional[int] taxon: The NCBI taxon ID of the organism whose links you wish to download. Omit for all organisms. Defaults to ``None``.
        :return: None
        """
        string.download("links", taxon, Path(raw_folder), version)

    @staticmethod
    def string_links_detailed(raw_folder: Path, version: str = "12.0", taxon: Optional[int] = None):
        """
        Download the `STRING <https://string-db.org/>`_ links dataset. Downloads the ``protein.links.detailed.{version}.txt.gz`` file.

        You can call this function from the CLI using::
        
            ppi_origami download string_links_detailed RAW_FOLDER --version 12.0 --taxon 9606

        :param pathlib.Path raw_folder: The folder to download the dataset to.
        :param str version: The version of the STRING database to download.
        :param Optional[int] taxon: The NCBI taxon ID of the organism whose links you wish to download. Omit for all organisms. Defaults to ``None``.
        :return: None
        """
        string.download("links_detailed", taxon, Path(raw_folder), version)

    @staticmethod
    def string_physical_links(raw_folder: Path, version: str = "12.0", taxon: Optional[int] = None):
        """
        Download the `STRING <https://string-db.org/>`_ physical links dataset. Downloads the ``protein.physical.links.{version}.txt.gz`` file. This dataset only includes PPIs with evidence of physical interactions.

        You can call this function from the CLI using::
        
            ppi_origami download string_physical_links RAW_FOLDER --version 12.0 --taxon 9606

        :param pathlib.Path raw_folder: The folder to download the dataset to.
        :param str version: The version of the STRING database to download.
        :param Optional[int] taxon: The NCBI taxon ID of the organism whose links you wish to download. Omit for all organisms.
        :return: None
        """
        string.download("physical_links", taxon, Path(raw_folder), version)

    @staticmethod
    def string_physical_detailed_links(raw_folder: Path, version: str = "12.0", taxon: Optional[int] = None):
        """
        Download the `STRING <https://string-db.org/>`_ physical links dataset. Downloads the ``protein.physical.links.detailed.{version}.txt.gz`` file. This dataset only includes PPIs with evidence of physical interactions, and provide more information than ``string_physical_links``.

        You can call this function from the CLI using::
        
            ppi_origami download string_physical_detailed_links RAW_FOLDER --version 12.0 --taxon 9606

        :param pathlib.Path raw_folder: The folder to download the dataset to.
        :param str version: The version of the STRING database to download.
        :param Optional[int] taxon: The NCBI taxon ID of the organism whose links you wish to download. Omit for all organisms.
        :return: None
        """
        string.download("physical_detailed_links", taxon, Path(raw_folder), version)

    @staticmethod
    def string_aliases(raw_folder: Path, version: str = "12.0"):
        """
        Download the `STRING <https://string-db.org/>`_ aliases dataset. These are mappings between STRING identifiers, and other identifiers (most importantly, UniProt).

        You can call this function from the CLI using::
        
            ppi_origami download string_aliases RAW_FOLDER --version 12.0

        :param pathlib.Path raw_folder: The folder to download the dataset to.
        :param str version: The version of the STRING database to download.
        :return: None
        """
        string.download("aliases", None, Path(raw_folder), version)

    @staticmethod
    def biogrid(raw_folder: Path, version: str = "4.4.224"):
        """
        Download the `BioGRID <https://thebiogrid.org/>`_ PPI dataset to the raw folder. Namely, it downloads the ``BIOGRID-ORGANISM-{version}.mitab.zip`` file from the `BioGRID release archive <https://downloads.thebiogrid.org/BioGRID/Release-Archive/>`_.

        You can call this function from the CLI using::
        
            ppi_origami download biogrid RAW_FOLDER --version 4.4.224

        :param pathlib.Path raw_folder: The folder to download the dataset to.
        :param str version: The BioGRID version to download, defaults to "4.4.224".
        :return: None
        """
        biogrid.download(Path(raw_folder), version)

    @staticmethod
    def hippie(raw_folder: Path, version: str = "current"):
        """
        Download the `HIPPIE <https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/>`_ PPI dataset to the raw folder. You may specify the version to download. A value of "current" results in the latest version being downloaded. To download older version, supply a version number like "2.2".

        You can call this function from the CLI using::
        
            ppi_origami download hippie RAW_FOLDER --version current

        :param pathlib.Path raw_folder: The folder to download the dataset to.
        :param str version: The version of the HIPPIE dataset to download, defaults to "current".
        :return: None
        """
        hippie.download(Path(raw_folder))

    @staticmethod
    def dscript(raw_folder: Path, taxon: int):
        """
        Download a `D-SCRIPT <https://dscript.csail.mit.edu/>`_ PPI dataset to the raw folder. D-SCRIPT specifies datasets for *H. sapiens*, *M. musculus*, *D. melanogaster*, *S. cerevisiae*, *C. elegans*, and *E. coli*.

        These organisms correspond to the `NCBI Taxon IDs <https://www.ncbi.nlm.nih.gov/taxonomy>`_: 9606, 10090, 7227, 4932, 6239, and 511145. These are the only valid values of the ``taxon`` argument.

        You can call this function from the CLI using::
        
            ppi_origami download dscript RAW_FOLDER TAXON

        :param pathlib.Path raw_folder: The folder to download the dataset to.
        :param int taxon: The NCBI taxon ID of the organism whose links you wish to download. Must be on of 9606, 10090, 7227, 4932, 6239, or 511145.
        :raises ValueError: A ValueError is raised when ``taxon`` is not one of 9606, 10090, 7227, 4932, 6239, or 511145.
        :return: None
        """
        dscript.download(Path(raw_folder), taxon)

    @staticmethod
    def oma(raw_folder: Path):
        """
        Download orthology data from `OMA <https://omabrowser.org/oma/home/>`_ to the raw folder. Specifically, it downloads the orthology data in a gzipped XML format (``oma-groups.orthoXML.xml.gz``), as well as mappings between OMA identifiers and UniProt identifiers (``oma-uniprot.txt.gz``).

        You can call this function from the CLI using::
        
            ppi_origami download oma RAW_FOLDER

        :param pathlib.Path raw_folder: The folder to download the dataset to.
        :return: None
        """
        oma.download(Path(raw_folder))


class Analysis(object):

    @staticmethod
    def verify_rapppid(dataset_path: Path, taxon_ids: Optional[List[int]] = None, sample_fraction: Optional[float] = None,
                       skip_protein_overlap: bool = False, skip_gotoh: bool = False, skip_sw: bool = False,
                       skip_taxa_check: bool = False
                       ):
        """
        Runs some sanity checks on a RAPPPID dataset. Specifically it can test:
        
        1. Overlap between the protein splits
        2. Overlap between proteins observed in interaction pairs belonging to different splits.
        3. Checks whether proteins belong to a set of expected taxa
        4. Checks whether sequences between splits are similar using either/both Gotoh or Smith-Waterman algorithms.

        You can call this function from the CLI using::
        
            ppi_origami analysis verify_rapppid DATASET_PATH --taxon_ids [9606,10090] --sample_fraction 0.05 --skip_protein_overlap False --skip_gotoh False --skip_sw False --skip_taxa_check False

        :param pathlib.Path dataset_path: The path to the RAPPPID dataset.
        :param Optional[List[int]] taxon_ids: A list of taxa which proteins will be required to be from. If proteins are detected from taxa other than those specified, the test will fail. Defaults to ``None``.
        :param Optional[float] sample_fraction: The fraction of dataset samples to test for sequence similarity and verifying protein organisms. Dataset samples are randomly chosen. Set to 1 or `None` to test the whole dataset. Defaults to ``None``.
        :param bool skip_protein_overlap: Skip the protein overlap test.
        :param bool skip_gotoh: Skip checking the sequence similarity with the Gotoh algorithm.
        :param bool skip_sw: Skip checking the sequence similarity with the Smith-Waterman algorithm.
        :param bool skip_taxa_check: Skip checking the protein's organism.
        :return: None
        """

        console = Console()

        # This stores all the data that will be saved in the ODS (Excel) spreadsheet.
        ods_data = OrderedDict()

        if skip_protein_overlap is False:
            ods_data.update({"Int'n Protein Overlap": []})
            ods_data.update({"Protein Split Overlap": []})
            ods_data.update({"Train Test Overlap": []})
            ods_data.update({"Train Val Overlap": []})
            ods_data.update({"Val Test Overlap": []})

        if skip_gotoh is False:
            ods_data.update({"Gotoh Sequence Similarity (Hist)": []})
            ods_data.update({"Gotoh Sequence Similarity": []})

        if skip_sw is False:
            ods_data.update({"SW Sequence Similarity (Hist)": []})
            ods_data.update({"SW Sequence Similarity": []})

        if skip_taxa_check is False:
            ods_data.update({"UPKB Taxa": []})

        if taxon_ids:
            taxon_ids = set(taxon_ids)

        observed_proteins = set()
        observed_proteins_split = {c: {split: set() for split in ["train", "val", "test"]} for c in [1, 2, 3]}

        rapppid_dataset = tb.open_file(str(dataset_path))

        interaction_tables = {
            3: {
                "train": rapppid_dataset.root.interactions.c3.c3_train,
                "test": rapppid_dataset.root.interactions.c3.c3_test,
                "val": rapppid_dataset.root.interactions.c3.c3_val
            },
            2: {
                "train": rapppid_dataset.root.interactions.c2.c2_train,
                "test": rapppid_dataset.root.interactions.c2.c2_test,
                "val": rapppid_dataset.root.interactions.c2.c2_val,
            },
            1: {
                "train": rapppid_dataset.root.interactions.c1.c1_train,
                "test": rapppid_dataset.root.interactions.c1.c1_test,
                "val": rapppid_dataset.root.interactions.c1.c1_val,
            }
        }

        split_tables = {
            "train": rapppid_dataset.root.splits.train,
            "test": rapppid_dataset.root.splits.test,
            "val": rapppid_dataset.root.splits.val,
        }

        if skip_protein_overlap is False:
            for c in [1, 2, 3]:
                for split in ["train", "test", "val"]:
                    for row in interaction_tables[c][split]:
                        observed_proteins.add(row['protein_id1'].decode('utf8'))
                        observed_proteins.add(row['protein_id2'].decode('utf8'))
                        observed_proteins_split[c][split].add(row['protein_id1'].decode('utf8'))
                        observed_proteins_split[c][split].add(row['protein_id2'].decode('utf8'))

            sequences = {}

            for row in rapppid_dataset.root.sequences:
                sequences[row[0].decode('utf8')] = row[1].decode('utf8')

            verify_pass, observed_overlap = utils.verify_observed_splits(observed_proteins_split)

            console.log("🧪 Verifying Interaction Protein Overlap...")
            if verify_pass:
                console.log("[green][b]PASS")
            else:
                console.log("[red][b]FAIL")

            table = Table()

            table.add_column("C")
            table.add_column("Split Pair")
            table.add_column("# Overlapping Proteins", justify='right')

            colors = ['', 'magenta', 'yellow', 'cyan']

            ods_data["Int'n Protein Overlap"].append(["C", "Split Pair", "# Overlapping Proteins"])
            for c in [1, 2, 3]:
                for i, split_pair in enumerate(['train/val', 'train/test', 'test/val']):
                    display_c = f"C{c}" if i == 0 else ""
                    display_overlap = humanize.intcomma(len(observed_overlap[c][split_pair]))
                    table.add_row(display_c, split_pair, display_overlap, style=colors[c])
                    ods_data["Int'n Protein Overlap"].append([display_c, split_pair, display_overlap])

                    if split_pair == "train/val":
                        sheet = "Train Val Overlap"
                    elif split_pair == "train/test":
                        sheet = "Train Test Overlap"
                    else:
                        sheet = "Val Test Overlap"

                    for protein in observed_overlap[c][split_pair]:
                        ods_data[sheet].append([f"C{c}", protein])

            console.log(table)

            splits = {split: set() for split in ['train', 'test', 'val']}

            for split in ['train', 'test', 'val']:
                for row in split_tables[split]:
                    try:
                        splits[split].add(row.decode('utf8'))
                    except AttributeError:
                        console.log(row)
                        splits[split].add(row)

            console.log("🧪 Verifying Protein Split Overlap...")

            split_overlap = {
                "train/val": splits['train'].intersection(splits['val']),
                "train/test": splits['train'].intersection(splits['test']),
                "test/val": splits['test'].intersection(splits['val'])
            }

            if len(split_overlap['train/val']) == 0 and len(split_overlap['train/test']) == 0 and \
                len(split_overlap['test/val']) == 0:
                console.log("[green][b]PASS")
            else:
                console.log("[red][b]FAIL")

            table = Table()

            table.add_column("Split Pair")
            table.add_column("# Overlapping Proteins", justify='right')

            ods_data["Protein Split Overlap"].append(["Split Pair", "# Overlapping Proteins"])
            for split_pair in ['train/val', 'train/test', 'test/val']:
                display_overlap = humanize.intcomma(len(split_overlap[split_pair]))
                table.add_row(split_pair, display_overlap)
                ods_data["Protein Split Overlap"].append([split_pair, display_overlap])

            console.log(table)

        if skip_gotoh is False or skip_sw is False:
            console.log("🧪 Verifying cross-split similarity of dataset sequences...")

            distances = {}

            if skip_gotoh is False:
                distances.update({'gotoh': []})
                ods_data["Gotoh Sequence Similarity"].append(["Gotoh Distance",
                                                              "Gotoh Distance (Normalized)"
                                                              "Train Protein Length",
                                                              "Test Protein Length",
                                                              "Train UPKB AC",
                                                              "Test UPKB AC",
                                                              "Train Protein Sequence",
                                                              "Test Protein Sequence"])

            if skip_sw is False:
                distances.update({'smith_waterman': []})

                ods_data["SW Sequence Similarity"].append(["Smith-Waterman Distance",
                                                          "Smith-Waterman (Normalized)"
                                                          "Train Protein Length",
                                                          "Test Protein Length",
                                                          "Train UPKB AC",
                                                          "Test UPKB AC",
                                                          "Train Protein Sequence",
                                                          "Test Protein Sequence"])

            if sample_fraction is not None:
                splits['train'] = sample(list(splits['train']), int(len(splits['train'])*sample_fraction))
                splits['test'] = sample(list(splits['test']), int(len(splits['test'])*sample_fraction))

            with Progress() as progress:
                train_proteins_pbar = progress.add_task("[red]Train Proteins...", total=len(splits['train']))

                total = 0
                for train_protein_id in list(splits['train']):
                    train_protein = sequences[train_protein_id]

                    test_proteins_pbar = progress.add_task("[green]Test Proteins...", total=len(splits['test']))
                    for test_protein_id in list(splits['test']):

                        test_protein = sequences[test_protein_id]

                        max_len = max(len(train_protein), len(test_protein))

                        if skip_gotoh is False:
                            gotoh_distance = float(gotoh(train_protein, test_protein))
                            distances['gotoh'].append(gotoh_distance / max_len)
                            ods_data["Gotoh Sequence Similarity"].append([gotoh_distance,
                                                                          gotoh_distance / max_len,
                                                                          len(train_protein),
                                                                          len(test_protein),
                                                                          train_protein_id,
                                                                          test_protein_id,
                                                                          train_protein,
                                                                          test_protein])

                        if skip_sw is False:
                            sw_distance = float(smith_waterman(train_protein, test_protein))
                            distances['smith_waterman'].append(sw_distance / max_len)
                            ods_data["SW Sequence Similarity"].append([sw_distance,
                                                                      sw_distance / max_len,
                                                                      len(train_protein),
                                                                      len(test_protein),
                                                                      train_protein_id,
                                                                      test_protein_id,
                                                                      train_protein,
                                                                      test_protein])
                            
                        progress.update(test_proteins_pbar, advance=1)
                        total += 1
                    progress.update(train_proteins_pbar, advance=1)
                    progress.update(test_proteins_pbar, visible=False)

            lower_bins = [0.1 * i for i in range(11)]
            hists = {key: np.histogram(distances[key], bins=lower_bins) for key in distances}
            upper_bins = [0.1 * (i+1) for i in range(11)]

            total_pairs = len(splits['train'])*len(splits['test'])

            for key in hists:
                table = Table(title=f'Train/Test Alignment Score ({key})')
                table.add_column("Bin")
                table.add_column("Counts", justify='right')
                table.add_column("%", justify='right')

                counts = hists[key][0]

                if key == "gotoh":
                    spreadsheet_name = "Gotoh Sequence Similarity (Hist)"
                else:
                    spreadsheet_name = "SW Sequence Similarity (Hist)"

                for i, count in enumerate(counts):
                    bin = f"{lower_bins[i]:.2}-{upper_bins[i]:.2}"
                    perc = 100*count / total_pairs
                    table.add_row(bin, humanize.intcomma(count), f"{perc:.3}%")
                    ods_data[spreadsheet_name].append([bin, humanize.intcomma(count), f"{perc:.3}%"])

                console.log(table)

        if skip_taxa_check is False:
            console.log("🧪 Verifying proteins are from expected organisms...")
            console.log("Collecting protein information from UniProt.org")
            if sample_fraction is not None:
                protein_sample = list(observed_proteins)
                protein_sample = set(sample(protein_sample, int(len(protein_sample)*sample_fraction)))
                upkb_data = analysis.get_uniprot_data(protein_sample)
            else:
                upkb_data = analysis.get_uniprot_data(observed_proteins)

            ods_data["UPKB Taxa"].append(['UPKB AC', 'NCBI Taxon ID'])
            if taxon_ids:
                console.log("Verifying Protein Taxa...")
                taxon_pass, observed_taxa, upkb_taxa = utils.verify_taxon(upkb_data, taxon_ids)
                if taxon_pass:
                    console.log("[green][b]PASS")
                else:
                    console.log("[red][b]FAIL")
                    console.log("Observed the following taxa:")
                    console.log(observed_taxa)

                for upkb_ac in upkb_taxa:
                    ods_data["UPKB Taxa"].append([upkb_ac, upkb_taxa[upkb_ac]])

        save_ods(f"report_{Path(dataset_path).stem}.ods", ods_data)


class Wizard(object):

    @staticmethod
    def build_rapppid_dataset(raw_path: Path, processed_path: Path, taxon: int, string_version: str, uniref_threshold: int) -> None:
        console = Console()

        console.log(f"Running ppi_origami download string_physical_detailed_links {raw_path} --taxon {taxon}")
        Download.string_physical_detailed_links(raw_path, string_version, taxon)

        console.log(f"Running ppi_origami download string_aliases {raw_path} {string_version}")
        Download.string_aliases(raw_path, string_version)

        console.log(f"Running ppi_origami download uniref {raw_path} {uniref_threshold}")
        Download.uniref(raw_path, uniref_threshold)

        console.log(f"Running ppi_origami process uniref {raw_path} {processed_path} {uniref_threshold}")
        Download.uniref(raw_path, uniref_threshold)

        console.log(f"Running ppi_origami download uniprot_seqs {processed_path} {taxon}")
        Download.uniprot_seqs(raw_path, taxon)

        console.log(f"Running ppi_origami download uniprot_id_mapping {raw_path}")
        Download.uniprot_id_mapping(raw_path)

        console.log(f"Running ppi_origami download uniprot_delac {raw_path}")
        Download.uniprot_delac(raw_path)

        console.log(f"Running ppi_origami download uniprot_sec_ac {raw_path}")
        Download.uniprot_sec_ac(raw_path)

        console.log(f"Running ppi_origami process string_upkb {raw_path} {processed_path} {string_version} {taxon}")
        Process.string_upkb(raw_path, processed_path, string_version, taxon)

        console.log(f"Running ppi_origami process common_format {processed_path} string {taxon} uniref{uniref_threshold} {string_version}")
        Process.common_format(processed_path, "string", taxon, "upkb", string_version)


class Pipeline(object):
    def __init__(self):
        self.process = Process()
        self.download = Download()
        self.analysis = Analysis()
        self.wizard = Wizard()

def main():
    fire.Fire(Pipeline, name="ppi_origami")

if __name__ == "__main__":
    main()
